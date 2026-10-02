//! The Opus encoder (RFC 6716 §5): mode, bandwidth and frame-size choice,
//! rate control, and packet assembly over the SILK and CELT encoders.

use crate::celt::encoder::{CeltEncoder, FrameConfig};
use crate::error::{Result, config};
use crate::packet::{self, Bandwidth, Mode, Toc};
use crate::range::RangeEncoder;
use crate::resample::Resampler;
use crate::silk::encoder::SilkEncoder;

/// The total algorithmic delay of the encoder in samples at 48 kHz: the
/// `pre_skip` a stream should carry (RFC 7845 §4.2). 2.5 ms of CELT overlap
/// plus 4 ms of lookahead, which also absorbs the SILK path's resampling.
pub const LOOKAHEAD_48K: usize = 312;

/// Lookahead the SILK analysis gets beyond each frame, in ms.
const SILK_LOOKAHEAD_MS: f64 = 4.5;

/// What the encoder is tuned for (RFC 6716 §2.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Application {
    /// Music and general audio: CELT except at very low rates.
    #[default]
    Audio,
    /// Speech: SILK and hybrid modes at the rates where they do best.
    Voip,
    /// CELT only, for the lowest delay.
    LowDelay,
}

/// Encoder settings.
#[derive(Debug, Clone, Copy)]
pub struct EncoderConfig {
    /// Input sample rate: 8, 12, 16, 24 or 48 kHz.
    pub sample_rate: u32,
    /// 1 or 2.
    pub channels: usize,
    pub application: Application,
    /// Target bit rate in bit/s, 6 000 to 510 000 (per stream).
    pub bitrate: u32,
    /// Variable bit rate (frame sizes follow the signal) or constant.
    pub vbr: bool,
    /// 0 (fastest) to 10.
    pub complexity: u8,
    /// Packet duration in samples at 48 kHz: 120, 240, 480, 960, 1920 or
    /// 2880 (2.5 to 60 ms).
    pub frame_size: usize,
    /// Force a mode instead of choosing by rate and application.
    pub mode: Option<Mode>,
    /// The widest bandwidth to code (also capped by the input rate).
    pub max_bandwidth: Option<Bandwidth>,
    /// Code LBRR (in-band FEC) data in SILK and hybrid frames.
    pub fec: bool,
    /// Expected packet loss, percent (sizes the LBRR frames).
    pub packet_loss_percent: u8,
}

impl Default for EncoderConfig {
    fn default() -> Self {
        Self {
            sample_rate: 48000,
            channels: 2,
            application: Application::Audio,
            bitrate: 96_000,
            vbr: true,
            complexity: 10,
            frame_size: 960,
            mode: None,
            max_bandwidth: None,
            fec: false,
            packet_loss_percent: 0,
        }
    }
}

/// A delay line of interleaved samples.
struct Delay {
    buf: Vec<f32>,
}

impl Delay {
    fn new(samples: usize) -> Self {
        Self { buf: vec![0.0; samples] }
    }
    /// Pushes `input` and returns as many samples, delayed.
    fn process(&mut self, input: &[f32]) -> Vec<f32> {
        self.buf.extend_from_slice(input);
        self.buf.drain(..input.len()).collect()
    }
}

/// An Opus encoder for one stream (mono or stereo).
pub struct Encoder {
    cfg: EncoderConfig,
    upsamplers: Vec<Resampler>,
    upsample_delay: usize,
    celt: CeltEncoder,
    silk: SilkEncoder,
    celt_delay: Delay,
    silk_down: Vec<Resampler>,
    silk_buf: Vec<f32>,
    silk_fs: usize,
    prev_mode: Option<Mode>,
    final_range: u32,
    /// VBR reservoir in bits: positive when earlier frames used less than
    /// their share.
    reservoir: i64,
}

fn valid_frame_size(n: usize) -> bool {
    matches!(n, 120 | 240 | 480 | 960 | 1920 | 2880)
}

impl Encoder {
    pub fn new(cfg: EncoderConfig) -> Result<Self> {
        if !crate::decoder::SAMPLE_RATES.contains(&cfg.sample_rate) {
            return Err(config(format!("input rate {} is not 8, 12, 16, 24 or 48 kHz", cfg.sample_rate)));
        }
        if !(1..=2).contains(&cfg.channels) {
            return Err(config(format!("{} channels: one Opus stream has 1 or 2", cfg.channels)));
        }
        if !valid_frame_size(cfg.frame_size) {
            return Err(config(format!("frame size {} is not 2.5, 5, 10, 20, 40 or 60 ms", cfg.frame_size)));
        }
        if !(6000..=510_000).contains(&cfg.bitrate) {
            return Err(config(format!("bit rate {} outside 6-510 kb/s", cfg.bitrate)));
        }
        let c = cfg.channels;
        let (upsamplers, upsample_delay) = if cfg.sample_rate == 48000 {
            (Vec::new(), 0)
        } else {
            // A linear-phase interpolator whose delay is a whole number of
            // 48 kHz samples, so the stream's pre-skip stays exact.
            let d = 48;
            let ms = d as f64 / 48.0;
            ((0..c).map(|_| Resampler::new(cfg.sample_rate as usize, 48000, ms)).collect(), d)
        };
        Ok(Self {
            cfg,
            upsamplers,
            upsample_delay,
            celt: CeltEncoder::new(c),
            silk: SilkEncoder::new(c),
            celt_delay: Delay::new((LOOKAHEAD_48K - crate::celt::tables::OVERLAP) * c),
            silk_down: Vec::new(),
            silk_buf: Vec::new(),
            silk_fs: 0,
            prev_mode: None,
            final_range: 0,
            reservoir: 0,
        })
    }

    /// The settings in use.
    pub fn config(&self) -> &EncoderConfig {
        &self.cfg
    }

    /// Samples at 48 kHz by which the decoded output lags the input: the
    /// stream's pre-skip.
    pub fn lookahead(&self) -> usize {
        LOOKAHEAD_48K + self.upsample_delay
    }

    /// Input samples per channel that make one packet.
    pub fn frame_samples(&self) -> usize {
        self.cfg.frame_size * self.cfg.sample_rate as usize / 48000
    }

    /// The range coder's final state of the last packet (of its last frame).
    pub fn final_range(&self) -> u32 {
        self.final_range
    }

    /// Changes the target bit rate.
    pub fn set_bitrate(&mut self, bitrate: u32) -> Result<()> {
        if !(6000..=510_000).contains(&bitrate) {
            return Err(config(format!("bit rate {bitrate} outside 6-510 kb/s")));
        }
        self.cfg.bitrate = bitrate;
        Ok(())
    }

    /// Switches between variable and constant bit rate.
    pub fn set_vbr(&mut self, vbr: bool) {
        self.cfg.vbr = vbr;
    }

    /// Sets the complexity (0–10).
    pub fn set_complexity(&mut self, complexity: u8) {
        self.cfg.complexity = complexity.min(10);
    }

    /// The mode for the current settings.
    fn choose_mode(&self) -> Mode {
        if let Some(m) = self.cfg.mode {
            return match m {
                Mode::Silk | Mode::Hybrid if self.cfg.frame_size < 480 => Mode::Celt,
                Mode::Hybrid if self.max_bw() < Bandwidth::SuperWide => Mode::Silk,
                m => m,
            };
        }
        let per_ch = self.cfg.bitrate / self.cfg.channels as u32;
        if self.cfg.frame_size < 480 || self.cfg.application == Application::LowDelay {
            return Mode::Celt;
        }
        match self.cfg.application {
            Application::Voip => {
                if per_ch < 20_000 || self.max_bw() <= Bandwidth::Wide {
                    if self.max_bw() <= Bandwidth::Wide || per_ch < 12_000 { Mode::Silk } else { Mode::Hybrid }
                } else if per_ch < 36_000 {
                    Mode::Hybrid
                } else {
                    Mode::Celt
                }
            }
            _ => {
                if per_ch < 12_000 {
                    Mode::Silk
                } else {
                    Mode::Celt
                }
            }
        }
    }

    fn max_bw(&self) -> Bandwidth {
        let by_rate = match self.cfg.sample_rate {
            8000 => Bandwidth::Narrow,
            12000 => Bandwidth::Medium,
            16000 => Bandwidth::Wide,
            24000 => Bandwidth::SuperWide,
            _ => Bandwidth::Full,
        };
        self.cfg.max_bandwidth.map_or(by_rate, |b| b.min(by_rate))
    }

    fn choose_bandwidth(&self, mode: Mode) -> Bandwidth {
        let per_ch = self.cfg.bitrate / self.cfg.channels as u32;
        let max = self.max_bw();
        let want = match mode {
            Mode::Silk => {
                if per_ch < 9000 {
                    Bandwidth::Narrow
                } else if per_ch < 12000 {
                    Bandwidth::Medium
                } else {
                    Bandwidth::Wide
                }
            }
            Mode::Hybrid => {
                if per_ch < 28_000 {
                    Bandwidth::SuperWide
                } else {
                    Bandwidth::Full
                }
            }
            Mode::Celt => {
                if per_ch >= 20_000 {
                    Bandwidth::Full
                } else if per_ch >= 14_000 {
                    Bandwidth::SuperWide
                } else if per_ch >= 10_000 {
                    Bandwidth::Wide
                } else {
                    Bandwidth::Narrow
                }
            }
        };
        let bw = want.min(max);
        match (mode, bw) {
            (Mode::Celt, Bandwidth::Medium) => Bandwidth::Wide,
            (Mode::Hybrid, b) if b < Bandwidth::SuperWide => Bandwidth::SuperWide,
            (Mode::Silk, b) if b > Bandwidth::Wide => Bandwidth::Wide,
            (_, b) => b,
        }
    }

    /// Encodes one packet from `pcm`: exactly [`Self::frame_samples`]
    /// interleaved samples per channel at the input rate, ±1.0 full scale.
    pub fn encode(&mut self, pcm: &[f32]) -> Result<Vec<u8>> {
        let c = self.cfg.channels;
        let need = self.frame_samples() * c;
        if pcm.len() != need {
            return Err(crate::Error::BadArgument(format!("{} samples given, a packet takes {need}", pcm.len())));
        }
        let x48: Vec<f32> = if self.upsamplers.is_empty() {
            pcm.to_vec()
        } else {
            let mut chans = Vec::with_capacity(c);
            for ch in 0..c {
                let mono: Vec<f32> = pcm.iter().skip(ch).step_by(c).copied().collect();
                let mut y = Vec::new();
                self.upsamplers[ch].process(&mono, &mut y);
                chans.push(y);
            }
            let n = chans[0].len();
            let mut inter = vec![0.0; n * c];
            for (ch, y) in chans.iter().enumerate() {
                for (i, v) in y.iter().enumerate() {
                    inter[i * c + ch] = *v;
                }
            }
            inter
        };
        let mode = self.choose_mode();
        let bw = self.choose_bandwidth(mode);
        let n = self.cfg.frame_size;
        if self.prev_mode.is_some_and(|m| m != mode) {
            // Without redundancy the decoder conceals the switch (§4.5,
            // Figure 19); start both layers clean, as it does.
            if mode != Mode::Silk {
                self.celt.reset();
            }
            if self.prev_mode == Some(Mode::Celt) {
                self.silk.reset();
            }
        }
        let packet = match mode {
            Mode::Celt => self.encode_celt(&x48, n, bw)?,
            Mode::Silk | Mode::Hybrid => self.encode_silk_hybrid(&x48, n, bw, mode)?,
        };
        self.prev_mode = Some(mode);
        Ok(packet)
    }

    /// The byte budget of a packet of `n` samples at the target rate, with
    /// the VBR adjustment `boost` (1.0 = average).
    fn packet_bytes(&mut self, n: usize, boost: f32) -> usize {
        let bits = i64::from(self.cfg.bitrate) * n as i64 / 48000;
        if !self.cfg.vbr {
            return ((bits + 4) / 8).max(3) as usize;
        }
        let want = (bits as f32 * boost) as i64;
        // Pull the long-term rate back to the target.
        let adjusted = want + self.reservoir / 8;
        let lo = bits / 2;
        let hi = bits * 2;
        let use_bits = adjusted.clamp(lo, hi);
        self.reservoir += bits - use_bits;
        ((use_bits + 4) / 8).max(3) as usize
    }

    fn encode_celt(&mut self, x48: &[f32], n: usize, bw: Bandwidth) -> Result<Vec<u8>> {
        let c = self.cfg.channels;
        let delayed = self.celt_delay.process(x48);
        let sub = n.min(960);
        let count = n / sub;
        let boost = if self.cfg.vbr { vbr_boost(&delayed, c) } else { 1.0 };
        let total = self.packet_bytes(n, boost);
        let overhead = if count == 1 { 1 } else { 2 };
        let per_frame = ((total.saturating_sub(overhead)) / count).clamp(2, packet::MAX_FRAME_BYTES);
        let toc = Toc {
            config: Toc::config_for(Mode::Celt, bw, sub).expect("CELT config"),
            stereo: c == 2,
            code: 0,
        };
        let cfg = FrameConfig { start: 0, end: bw.celt_end_band(), bitrate: self.cfg.bitrate as i32 };
        let mut frames = Vec::with_capacity(count);
        for k in 0..count {
            let mut enc = RangeEncoder::new(per_frame);
            self.celt.encode(&delayed[k * sub * c..(k + 1) * sub * c], sub, &mut enc, cfg);
            self.final_range = enc.range();
            frames.push(enc.finish());
        }
        let refs: Vec<&[u8]> = frames.iter().map(|f| f.as_slice()).collect();
        packet::build(toc, &refs, None)
    }

    /// Downsamples to the SILK rate and keeps `SILK_LOOKAHEAD_MS` of
    /// lookahead; returns the frame plus lookahead, interleaved.
    fn silk_input(&mut self, x48: &[f32], fs_khz: usize) -> Vec<f32> {
        let c = self.cfg.channels;
        if self.silk_fs != fs_khz {
            // The SILK path's delay: X ms of lookahead buffer, the
            // downsampler, the decoder's one-sample stereo delay and its
            // resampler (Table 54) add up to the encoder's lookahead.
            let table54 = match fs_khz {
                8 => 0.538,
                12 => 0.692,
                _ => 0.706,
            };
            let total_ms = LOOKAHEAD_48K as f64 / 48.0;
            let rs_ms = total_ms - SILK_LOOKAHEAD_MS - 1.0 / fs_khz as f64 - table54;
            self.silk_down = (0..c).map(|_| Resampler::new(48000, fs_khz * 1000, rs_ms)).collect();
            let ahead = (SILK_LOOKAHEAD_MS * fs_khz as f64).round() as usize;
            self.silk_buf = vec![0.0; ahead * c];
            self.silk_fs = fs_khz;
        }
        let mut chans = Vec::with_capacity(c);
        for ch in 0..c {
            let mono: Vec<f32> = x48.iter().skip(ch).step_by(c).copied().collect();
            let mut y = Vec::new();
            self.silk_down[ch].process(&mono, &mut y);
            chans.push(y);
        }
        for i in 0..chans[0].len() {
            for y in &chans {
                self.silk_buf.push(y[i]);
            }
        }
        self.silk_buf.clone()
    }

    fn encode_silk_hybrid(&mut self, x48: &[f32], n: usize, bw: Bandwidth, mode: Mode) -> Result<Vec<u8>> {
        let c = self.cfg.channels;
        let fs_khz = if mode == Mode::Hybrid { 16 } else { crate::silk::decoder::fs_khz(bw) };
        let silk_in = self.silk_input(x48, fs_khz);
        let delayed = self.celt_delay.process(x48);
        let frame_len_total = n * fs_khz / 48;
        let boost = if self.cfg.vbr { vbr_boost(&delayed, c) } else { 1.0 };
        let total_bytes = self.packet_bytes(n, boost);
        let result = if mode == Mode::Silk {
            // One Opus frame of up to 60 ms.
            let budget = ((total_bytes - 1).min(packet::MAX_FRAME_BYTES) * 8) as i32;
            let mut enc = RangeEncoder::new(packet::MAX_FRAME_BYTES);
            self.silk.encode(&mut enc, &silk_in, frame_len_total, fs_khz, n / 48, budget, self.cfg.fec);
            // Size the frame to the bits used: fewer than 17 bits left over
            // tell the decoder there is no redundancy (§4.5.1.1).
            let used = ((enc.tell() + 7) >> 3) as usize;
            enc.shrink(used.max(1));
            self.final_range = enc.range();
            let frame = enc.finish();
            let toc = Toc { config: Toc::config_for(Mode::Silk, bw, n).expect("SILK config"), stereo: c == 2, code: 0 };
            let pad = if self.cfg.vbr { None } else { Some(total_bytes) };
            packet::build(toc, &[&frame], pad)
        } else {
            let sub = n.min(960);
            let count = n / sub;
            let overhead = if count == 1 { 1 } else { 2 };
            let per_frame = ((total_bytes.saturating_sub(overhead)) / count).clamp(8, packet::MAX_FRAME_BYTES);
            // The SILK layer's share (§2.1.1: hybrid splits the rate).
            let per_ch = self.cfg.bitrate as f64 / c as f64;
            let silk_rate = (0.55 * per_ch).clamp(8000.0, 22000.0) * c as f64;
            let silk_bits = ((silk_rate * sub as f64 / 48000.0) as i32).min(per_frame as i32 * 8 - 40);
            let sub_len = sub * fs_khz / 48;
            let ahead_len = silk_in.len() / c - frame_len_total;
            let mut frames = Vec::with_capacity(count);
            let cfg = FrameConfig { start: 17, end: bw.celt_end_band(), bitrate: self.cfg.bitrate as i32 };
            for k in 0..count {
                let mut enc = RangeEncoder::new(per_frame);
                let lo = k * sub_len * c;
                let hi = (k * sub_len + sub_len + ahead_len) * c;
                self.silk.encode(&mut enc, &silk_in[lo..hi.min(silk_in.len())], sub_len, fs_khz, sub / 48, silk_bits, self.cfg.fec);
                // The redundancy flag (§4.5.1.1), always off.
                if enc.tell() + 37 <= (per_frame * 8) as i32 {
                    enc.bit_logp(false, 12);
                }
                self.celt.encode(&delayed[k * sub * c..(k + 1) * sub * c], sub, &mut enc, cfg);
                self.final_range = enc.range();
                frames.push(enc.finish());
            }
            let toc = Toc { config: Toc::config_for(Mode::Hybrid, bw, sub).expect("hybrid config"), stereo: c == 2, code: 0 };
            let refs: Vec<&[u8]> = frames.iter().map(|f| f.as_slice()).collect();
            packet::build(toc, &refs, None)
        };
        self.silk_buf.drain(..frame_len_total * c);
        result
    }
}

/// A VBR weight for a frame: more bits for transients and loud, rich
/// frames, fewer for near-silence.
fn vbr_boost(x: &[f32], c: usize) -> f32 {
    let n = x.len() / c;
    if n == 0 {
        return 1.0;
    }
    let e: f32 = x.iter().map(|v| v * v).sum::<f32>() / x.len() as f32;
    let db = 10.0 * (e + 1e-12).log10();
    // Quiet frames need fewer bits.
    let mut w = if db < -70.0 { 0.3 } else if db < -50.0 { 0.6 } else { 1.0 };
    // A sharp onset inside the frame.
    let q = n / 4;
    if q > 0 {
        let seg = |k: usize| -> f32 { x[k * q * c..(k + 1) * q * c].iter().map(|v| v * v).sum::<f32>() + 1e-9 };
        let mut prev = seg(0);
        for k in 1..4 {
            let s = seg(k);
            if s > 10.0 * prev {
                w *= 1.4;
                break;
            }
            prev = s;
        }
    }
    w
}
