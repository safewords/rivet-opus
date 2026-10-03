//! The Opus decoder (RFC 6716 §4): SILK, CELT and hybrid frames, mode
//! transitions with and without redundancy (§4.5), packet loss concealment
//! (§4.4) and LBRR forward error correction (§4.2.5).

use crate::celt::decoder::CeltDecoder;
use crate::celt::mode::mode as celt_mode;
use crate::error::{Error, Result, config};
use crate::packet::{self, Bandwidth, Mode};
use crate::range::RangeDecoder;
use crate::resample::Resampler;
use crate::silk::SilkDecoder;

/// The sample rates Opus decodes to.
pub const SAMPLE_RATES: [u32; 5] = [8000, 12000, 16000, 24000, 48000];

/// The SILK resampler delay allocations of Table 54, in ms.
fn silk_delay_ms(fs_khz: usize) -> f64 {
    match fs_khz {
        8 => 0.538,
        12 => 0.692,
        _ => 0.706,
    }
}

/// An Opus decoder for one stream (mono or stereo).
pub struct Decoder {
    fs: u32,
    channels: usize,
    celt: CeltDecoder,
    silk: SilkDecoder,
    resamplers: Vec<Resampler>,
    prev_mode: Option<Mode>,
    prev_redundancy: bool,
    last_frame_size: usize,
    last_bw: Bandwidth,
    last_stereo: bool,
    final_range: u32,
    last_duration: usize,
}

impl Decoder {
    /// A decoder producing `channels` (1 or 2) interleaved channels at
    /// `sample_rate` (8, 12, 16, 24 or 48 kHz).
    pub fn new(sample_rate: u32, channels: usize) -> Result<Self> {
        if !SAMPLE_RATES.contains(&sample_rate) {
            return Err(config(format!("sample rate {sample_rate} is not one of 8, 12, 16, 24 or 48 kHz")));
        }
        if !(1..=2).contains(&channels) {
            return Err(config(format!("{channels} channels: a single Opus stream has 1 or 2")));
        }
        Ok(Self {
            fs: sample_rate,
            channels,
            celt: CeltDecoder::new(channels, (48000 / sample_rate) as usize),
            silk: SilkDecoder::new(),
            resamplers: Vec::new(),
            prev_mode: None,
            prev_redundancy: false,
            last_frame_size: 960,
            last_bw: Bandwidth::Full,
            last_stereo: channels == 2,
            final_range: 0,
            last_duration: 960,
        })
    }

    /// The output sample rate.
    pub fn sample_rate(&self) -> u32 {
        self.fs
    }

    /// The number of output channels.
    pub fn channels(&self) -> usize {
        self.channels
    }

    /// The range coder's final state after the last packet (RFC 6716 §6:
    /// a conforming decoder reproduces the reference's), XORed with that of
    /// the redundant frame, if any.
    pub fn final_range(&self) -> u32 {
        self.final_range
    }

    /// Duration in samples per channel (at the output rate) of the last
    /// packet decoded or concealed.
    pub fn last_packet_duration(&self) -> usize {
        self.last_duration * self.fs as usize / 48000
    }

    /// Whether intensity-stereo bands are decoded without their 180° phase
    /// inversion (RFC 8251 §10, for output that will be downmixed). Off by
    /// default for stereo output, on for mono.
    pub fn set_phase_inversion_disabled(&mut self, disabled: bool) {
        self.celt.disable_inv = disabled;
    }

    /// Resets the decoder to its initial state.
    pub fn reset(&mut self) {
        self.celt.reset();
        self.silk.reset();
        for r in &mut self.resamplers {
            r.reset();
        }
        self.prev_mode = None;
        self.prev_redundancy = false;
        self.final_range = 0;
    }

    /// Decodes a packet into interleaved samples at the output rate (±1.0
    /// full scale). `None` (a lost packet) conceals one packet's worth of
    /// audio, the duration of the last packet.
    pub fn decode(&mut self, packet: Option<&[u8]>) -> Result<Vec<f32>> {
        match packet {
            None | Some([]) => {
                let n = self.last_duration;
                Ok(self.conceal(n))
            }
            Some(data) => {
                let p = packet::parse(data)?;
                let toc = p.toc;
                let mut out = Vec::with_capacity(p.samples() * self.channels * self.fs as usize / 48000);
                for frame in &p.frames {
                    let pcm = self.decode_frame(Some(frame), toc.mode(), toc.bandwidth(), toc.frame_size(), toc.stereo, false);
                    out.extend_from_slice(&pcm);
                }
                self.last_duration = p.samples();
                Ok(out)
            }
        }
    }

    /// Recovers the packet before `next` from `next`'s LBRR data (forward
    /// error correction, §4.2.5): one frame of `next`'s frame size. When
    /// `next` carries no LBRR data for it, this is concealment.
    pub fn decode_fec(&mut self, next: &[u8]) -> Result<Vec<f32>> {
        let p = packet::parse(next)?;
        let toc = p.toc;
        let Some(first) = p.frames.first() else {
            return Err(Error::InvalidPacket("no frames".into()));
        };
        let n = toc.frame_size();
        if toc.mode() == Mode::Celt || first.is_empty() || self.prev_mode.is_none() {
            return Ok(self.conceal(n));
        }
        let pcm = self.decode_frame(Some(first), toc.mode(), toc.bandwidth(), n, toc.stereo, true);
        self.last_duration = n;
        Ok(pcm)
    }

    /// Conceals `n` samples (48 kHz) with the PLC of the last mode.
    fn conceal(&mut self, n: usize) -> Vec<f32> {
        let mut out = Vec::new();
        let mut left = n;
        while left > 0 {
            let chunk = match self.prev_mode {
                Some(Mode::Silk) | Some(Mode::Hybrid) => {
                    if left >= 960 { 960 } else { 480 }
                }
                _ => {
                    let mut c = 960;
                    while c > left && c > 120 {
                        c /= 2;
                    }
                    c
                }
            };
            let pcm = self.decode_frame(None, self.prev_mode.unwrap_or(Mode::Celt), self.last_bw, chunk, self.last_stereo, false);
            let take = chunk.min(left) * self.channels * self.fs as usize / 48000;
            out.extend_from_slice(&pcm[..take.min(pcm.len())]);
            left = left.saturating_sub(chunk);
        }
        self.last_duration = n;
        out
    }

    fn silk_to_pcm(&mut self, chans: Vec<Vec<f32>>, out: &mut [f32]) {
        let fs_khz = self.silk.fs_khz();
        let in_rate = fs_khz * 1000;
        if self.resamplers.first().is_none_or(|r| r.in_rate() != in_rate) {
            self.resamplers =
                (0..self.channels).map(|_| Resampler::new(in_rate, self.fs as usize, silk_delay_ms(fs_khz))).collect();
        }
        let cc = self.channels;
        for (c, ch) in chans.iter().enumerate() {
            let mut y = Vec::with_capacity(out.len() / cc);
            self.resamplers[c].process(ch, &mut y);
            for (i, v) in y.iter().enumerate() {
                if i * cc + c < out.len() {
                    out[i * cc + c] = *v;
                }
            }
        }
    }

    /// Decodes (or conceals, `data == None`) one frame of `n48` samples.
    fn decode_frame(
        &mut self,
        data: Option<&[u8]>,
        mode: Mode,
        bw: Bandwidth,
        n48: usize,
        stereo: bool,
        fec: bool,
    ) -> Vec<f32> {
        let cc = self.channels;
        let ds = 48000 / self.fs as usize;
        let n = n48 / ds;
        let f2_5 = 120 / ds;
        let f5 = 240 / ds;
        let mut pcm = vec![0.0f32; n * cc];
        let data = data.filter(|d| !d.is_empty());
        let lost = data.is_none();
        if lost && self.prev_mode.is_none() {
            return pcm;
        }
        let stream_c = if stereo { 2 } else { 1 };
        let mut ec_store = data.map(RangeDecoder::new);
        // A transition the stream does not cover with redundancy gets the
        // concealment of the old mode, cross-faded in (§4.5, Figure 19).
        let transition = !lost
            && !fec
            && self.prev_mode.is_some_and(|pm| {
                (mode == Mode::Celt && pm != Mode::Celt && !self.prev_redundancy) || (mode != Mode::Celt && pm == Mode::Celt)
            });
        // Into CELT, the old mode's concealment comes first; out of CELT it
        // is only needed when the frame carries no redundancy (decided
        // below).
        let mut pcm_transition = if transition && mode == Mode::Celt {
            Some(self.decode_frame(None, self.prev_mode.unwrap_or(Mode::Silk), self.last_bw, 480, self.last_stereo, false))
        } else {
            None
        };
        // SILK.
        if mode != Mode::Celt {
            if self.prev_mode == Some(Mode::Celt) {
                self.silk.reset();
                for r in &mut self.resamplers {
                    r.reset();
                }
            }
            let silk_bw = if mode == Mode::Hybrid { Bandwidth::Wide } else { bw };
            let frame_ms = n48 / 48;
            let chans = self.silk.decode(ec_store.as_mut(), stereo, silk_bw, frame_ms.max(10), cc, fec);
            self.silk_to_pcm(chans, &mut pcm);
        }
        // §4.5.1: redundancy.
        let mut redundancy = false;
        let mut celt_to_silk = false;
        let mut redundancy_bytes = 0usize;
        let mut len = data.map_or(0, |d| d.len());
        if let (Some(ec), Some(_)) = (ec_store.as_mut(), data)
            && mode != Mode::Celt
            && !fec
        {
            if mode == Mode::Hybrid {
                if ec.tell() + 17 + 20 <= 8 * len as i32 {
                    redundancy = ec.bit_logp(12);
                }
            } else {
                redundancy = ec.tell() + 17 <= 8 * len as i32;
            }
            if redundancy {
                celt_to_silk = ec.bit_logp(1);
                redundancy_bytes = if mode == Mode::Hybrid {
                    ec.uint(256) as usize + 2
                } else {
                    len - ((ec.tell() as usize + 7) >> 3)
                };
                if redundancy_bytes > len || ((len - redundancy_bytes) * 8) < ec.tell() as usize {
                    // An impossible size: treat as no redundancy.
                    len = 0;
                    redundancy_bytes = 0;
                    redundancy = false;
                } else {
                    len -= redundancy_bytes;
                    ec.shrink(redundancy_bytes);
                }
            }
        }
        if transition && mode != Mode::Celt && !redundancy {
            pcm_transition = Some(self.decode_frame(None, Mode::Celt, self.last_bw, 240, self.last_stereo, false));
        }
        let start_band = if mode == Mode::Hybrid { 17 } else { 0 };
        let end_band = bw.celt_end_band();
        let mut redundant_audio = Vec::new();
        let mut redundant_rng = 0u32;
        let red_slice = data.map(|d| &d[len..len + redundancy_bytes]);
        if redundancy && celt_to_silk {
            let mut red = vec![0.0f32; f5 * cc];
            let mut rec = RangeDecoder::new(red_slice.unwrap_or(&[]));
            self.celt.decode(&mut rec, 240, stream_c, 0, end_band, &mut red, false);
            redundant_rng = rec.range();
            redundant_audio = red;
        }
        // CELT.
        if mode != Mode::Silk {
            if self.prev_mode.is_some_and(|pm| pm != mode) && !self.prev_redundancy {
                self.celt.reset();
            }
            let accumulate = mode == Mode::Hybrid;
            match (ec_store.as_mut(), fec) {
                (Some(ec), false) => {
                    self.celt.decode(ec, n48.min(960), stream_c, start_band, end_band, &mut pcm, accumulate);
                }
                _ => self.celt.decode_lost(n48.min(960), &mut pcm, accumulate),
            }
        } else if !lost && self.prev_mode == Some(Mode::Hybrid) && !(redundancy && celt_to_silk && self.prev_redundancy) {
            // §4.5: Hybrid to SILK lets the CELT overlap fade out by decoding
            // a 2.5 ms silence frame.
            let silence = [0xFFu8, 0xFF];
            let mut sec = RangeDecoder::new(&silence);
            let mut tail = vec![0.0f32; f2_5 * cc];
            self.celt.decode(&mut sec, 120, stream_c, 0, end_band, &mut tail, false);
            for (o, t) in pcm.iter_mut().zip(&tail) {
                *o += t;
            }
        }
        let window = &celt_mode().window;
        let fade = |in1: &[f32], in2: &[f32], out: &mut [f32]| {
            for i in 0..f2_5 {
                let w = window[i * ds] * window[i * ds];
                for c in 0..cc {
                    let k = i * cc + c;
                    out[k] = w * in2[k] + (1.0 - w) * in1[k];
                }
            }
        };
        if redundancy && !celt_to_silk {
            self.celt.reset();
            let mut red = vec![0.0f32; f5 * cc];
            let mut rec = RangeDecoder::new(red_slice.unwrap_or(&[]));
            self.celt.decode(&mut rec, 240, stream_c, 0, end_band, &mut red, false);
            redundant_rng = rec.range();
            if n >= f2_5 {
                let at = (n - f2_5) * cc;
                let a = pcm[at..].to_vec();
                fade(&a, &red[f2_5 * cc..], &mut pcm[at..]);
            }
        }
        if redundancy && celt_to_silk && n >= f5 {
            pcm[..f2_5 * cc].copy_from_slice(&redundant_audio[..f2_5 * cc]);
            let b = pcm[f2_5 * cc..f5 * cc].to_vec();
            fade(&redundant_audio[f2_5 * cc..], &b, &mut pcm[f2_5 * cc..f5 * cc]);
        }
        if let Some(t) = pcm_transition {
            if n >= f5 {
                pcm[..f2_5 * cc].copy_from_slice(&t[..f2_5 * cc]);
                let b = pcm[f2_5 * cc..f5 * cc].to_vec();
                fade(&t[f2_5 * cc..f5 * cc], &b, &mut pcm[f2_5 * cc..f5 * cc]);
            } else {
                let b = pcm[..f2_5 * cc].to_vec();
                fade(&t[..f2_5 * cc], &b, &mut pcm[..f2_5 * cc]);
            }
        }
        if !lost {
            self.final_range = match ec_store.as_ref() {
                Some(ec) if data.is_some_and(|d| d.len() > 1) => ec.range() ^ redundant_rng,
                _ => 0,
            };
            self.prev_mode = Some(mode);
            self.prev_redundancy = redundancy && !celt_to_silk;
            self.last_bw = bw;
            self.last_stereo = stereo;
            self.last_frame_size = n48;
        }
        pcm
    }
}
