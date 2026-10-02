//! The SILK layer (RFC 6716 §4.2 and §5.2).

pub(crate) mod decoder;
pub(crate) mod encoder;
pub(crate) mod nlsf;
pub(crate) mod tables;

use crate::packet::Bandwidth;
use crate::range::RangeDecoder;
use decoder::{ChannelState, CondCoding, StereoState, decode_indices, decode_pulses, decode_stereo_weights, excitation_q23, fs_khz};
use tables::*;

/// The signal type of a SILK frame (Table 10).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SignalType {
    #[default]
    Inactive = 0,
    Unvoiced = 1,
    Voiced = 2,
}

impl SignalType {
    pub fn from_index(i: usize) -> Self {
        match i {
            0 => SignalType::Inactive,
            1 => SignalType::Unvoiced,
            _ => SignalType::Voiced,
        }
    }
}

/// The coded parameters of one SILK frame, as they appear in the
/// bitstream (§4.2.7).
#[derive(Clone, Debug, Default)]
pub struct FrameIndices {
    pub signal_type: SignalType,
    pub qoff: usize,
    /// Subframe 0: the absolute 6-bit index when coded independently, else
    /// the delta index; other subframes: delta indices.
    pub gains: [i32; 4],
    pub nlsf_i1: usize,
    pub nlsf_i2: [i32; 16],
    /// 4 when not coded (10 ms frames).
    pub interp_q2: i32,
    /// The primary pitch lag.
    pub lag: i32,
    pub contour: usize,
    pub periodicity: usize,
    pub ltp: [usize; 4],
    pub ltp_scale: usize,
    pub seed: u32,
}

/// The SILK decoder of one Opus stream.
pub struct SilkDecoder {
    ch: [ChannelState; 2],
    stereo: StereoState,
    fs_khz: usize,
    stream_channels: usize,
    side_prev_uncoded: bool,
}

/// A decoded (or parsed) frame of the LBRR or regular sequence.
struct Parsed {
    ix: FrameIndices,
    raw: Vec<i32>,
    cond: CondCoding,
}

impl SilkDecoder {
    pub fn new() -> Self {
        Self {
            ch: [ChannelState::new(16), ChannelState::new(16)],
            stereo: StereoState::default(),
            fs_khz: 0,
            stream_channels: 1,
            side_prev_uncoded: false,
        }
    }

    /// The decoder reset of §4.5.2.
    pub fn reset(&mut self) {
        let fs = if self.fs_khz == 0 { 16 } else { self.fs_khz };
        self.ch[0].reset(fs);
        self.ch[1].reset(fs);
        // RFC 8251 §3: the stereo state is reset too.
        self.stereo = StereoState::default();
        self.side_prev_uncoded = false;
    }

    /// The internal rate of the last frame, kHz.
    pub fn fs_khz(&self) -> usize {
        self.fs_khz
    }

    fn prepare(&mut self, stereo: bool, bw: Bandwidth) {
        let fs = fs_khz(bw);
        if fs != self.fs_khz {
            self.ch[0].reset(fs);
            self.ch[1].reset(fs);
            self.stereo = StereoState::default();
            self.side_prev_uncoded = false;
            self.fs_khz = fs;
        }
        let c = if stereo { 2 } else { 1 };
        if c == 2 && self.stream_channels == 1 {
            self.ch[1].reset(fs);
            self.stereo.prev_w = [0, 0];
            self.side_prev_uncoded = false;
        }
        self.stream_channels = c;
    }

    /// Decodes the SILK layer of one Opus frame of `frame_ms` ms (10, 20,
    /// 40 or 60) from `ec`, or conceals it when `ec` is `None`. With `fec`,
    /// the frame's LBRR data is decoded instead of its regular frames (the
    /// redundant copy of the previous packet's audio). Returns
    /// `out_channels` channels at the internal rate.
    pub fn decode(
        &mut self,
        ec: Option<&mut RangeDecoder>,
        stereo: bool,
        bw: Bandwidth,
        frame_ms: usize,
        out_channels: usize,
        fec: bool,
    ) -> Vec<Vec<f32>> {
        self.prepare(stereo, bw);
        let fs = self.fs_khz;
        let c = self.stream_channels;
        let nb_subfr = if frame_ms == 10 { 2 } else { 4 };
        let nf = (frame_ms / 20).max(1);
        let mut outs: Vec<Vec<f32>> = vec![Vec::new(); out_channels];
        let Some(ec) = ec else {
            for _ in 0..nf {
                let mid = self.ch[0].conceal(nb_subfr);
                let side = if c == 2 { self.ch[1].conceal(nb_subfr) } else { vec![0.0; mid.len()] };
                let w = self.stereo.prev_w;
                self.emit(&mid, &side, c, w, out_channels, &mut outs);
            }
            return outs;
        };
        // §4.2.3: header bits.
        let mut vad = [[false; 3]; 2];
        let mut lbrr = [false; 2];
        for ch in 0..c {
            for v in vad[ch].iter_mut().take(nf) {
                *v = ec.bit_logp(1);
            }
            lbrr[ch] = ec.bit_logp(1);
        }
        // §4.2.4: per-frame LBRR flags.
        let mut lbrr_flags = [[false; 3]; 2];
        for ch in 0..c {
            if lbrr[ch] {
                if nf == 1 {
                    lbrr_flags[ch][0] = true;
                } else {
                    let v = if nf == 2 { ec.icdf(&LBRR_FLAGS_2_ICDF, 8) } else { ec.icdf(&LBRR_FLAGS_3_ICDF, 8) } + 1;
                    for (i, f) in lbrr_flags[ch].iter_mut().enumerate().take(nf) {
                        *f = (v >> i) & 1 != 0;
                    }
                }
            }
        }
        // §4.2.5: LBRR frames.
        let frame_len = 5 * fs * nb_subfr;
        let mut lbrr_frames: Vec<[Option<Parsed>; 2]> = (0..nf).map(|_| [None, None]).collect();
        let mut lbrr_w = vec![[0i32; 2]; nf];
        let mut lbrr_mid_only = vec![false; nf];
        let mut prev_type = [SignalType::Inactive; 2];
        let mut prev_lag = [0i32; 2];
        for i in 0..nf {
            for ch in 0..c {
                if !lbrr_flags[ch][i] {
                    continue;
                }
                if ch == 0 && c == 2 {
                    lbrr_w[i] = decode_stereo_weights(ec);
                    if !lbrr_flags[1][i] {
                        lbrr_mid_only[i] = ec.icdf(&MID_ONLY_ICDF, 8) == 1;
                    }
                }
                let cond = if i > 0 && lbrr_flags[ch][i - 1] { CondCoding::Conditional } else { CondCoding::Independent };
                let ix = decode_indices(ec, fs, nb_subfr, true, cond, prev_type[ch], prev_lag[ch]);
                let raw = decode_pulses(ec, ix.signal_type, ix.qoff, frame_len);
                prev_type[ch] = ix.signal_type;
                prev_lag[ch] = ix.lag;
                lbrr_frames[i][ch] = Some(Parsed { ix, raw, cond });
            }
        }
        if fec {
            for (i, frames) in lbrr_frames.iter().enumerate() {
                let mid = match &frames[0] {
                    Some(p) => self.synth_parsed(0, p, nb_subfr, frame_len),
                    None => self.ch[0].conceal(nb_subfr),
                };
                let side = if c == 2 {
                    match &frames[1] {
                        Some(p) => self.synth_parsed(1, p, nb_subfr, frame_len),
                        None if lbrr_mid_only[i] && frames[0].is_some() => vec![0.0; mid.len()],
                        None => self.ch[1].conceal(nb_subfr),
                    }
                } else {
                    vec![0.0; mid.len()]
                };
                let w = if frames[0].is_some() { lbrr_w[i] } else { self.stereo.prev_w };
                self.emit(&mid, &side, c, w, out_channels, &mut outs);
            }
            return outs;
        }
        // §4.2.6: regular frames.
        for i in 0..nf {
            let mut w = [0i32; 2];
            let mut mid_only = false;
            let mut mid = Vec::new();
            let mut side = Vec::new();
            for ch in 0..c {
                if ch == 0 && c == 2 {
                    w = decode_stereo_weights(ec);
                    if !vad[1][i] {
                        mid_only = ec.icdf(&MID_ONLY_ICDF, 8) == 1;
                    }
                }
                if ch == 1 && mid_only {
                    side = vec![0.0; frame_len];
                    continue;
                }
                let mut cond = if i == 0 { CondCoding::Independent } else { CondCoding::Conditional };
                if ch == 1 && self.side_prev_uncoded {
                    // §4.2.7.9: after an uncoded side frame the side channel
                    // starts from a clean state.
                    self.ch[1].reset(fs);
                    if i > 0 {
                        cond = CondCoding::IndependentNoLtpScaling;
                    }
                }
                let st = &self.ch[ch];
                let ix = decode_indices(ec, fs, nb_subfr, vad[ch][i], cond, st.prev_signal_type, st.prev_lag);
                let raw = decode_pulses(ec, ix.signal_type, ix.qoff, frame_len);
                let out = self.synth_parsed(ch, &Parsed { ix, raw, cond }, nb_subfr, frame_len);
                if ch == 0 {
                    mid = out;
                } else {
                    side = out;
                    self.side_prev_uncoded = false;
                }
            }
            if c == 2 && mid_only {
                self.side_prev_uncoded = true;
            }
            if c == 1 {
                side = vec![0.0; mid.len()];
            }
            self.emit(&mid, &side, c, w, out_channels, &mut outs);
        }
        outs
    }

    fn synth_parsed(&mut self, ch: usize, p: &Parsed, nb_subfr: usize, frame_len: usize) -> Vec<f32> {
        let st = &mut self.ch[ch];
        let params = st.dequantize(&p.ix, nb_subfr, p.cond);
        let e = excitation_q23(&p.raw, p.ix.signal_type, p.ix.qoff, p.ix.seed, frame_len);
        st.synthesize(&params, &e, nb_subfr)
    }

    fn emit(&mut self, mid: &[f32], side: &[f32], c: usize, w: [i32; 2], out_channels: usize, outs: &mut [Vec<f32>]) {
        if out_channels == 2 && c == 2 {
            let (l, r) = self.stereo.unmix(mid, side, w, self.fs_khz);
            outs[0].extend_from_slice(&l);
            outs[1].extend_from_slice(&r);
        } else {
            if c == 2 {
                self.stereo.prev_w = w;
            }
            let m = self.stereo.delay_mono(mid);
            for o in outs.iter_mut() {
                o.extend_from_slice(&m);
            }
        }
    }
}
