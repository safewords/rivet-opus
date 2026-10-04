//! The SILK decoder (RFC 6716 §4.2): header and LBRR flags, the frame
//! parameters, the excitation, LTP and LPC synthesis, stereo unmixing, and
//! packet loss concealment.

use super::nlsf;
use super::tables::*;
use super::{FrameIndices, SignalType};
use crate::packet::Bandwidth;
use crate::range::RangeDecoder;

/// Samples of clamped output kept for the LTP rewhitening (§4.2.7.9.1).
const OUT_HIST: usize = 320;

/// Conditional coding of a SILK frame's parameters (§4.2.7.4, §4.2.7.6).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CondCoding {
    /// Absolute gain, absolute pitch, LTP scaling coded.
    Independent,
    /// Absolute gain and pitch, no LTP scaling (a side frame after an
    /// uncoded one within the Opus frame).
    IndependentNoLtpScaling,
    /// Delta gain, delta pitch when the previous frame was voiced.
    Conditional,
}

/// The rate of a SILK bandwidth in kHz.
pub fn fs_khz(bw: Bandwidth) -> usize {
    match bw {
        Bandwidth::Narrow => 8,
        Bandwidth::Medium => 12,
        _ => 16,
    }
}

/// Pitch lag limits and low-part table for a rate (Table 30).
fn pitch_params(fs_khz: usize) -> (i32, i32, i32, &'static [u8]) {
    match fs_khz {
        8 => (4, 16, 144, &PITCH_LOW_NB_ICDF),
        12 => (6, 24, 216, &PITCH_LOW_MB_ICDF),
        _ => (8, 32, 288, &PITCH_LOW_WB_ICDF),
    }
}

/// `silk_log2lin` (§4.2.7.4): approximately `2^(x/128)`.
pub fn log2lin(x: i32) -> i32 {
    let i = x >> 7;
    let f = x & 127;
    (1 << i) + ((((-174 * f * (128 - f)) >> 16) + f) * ((1 << i) >> 7))
}

/// The Q16 gain of a gain index (§4.2.7.4).
pub fn gain_q16(log_gain: i32) -> i32 {
    log2lin(((0x1D_1C71i64 * i64::from(log_gain)) >> 16) as i32 + 2090)
}

/// Decodes the stereo prediction weights (§4.2.7.1), Q13.
pub fn decode_stereo_weights(ec: &mut RangeDecoder) -> [i32; 2] {
    let n = ec.icdf(&STEREO_STAGE1_ICDF, 8) as i32;
    let i0 = ec.icdf(&STEREO_STAGE2_ICDF, 8) as i32;
    let i1 = ec.icdf(&STEREO_STAGE3_ICDF, 8) as i32;
    let i2 = ec.icdf(&STEREO_STAGE2_ICDF, 8) as i32;
    let i3 = ec.icdf(&STEREO_STAGE3_ICDF, 8) as i32;
    stereo_weights(n, i0, i1, i2, i3)
}

/// Table 7 interpolation of the stereo weights from their indices.
pub fn stereo_weights(n: i32, i0: i32, i1: i32, i2: i32, i3: i32) -> [i32; 2] {
    let w = &STEREO_WEIGHTS_Q13;
    let wi0 = (i0 + 3 * (n / 5)) as usize;
    let wi1 = (i2 + 3 * (n % 5)) as usize;
    let w1 = w[wi1] + (((w[wi1 + 1] - w[wi1]) * 6554) >> 16) * (2 * i3 + 1);
    let w0 = w[wi0] + (((w[wi0 + 1] - w[wi0]) * 6554) >> 16) * (2 * i1 + 1) - w1;
    [w0, w1]
}

/// Decodes one SILK frame's parameters (§4.2.7.3–§4.2.7.7).
pub fn decode_indices(
    ec: &mut RangeDecoder,
    fs_khz: usize,
    nb_subfr: usize,
    active: bool,
    cond: CondCoding,
    prev_signal_type: SignalType,
    prev_lag: i32,
) -> FrameIndices {
    let mut ix = FrameIndices::default();
    let ftype = if active {
        ec.icdf(&FRAME_TYPE_ACTIVE_ICDF, 8) + 2
    } else {
        ec.icdf(&FRAME_TYPE_INACTIVE_ICDF, 8)
    };
    ix.signal_type = SignalType::from_index(ftype >> 1);
    ix.qoff = ftype & 1;
    for k in 0..nb_subfr {
        if k == 0 && cond != CondCoding::Conditional {
            let msb = ec.icdf(&GAIN_MSB_ICDF[ix.signal_type as usize], 8) as i32;
            let lsb = ec.icdf(&GAIN_LSB_ICDF, 8) as i32;
            ix.gains[k] = (msb << 3) | lsb;
        } else {
            ix.gains[k] = ec.icdf(&GAIN_DELTA_ICDF, 8) as i32;
        }
    }
    let wb = fs_khz == 16;
    let s1 = usize::from(wb) * 2 + usize::from(ix.signal_type == SignalType::Voiced);
    ix.nlsf_i1 = ec.icdf(&NLSF_STAGE1_ICDF[s1], 8);
    for k in 0..nlsf::order(wb) {
        let mut v = ec.icdf(nlsf::stage2_icdf(wb, ix.nlsf_i1, k), 8) as i32 - 4;
        if v == -4 {
            v -= ec.icdf(&NLSF_EXT_ICDF, 8) as i32;
        } else if v == 4 {
            v += ec.icdf(&NLSF_EXT_ICDF, 8) as i32;
        }
        ix.nlsf_i2[k] = v;
    }
    ix.interp_q2 = if nb_subfr == 4 {
        ec.icdf(&NLSF_INTERP_ICDF, 8) as i32
    } else {
        4
    };
    if ix.signal_type == SignalType::Voiced {
        let (scale, min_lag, _, low_icdf) = pitch_params(fs_khz);
        let mut absolute = true;
        if cond == CondCoding::Conditional && prev_signal_type == SignalType::Voiced {
            let delta = ec.icdf(&PITCH_DELTA_ICDF, 8) as i32;
            if delta > 0 {
                ix.lag = prev_lag + delta - 9;
                absolute = false;
            }
        }
        if absolute {
            let high = ec.icdf(&PITCH_HIGH_ICDF, 8) as i32;
            let low = ec.icdf(low_icdf, 8) as i32;
            ix.lag = high * scale + low + min_lag;
        }
        ix.contour = ec.icdf(contour_icdf(fs_khz, nb_subfr), 8);
        ix.periodicity = ec.icdf(&PERIODICITY_ICDF, 8);
        for k in 0..nb_subfr {
            ix.ltp[k] = ec.icdf(LTP_FILTER_ICDF[ix.periodicity], 8);
        }
        ix.ltp_scale = if cond == CondCoding::Independent {
            ec.icdf(&LTP_SCALE_ICDF, 8)
        } else {
            0
        };
    }
    ix.seed = ec.icdf(&SEED_ICDF, 8) as u32;
    ix
}

/// The pitch contour PDF for a rate and frame size (Table 32).
pub fn contour_icdf(fs_khz: usize, nb_subfr: usize) -> &'static [u8] {
    match (fs_khz == 8, nb_subfr == 4) {
        (true, false) => &CONTOUR_NB10_ICDF,
        (true, true) => &CONTOUR_NB20_ICDF,
        (false, false) => &CONTOUR_WB10_ICDF,
        (false, true) => &CONTOUR_WB20_ICDF,
    }
}

/// The per-subframe lag offsets of a contour index (Tables 33–36).
pub fn contour_offsets(fs_khz: usize, nb_subfr: usize, idx: usize) -> [i32; 4] {
    let mut o = [0; 4];
    match (fs_khz == 8, nb_subfr == 4) {
        (true, false) => o[..2].copy_from_slice(&CONTOUR_NB10[idx]),
        (true, true) => o.copy_from_slice(&CONTOUR_NB20[idx]),
        (false, false) => o[..2].copy_from_slice(&CONTOUR_WB10[idx]),
        (false, true) => o.copy_from_slice(&CONTOUR_WB20[idx]),
    }
    o
}

/// The LTP taps (Q7) of a periodicity and filter index (Tables 39–41).
pub fn ltp_taps(periodicity: usize, idx: usize) -> [i32; 5] {
    match periodicity {
        0 => LTP_TAPS_0[idx],
        1 => LTP_TAPS_1[idx],
        _ => LTP_TAPS_2[idx],
    }
}

/// The Q14 LTP scale factors (§4.2.7.6.3).
pub const LTP_SCALES_Q14: [i32; 3] = [15565, 12288, 8192];

fn shell_split(ec: &mut RangeDecoder, total: i32, table: &[&[u8]]) -> (i32, i32) {
    if total == 0 {
        return (0, 0);
    }
    let l = ec.icdf(table[total as usize - 1], 8) as i32;
    (l, total - l)
}

/// §4.2.7.8.3: the pulse locations of one 16-sample block.
fn shell_decode(ec: &mut RangeDecoder, total: i32, out: &mut [i32]) {
    let (a, b) = shell_split(ec, total, &SHELL16_ICDF);
    for (half, cnt) in [(0usize, a), (8, b)] {
        let (c, d) = shell_split(ec, cnt, &SHELL8_ICDF);
        for (q, cnt4) in [(half, c), (half + 4, d)] {
            let (e, f) = shell_split(ec, cnt4, &SHELL4_ICDF);
            for (r, cnt2) in [(q, e), (q + 2, f)] {
                let (g, h) = shell_split(ec, cnt2, &SHELL2_ICDF);
                out[r] = g;
                out[r + 1] = h;
            }
        }
    }
}

/// §4.2.7.8: the raw excitation (signed pulse magnitudes) of a frame of
/// `frame_len` samples (a multiple of 16 is decoded).
pub fn decode_pulses(
    ec: &mut RangeDecoder,
    signal_type: SignalType,
    qoff: usize,
    frame_len: usize,
) -> Vec<i32> {
    let blocks = frame_len.div_ceil(16);
    let rate_level = ec.icdf(
        &RATE_LEVEL_ICDF[usize::from(signal_type == SignalType::Voiced)],
        8,
    );
    let mut counts = vec![0i32; blocks];
    let mut lsbs = vec![0usize; blocks];
    for b in 0..blocks {
        let mut c = ec.icdf(&PULSE_COUNT_ICDF[rate_level], 8);
        while c == 17 {
            lsbs[b] += 1;
            let level = if lsbs[b] == 10 { 10 } else { 9 };
            c = ec.icdf(&PULSE_COUNT_ICDF[level], 8);
        }
        counts[b] = c as i32;
    }
    let mut e = vec![0i32; blocks * 16];
    for b in 0..blocks {
        if counts[b] > 0 {
            shell_decode(ec, counts[b], &mut e[b * 16..b * 16 + 16]);
        }
    }
    for b in 0..blocks {
        if lsbs[b] > 0 {
            for v in &mut e[b * 16..b * 16 + 16] {
                for _ in 0..lsbs[b] {
                    *v = 2 * *v + ec.icdf(&LSB_ICDF, 8) as i32;
                }
            }
        }
    }
    let group = signal_type as usize * 2 + qoff;
    for b in 0..blocks {
        let icdf = &SIGN_ICDF[group][counts[b].min(6) as usize];
        for v in &mut e[b * 16..b * 16 + 16] {
            if *v != 0 && ec.icdf(icdf, 8) == 0 {
                *v = -*v;
            }
        }
    }
    e
}

/// §4.2.7.8.6: the excitation in Q23 from the raw pulses and the seed.
pub fn excitation_q23(
    raw: &[i32],
    signal_type: SignalType,
    qoff: usize,
    seed: u32,
    len: usize,
) -> Vec<i32> {
    let offset = QUANT_OFFSETS_Q23[signal_type as usize][qoff];
    let mut seed = seed;
    raw[..len]
        .iter()
        .map(|&r| {
            let mut e = (r << 8) - r.signum() * 20 + offset;
            seed = seed.wrapping_mul(196_314_165).wrapping_add(907_633_515);
            if seed & 0x8000_0000 != 0 {
                e = -e;
            }
            seed = seed.wrapping_add(r as u32);
            e
        })
        .collect()
}

/// The dequantized parameters of one frame, ready for synthesis.
#[derive(Clone, Debug)]
pub struct FrameParams {
    pub signal_type: SignalType,
    pub gains_q16: [i32; 4],
    /// LPC coefficients (Q12) for the first and second halves.
    pub a_q12: [[i32; 16]; 2],
    /// Whether the first half uses the interpolated LSFs.
    pub interp: bool,
    pub pitch_lags: [i32; 4],
    pub ltp_taps: [[i32; 5]; 4],
    pub ltp_scale_q14: i32,
}

/// One channel's SILK state.
#[derive(Clone)]
pub struct ChannelState {
    pub fs_khz: usize,
    pub lpc_order: usize,
    pub prev_nlsf: [i32; 16],
    pub first_frame_after_reset: bool,
    pub last_gain_index: i32,
    pub prev_signal_type: SignalType,
    pub prev_lag: i32,
    pub(crate) out_hist: Vec<f32>,
    pub(crate) lpc_hist: [f32; 16],
    // Concealment state.
    plc_a: [f32; 16],
    plc_exc: Vec<f32>,
    plc_lag: usize,
    plc_gain: f32,
    plc_voiced: bool,
    plc_count: u32,
    plc_seed: u32,
}

impl ChannelState {
    pub fn new(fs_khz: usize) -> Self {
        let mut s = Self {
            fs_khz,
            lpc_order: if fs_khz == 16 { 16 } else { 10 },
            prev_nlsf: [0; 16],
            first_frame_after_reset: true,
            last_gain_index: 10,
            prev_signal_type: SignalType::Inactive,
            prev_lag: 100,
            out_hist: vec![0.0; OUT_HIST],
            lpc_hist: [0.0; 16],
            plc_a: [0.0; 16],
            plc_exc: Vec::new(),
            plc_lag: 0,
            plc_gain: 0.0,
            plc_voiced: false,
            plc_count: 0,
            plc_seed: 22222,
        };
        s.reset(fs_khz);
        s
    }

    /// Back to the state after a decoder reset (§4.5.2) at `fs_khz`.
    pub fn reset(&mut self, fs_khz: usize) {
        self.fs_khz = fs_khz;
        self.lpc_order = if fs_khz == 16 { 16 } else { 10 };
        self.prev_nlsf = [0; 16];
        self.first_frame_after_reset = true;
        self.last_gain_index = 10;
        self.prev_signal_type = SignalType::Inactive;
        self.prev_lag = 100;
        self.out_hist.iter_mut().for_each(|v| *v = 0.0);
        self.lpc_hist = [0.0; 16];
        self.plc_a = [0.0; 16];
        self.plc_exc.clear();
        self.plc_lag = 0;
        self.plc_gain = 0.0;
        self.plc_voiced = false;
        self.plc_count = 0;
    }

    /// Dequantizes the parameters of a decoded frame and advances the
    /// parameter state (§4.2.7.4–§4.2.7.6).
    pub fn dequantize(
        &mut self,
        ix: &FrameIndices,
        nb_subfr: usize,
        cond: CondCoding,
    ) -> FrameParams {
        let mut gains_q16 = [0i32; 4];
        let mut prev = self.last_gain_index;
        for k in 0..nb_subfr {
            if k == 0 && cond != CondCoding::Conditional {
                prev = ix.gains[0].max(prev - 16);
            } else {
                let d = ix.gains[k];
                prev = (2 * d - 16).max(prev + d - 4).clamp(0, 63);
            }
            gains_q16[k] = gain_q16(prev);
        }
        self.last_gain_index = prev;
        let wb = self.fs_khz == 16;
        let d = self.lpc_order;
        let mut n2 = nlsf::reconstruct(wb, ix.nlsf_i1, &ix.nlsf_i2);
        nlsf::stabilize(&mut n2[..d], wb);
        let w_q2 = if self.first_frame_after_reset {
            4
        } else {
            ix.interp_q2
        };
        let interp = nb_subfr == 4 && w_q2 < 4;
        let mut a_q12 = [[0i32; 16]; 2];
        let a2 = nlsf::nlsf_to_lpc(&n2[..d], wb);
        a_q12[1][..d].copy_from_slice(&a2);
        if interp {
            let n1 = nlsf::interpolate(&self.prev_nlsf, &n2, w_q2, d);
            let a1 = nlsf::nlsf_to_lpc(&n1[..d], wb);
            a_q12[0][..d].copy_from_slice(&a1);
        } else {
            a_q12[0] = a_q12[1];
        }
        self.prev_nlsf = n2;
        self.first_frame_after_reset = false;
        let mut pitch_lags = [0i32; 4];
        let mut ltp = [[0i32; 5]; 4];
        let mut ltp_scale_q14 = LTP_SCALES_Q14[0];
        if ix.signal_type == SignalType::Voiced {
            let (_, min_lag, max_lag, _) = pitch_params(self.fs_khz);
            let off = contour_offsets(self.fs_khz, nb_subfr, ix.contour);
            for k in 0..nb_subfr {
                pitch_lags[k] = (ix.lag + off[k]).clamp(min_lag, max_lag);
                ltp[k] = ltp_taps(ix.periodicity, ix.ltp[k]);
            }
            ltp_scale_q14 = LTP_SCALES_Q14[ix.ltp_scale];
            self.prev_lag = ix.lag;
        }
        self.prev_signal_type = ix.signal_type;
        FrameParams {
            signal_type: ix.signal_type,
            gains_q16,
            a_q12,
            interp,
            pitch_lags,
            ltp_taps: ltp,
            ltp_scale_q14,
        }
    }

    /// §4.2.7.9: LTP and LPC synthesis of one frame from its excitation;
    /// returns the clamped output.
    pub fn synthesize(&mut self, p: &FrameParams, e_q23: &[i32], nb_subfr: usize) -> Vec<f32> {
        self.synthesize_with(p, nb_subfr, &mut |i, _, _, _| e_q23[i])
    }

    /// [`Self::synthesize`] with the excitation chosen sample by sample:
    /// `choose(i, ltp, lpc, gain)` gets the frame position, the LTP
    /// prediction of the residual, the LPC prediction of the output and the
    /// subframe gain (Q16 as a float), and returns `e_Q23[i]`. The encoder's
    /// closed-loop quantizer runs through this, so it tracks the decoder
    /// exactly.
    pub fn synthesize_with(
        &mut self,
        p: &FrameParams,
        nb_subfr: usize,
        choose: &mut dyn FnMut(usize, f32, f32, f32) -> i32,
    ) -> Vec<f32> {
        let d = self.lpc_order;
        let n = 5 * self.fs_khz;
        let frame_len = n * nb_subfr;
        let h = OUT_HIST;
        // out[h + i] for frame sample i; lpc[16 + i].
        let mut out = vec![0.0f32; h + frame_len];
        out[..h].copy_from_slice(&self.out_hist);
        let mut lpc = vec![0.0f32; 16 + frame_len];
        lpc[..16].copy_from_slice(&self.lpc_hist);
        let mut exc_all = vec![0.0f32; frame_len];
        for s in 0..nb_subfr {
            let j = s * n;
            let a: Vec<f32> = p.a_q12[usize::from(!(s < 2 && p.interp))][..d]
                .iter()
                .map(|&v| v as f32 / 4096.0)
                .collect();
            let gain = p.gains_q16[s] as f32;
            let mut res = vec![0.0f32; n];
            if p.signal_type == SignalType::Voiced {
                let lag = p.pitch_lags[s] as usize;
                let (out_end, scale) = if s >= 2 && p.interp {
                    (j as isize - ((s - 2) * n) as isize, 16384.0f32)
                } else {
                    (j as isize - (s * n) as isize, p.ltp_scale_q14 as f32)
                };
                let start = j as isize - lag as isize - 2;
                let mut hist = vec![0.0f32; lag + 2];
                for i in start..j as isize {
                    let idx = (i - start) as usize;
                    if i < out_end {
                        let oi = (h as isize + i) as usize;
                        let mut acc = out[oi];
                        for (k, ak) in a.iter().enumerate() {
                            acc -= out[oi - k - 1] * ak;
                        }
                        hist[idx] = 4.0 * scale / gain * acc.clamp(-1.0, 1.0);
                    } else {
                        let li = (16 + i) as usize;
                        let mut acc = lpc[li];
                        for (k, ak) in a.iter().enumerate() {
                            acc -= lpc[li - k - 1] * ak;
                        }
                        hist[idx] = 65536.0 / gain * acc;
                    }
                }
                // res over [start, j + n).
                let mut r = hist;
                r.resize(lag + 2 + n, 0.0);
                let b: Vec<f32> = p.ltp_taps[s].iter().map(|&v| v as f32 / 128.0).collect();
                for i in 0..n {
                    let ri = lag + 2 + i;
                    let mut ltp = 0.0f32;
                    for (k, bk) in b.iter().enumerate() {
                        ltp += r[ri - lag + 2 - k] * bk;
                    }
                    let li = 16 + j + i;
                    let mut pred = 0.0f32;
                    for (k, ak) in a.iter().enumerate() {
                        pred += lpc[li - k - 1] * ak;
                    }
                    let e = choose(j + i, ltp, pred, gain);
                    let v = e as f32 / 8_388_608.0 + ltp;
                    r[ri] = v;
                    res[i] = v;
                    let exc = gain / 65536.0 * v;
                    exc_all[j + i] = exc;
                    let o = exc + pred;
                    lpc[li] = o;
                    out[h + j + i] = o.clamp(-1.0, 1.0);
                }
            } else {
                for i in 0..n {
                    let li = 16 + j + i;
                    let mut pred = 0.0f32;
                    for (k, ak) in a.iter().enumerate() {
                        pred += lpc[li - k - 1] * ak;
                    }
                    let e = choose(j + i, 0.0, pred, gain);
                    let v = e as f32 / 8_388_608.0;
                    res[i] = v;
                    let exc = gain / 65536.0 * v;
                    exc_all[j + i] = exc;
                    let o = exc + pred;
                    lpc[li] = o;
                    out[h + j + i] = o.clamp(-1.0, 1.0);
                }
            }
        }
        self.out_hist
            .copy_from_slice(&out[frame_len..frame_len + h]);
        self.lpc_hist
            .copy_from_slice(&lpc[frame_len..frame_len + 16]);
        // Remember what the concealment needs.
        let a_last = &p.a_q12[1];
        for k in 0..16 {
            self.plc_a[k] = a_last[k] as f32 / 4096.0;
        }
        self.plc_voiced = p.signal_type == SignalType::Voiced;
        self.plc_lag = if self.plc_voiced {
            p.pitch_lags[nb_subfr - 1] as usize
        } else {
            0
        };
        self.plc_exc = exc_all;
        self.plc_gain = 1.0;
        self.plc_count = 0;
        out[h..].to_vec()
    }

    /// Concealment of a lost frame of `nb_subfr` subframes (§4.4): the last
    /// excitation repeated at the pitch period (voiced) or as noise
    /// (unvoiced), decaying, through the last LPC filter.
    pub fn conceal(&mut self, nb_subfr: usize) -> Vec<f32> {
        let d = self.lpc_order;
        let n = 5 * self.fs_khz * nb_subfr;
        self.plc_count += 1;
        let decay = if self.plc_voiced { 0.92f32 } else { 0.8 };
        let mut out = vec![0.0f32; n];
        if self.plc_exc.is_empty() {
            self.out_hist.copy_within(n.min(OUT_HIST).., 0);
            return out;
        }
        let energy: f32 =
            self.plc_exc.iter().map(|v| v * v).sum::<f32>() / self.plc_exc.len() as f32;
        let rms = energy.sqrt();
        let lag = self.plc_lag.clamp(1, self.plc_exc.len());
        let hist = self.plc_exc.clone();
        let mut exc = Vec::with_capacity(n);
        let g0 = self.plc_gain;
        let g1 = g0 * decay.powi(nb_subfr as i32 / 2 + 1);
        for i in 0..n {
            let g = g0 + (g1 - g0) * i as f32 / n as f32;
            let v = if self.plc_voiced {
                hist[hist.len() - lag + (i % lag)]
            } else {
                self.plc_seed = self
                    .plc_seed
                    .wrapping_mul(196_314_165)
                    .wrapping_add(907_633_515);
                ((self.plc_seed >> 16) as f32 / 32768.0 - 1.0) * rms * 1.7
            };
            exc.push(v * g);
        }
        self.plc_gain = g1;
        let mut lpc = vec![0.0f32; 16 + n];
        lpc[..16].copy_from_slice(&self.lpc_hist);
        let mut bw = 1.0f32;
        let a: Vec<f32> = (0..d)
            .map(|k| {
                bw *= 0.99;
                self.plc_a[k] * bw
            })
            .collect();
        for i in 0..n {
            let mut v = exc[i];
            for (k, ak) in a.iter().enumerate() {
                v += lpc[16 + i - k - 1] * ak;
            }
            lpc[16 + i] = v;
            out[i] = v.clamp(-1.0, 1.0);
        }
        self.lpc_hist.copy_from_slice(&lpc[n..n + 16]);
        let mut full = self.out_hist.clone();
        full.extend_from_slice(&out);
        self.out_hist
            .copy_from_slice(&full[full.len() - OUT_HIST..]);
        // Keep the periodic part going for the next lost frame.
        let mut ex = self.plc_exc.clone();
        ex.extend_from_slice(&exc);
        let keep = self.plc_exc.len();
        self.plc_exc = ex[ex.len() - keep..].to_vec();
        out
    }
}

/// The stereo unmixing state (§4.2.8).
#[derive(Clone, Default)]
pub struct StereoState {
    pub prev_w: [i32; 2],
    /// mid[i-2], mid[i-1].
    mid_hist: [f32; 2],
    side_hist: f32,
}

impl StereoState {
    /// Converts one frame of mid/side to left/right with a one-sample delay.
    pub fn unmix(
        &mut self,
        mid: &[f32],
        side: &[f32],
        w: [i32; 2],
        fs_khz: usize,
    ) -> (Vec<f32>, Vec<f32>) {
        let n2 = mid.len();
        let n1 = 8 * fs_khz;
        let mut m = Vec::with_capacity(n2 + 2);
        m.extend_from_slice(&self.mid_hist);
        m.extend_from_slice(mid);
        let mut s = Vec::with_capacity(n2 + 1);
        s.push(self.side_hist);
        s.extend_from_slice(side);
        let (pw0, pw1) = (
            self.prev_w[0] as f32 / 8192.0,
            self.prev_w[1] as f32 / 8192.0,
        );
        let (dw0, dw1) = (
            (w[0] - self.prev_w[0]) as f32 / (8192.0 * n1 as f32),
            (w[1] - self.prev_w[1]) as f32 / (8192.0 * n1 as f32),
        );
        let mut left = vec![0.0f32; n2];
        let mut right = vec![0.0f32; n2];
        for i in 0..n2 {
            let t = i.min(n1) as f32;
            let w0 = pw0 + t * dw0;
            let w1 = pw1 + t * dw1;
            let p0 = (m[i] + 2.0 * m[i + 1] + m[i + 2]) / 4.0;
            let mi = m[i + 1];
            let si = s[i];
            left[i] = ((1.0 + w1) * mi + si + w0 * p0).clamp(-1.0, 1.0);
            right[i] = ((1.0 - w1) * mi - si - w0 * p0).clamp(-1.0, 1.0);
        }
        self.mid_hist = [m[n2], m[n2 + 1]];
        self.side_hist = s[n2];
        self.prev_w = w;
        (left, right)
    }

    /// The mono path: the mid channel with the same one-sample delay.
    pub fn delay_mono(&mut self, mid: &[f32]) -> Vec<f32> {
        let n = mid.len();
        let mut out = Vec::with_capacity(n);
        out.push(self.mid_hist[1]);
        out.extend_from_slice(&mid[..n - 1]);
        self.mid_hist = [
            if n >= 2 { mid[n - 2] } else { self.mid_hist[1] },
            mid[n - 1],
        ];
        self.side_hist = 0.0;
        out
    }
}
