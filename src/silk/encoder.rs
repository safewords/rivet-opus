//! The SILK encoder (RFC 6716 §5.2).
//!
//! An original design that follows the RFC's description in outline:
//! - voice activity from a tracked noise floor (§5.2.3.1);
//! - pitch from the normalized autocorrelation of the LPC residual, refined
//!   per subframe through the contour codebooks (§5.2.3.2);
//! - LPC by the autocorrelation method on a windowed block with lookahead,
//!   converted to normalized LSFs and quantized by a two-stage search that
//!   reconstructs exactly as the decoder does (§5.2.3.5);
//! - 5-tap LTP fitted per subframe and quantized to the codebook with the
//!   lowest weighted rate-distortion cost (§5.2.3.6);
//! - a closed-loop scalar quantizer that runs the decoder's own synthesis
//!   (§5.2.3.8), so encoder and decoder never drift;
//! - a rate loop that scales the subframe gains until the frame fits its
//!   bit budget (§5.2.3.9);
//! - stereo as adaptive mid/side with the decoder's prediction (§5.2.2);
//! - LBRR (in-band FEC) frames at a coarser quantization (§4.2.5).
//!
//! Noise shaping is limited to what the closed loop gives (white output
//! noise); the delayed-decision quantizer and the pre-filter are not
//! implemented.

use super::decoder::{
    ChannelState, CondCoding, FrameParams, contour_icdf, contour_offsets, gain_q16, ltp_taps,
};
use super::nlsf;
use super::tables::*;
use super::{FrameIndices, SignalType};
use crate::range::RangeEncoder;

/// Cost in bits of symbol `s` of an inverse-CDF table (ft 256).
fn sym_bits(icdf: &[u8], s: usize) -> f32 {
    let hi = if s == 0 { 256.0 } else { f32::from(icdf[s - 1]) };
    let p = (hi - f32::from(icdf[s])) / 256.0;
    if p <= 0.0 { 30.0 } else { -p.log2() }
}

fn pitch_params(fs_khz: usize) -> (i32, i32, i32, &'static [u8]) {
    match fs_khz {
        8 => (4, 16, 144, &PITCH_LOW_NB_ICDF),
        12 => (6, 24, 216, &PITCH_LOW_MB_ICDF),
        _ => (8, 32, 288, &PITCH_LOW_WB_ICDF),
    }
}

/// Encodes one frame's parameters (the inverse of `decode_indices`).
pub fn encode_indices(
    enc: &mut RangeEncoder,
    ix: &FrameIndices,
    fs_khz: usize,
    nb_subfr: usize,
    active: bool,
    cond: CondCoding,
    prev_signal_type: SignalType,
    prev_lag: i32,
) {
    let ftype = ix.signal_type as usize * 2 + ix.qoff;
    if active {
        enc.icdf(ftype - 2, &FRAME_TYPE_ACTIVE_ICDF, 8);
    } else {
        enc.icdf(ftype, &FRAME_TYPE_INACTIVE_ICDF, 8);
    }
    for k in 0..nb_subfr {
        if k == 0 && cond != CondCoding::Conditional {
            enc.icdf((ix.gains[0] >> 3) as usize, &GAIN_MSB_ICDF[ix.signal_type as usize], 8);
            enc.icdf((ix.gains[0] & 7) as usize, &GAIN_LSB_ICDF, 8);
        } else {
            enc.icdf(ix.gains[k] as usize, &GAIN_DELTA_ICDF, 8);
        }
    }
    let wb = fs_khz == 16;
    let s1 = usize::from(wb) * 2 + usize::from(ix.signal_type == SignalType::Voiced);
    enc.icdf(ix.nlsf_i1, &NLSF_STAGE1_ICDF[s1], 8);
    for k in 0..nlsf::order(wb) {
        let v = ix.nlsf_i2[k];
        enc.icdf((v.clamp(-4, 4) + 4) as usize, nlsf::stage2_icdf(wb, ix.nlsf_i1, k), 8);
        if v <= -4 {
            enc.icdf((-4 - v) as usize, &NLSF_EXT_ICDF, 8);
        } else if v >= 4 {
            enc.icdf((v - 4) as usize, &NLSF_EXT_ICDF, 8);
        }
    }
    if nb_subfr == 4 {
        enc.icdf(ix.interp_q2 as usize, &NLSF_INTERP_ICDF, 8);
    }
    if ix.signal_type == SignalType::Voiced {
        let (scale, min_lag, _, low_icdf) = pitch_params(fs_khz);
        let mut absolute = true;
        if cond == CondCoding::Conditional && prev_signal_type == SignalType::Voiced {
            let delta = ix.lag - prev_lag + 9;
            if (1..=20).contains(&delta) {
                enc.icdf(delta as usize, &PITCH_DELTA_ICDF, 8);
                absolute = false;
            } else {
                enc.icdf(0, &PITCH_DELTA_ICDF, 8);
            }
        }
        if absolute {
            let v = ix.lag - min_lag;
            enc.icdf((v / scale) as usize, &PITCH_HIGH_ICDF, 8);
            enc.icdf((v % scale) as usize, low_icdf, 8);
        }
        enc.icdf(ix.contour, contour_icdf(fs_khz, nb_subfr), 8);
        enc.icdf(ix.periodicity, &PERIODICITY_ICDF, 8);
        for k in 0..nb_subfr {
            enc.icdf(ix.ltp[k], LTP_FILTER_ICDF[ix.periodicity], 8);
        }
        if cond == CondCoding::Independent {
            enc.icdf(ix.ltp_scale, &LTP_SCALE_ICDF, 8);
        }
    }
    enc.icdf(ix.seed as usize, &SEED_ICDF, 8);
}

fn shell_split(enc: &mut RangeEncoder, left: i32, total: i32, table: &[&[u8]]) {
    if total > 0 {
        enc.icdf(left as usize, table[total as usize - 1], 8);
    }
}

fn shell_encode(enc: &mut RangeEncoder, m: &[i32]) {
    let sum = |a: usize, b: usize| m[a..b].iter().sum::<i32>();
    shell_split(enc, sum(0, 8), sum(0, 16), &SHELL16_ICDF);
    for half in [0usize, 8] {
        shell_split(enc, sum(half, half + 4), sum(half, half + 8), &SHELL8_ICDF);
        for q in [half, half + 4] {
            shell_split(enc, sum(q, q + 2), sum(q, q + 4), &SHELL4_ICDF);
            for r in [q, q + 2] {
                shell_split(enc, m[r], m[r] + m[r + 1], &SHELL2_ICDF);
            }
        }
    }
}

/// Per-block decomposition: LSB count and the pulse count of the top part.
fn block_shape(q: &[i32]) -> (usize, i32) {
    let mut s = 0usize;
    loop {
        let count: i32 = q.iter().map(|v| v.abs() >> s).sum();
        if count <= 16 || s == 10 {
            return (s, count.min(16));
        }
        s += 1;
    }
}

fn level_after(rate_level: usize, t: usize) -> usize {
    if t == 0 { rate_level } else if t == 10 { 10 } else { 9 }
}

/// Encodes the excitation pulses (the inverse of `decode_pulses`); `q` has
/// a whole number of 16-sample blocks.
pub fn encode_pulses(enc: &mut RangeEncoder, q: &[i32], signal_type: SignalType, qoff: usize) {
    let blocks = q.len() / 16;
    let shapes: Vec<(usize, i32)> = (0..blocks).map(|b| block_shape(&q[b * 16..b * 16 + 16])).collect();
    let rl_icdf = &RATE_LEVEL_ICDF[usize::from(signal_type == SignalType::Voiced)];
    let mut best = (f32::MAX, 0usize);
    for level in 0..9 {
        let mut bits = sym_bits(rl_icdf, level);
        for &(s, count) in &shapes {
            for t in 0..s {
                bits += sym_bits(&PULSE_COUNT_ICDF[level_after(level, t)], 17);
            }
            bits += sym_bits(&PULSE_COUNT_ICDF[level_after(level, s)], count as usize);
        }
        if bits < best.0 {
            best = (bits, level);
        }
    }
    let level = best.1;
    enc.icdf(level, rl_icdf, 8);
    for &(s, count) in &shapes {
        for t in 0..s {
            enc.icdf(17, &PULSE_COUNT_ICDF[level_after(level, t)], 8);
        }
        enc.icdf(count as usize, &PULSE_COUNT_ICDF[level_after(level, s)], 8);
    }
    for (b, &(s, count)) in shapes.iter().enumerate() {
        if count > 0 {
            let m: Vec<i32> = q[b * 16..b * 16 + 16].iter().map(|v| v.abs() >> s).collect();
            shell_encode(enc, &m);
        }
    }
    for (b, &(s, _)) in shapes.iter().enumerate() {
        if s > 0 {
            for &v in &q[b * 16..b * 16 + 16] {
                for bit in (0..s).rev() {
                    enc.icdf(((v.abs() >> bit) & 1) as usize, &LSB_ICDF, 8);
                }
            }
        }
    }
    let group = signal_type as usize * 2 + qoff;
    for (b, &(_, count)) in shapes.iter().enumerate() {
        let icdf = &SIGN_ICDF[group][count.min(6) as usize];
        for &v in &q[b * 16..b * 16 + 16] {
            if v != 0 {
                enc.icdf(usize::from(v > 0), icdf, 8);
            }
        }
    }
}

/// Levinson-Durbin on an autocorrelation; returns `a[1..=d]` with
/// `A(z) = 1 - sum a_k z^-k`, and the residual energy.
fn levinson(r: &[f64], d: usize) -> (Vec<f64>, f64) {
    let mut a = vec![0.0f64; d + 1];
    let mut err = r[0];
    for i in 1..=d {
        if err <= 1e-12 {
            break;
        }
        let mut acc = r[i];
        for j in 1..i {
            acc -= a[j] * r[i - j];
        }
        let k = (acc / err).clamp(-0.9999, 0.9999);
        let prev = a.clone();
        a[i] = k;
        for j in 1..i {
            a[j] = prev[j] - k * prev[i - j];
        }
        err *= 1.0 - k * k;
    }
    (a[1..].to_vec(), err)
}

/// LPC analysis of `x` (windowed internally).
fn lpc_analysis(x: &[f32], d: usize) -> Vec<f64> {
    let n = x.len();
    let w: Vec<f64> = (0..n)
        .map(|i| {
            let t = (i as f64 + 0.5) / n as f64;
            f64::from(x[i]) * (std::f64::consts::PI * t).sin()
        })
        .collect();
    let mut r = vec![0.0f64; d + 1];
    for (k, rk) in r.iter_mut().enumerate() {
        *rk = (k..n).map(|i| w[i] * w[i - k]).sum();
    }
    // White noise correction and a lag window.
    r[0] *= 1.0 + 1e-4;
    r[0] += 1e-9;
    for (k, rk) in r.iter_mut().enumerate().skip(1) {
        let g = 2.0 * std::f64::consts::PI * 60.0 * k as f64 / 16000.0;
        *rk *= (-0.5 * g * g).exp();
    }
    let (mut a, _) = levinson(&r, d);
    let mut g = 1.0;
    for v in a.iter_mut() {
        g *= 0.995;
        *v *= g;
    }
    a
}

/// Normalized LSFs (Q15) of an LPC filter, or `None` if the roots could not
/// be separated.
fn lpc_to_nlsf(a: &[f64]) -> Option<Vec<i32>> {
    let d = a.len();
    let mut c = vec![0.0f64; d + 2];
    c[0] = 1.0;
    for k in 1..=d {
        c[k] = -a[k - 1];
    }
    let p: Vec<f64> = (0..=d + 1).map(|k| c[k] + c[d + 1 - k]).collect();
    let q: Vec<f64> = (0..=d + 1).map(|k| c[k] - c[d + 1 - k]).collect();
    let half = (d + 1) as f64 / 2.0;
    let fp = |w: f64| -> f64 { p.iter().enumerate().map(|(k, v)| v * (w * (half - k as f64)).cos()).sum() };
    let fq = |w: f64| -> f64 { q.iter().enumerate().map(|(k, v)| v * (w * (half - k as f64)).sin()).sum() };
    let grid = 2048;
    let mut roots: Vec<(f64, bool)> = Vec::with_capacity(d);
    for (is_p, f) in [(true, &fp as &dyn Fn(f64) -> f64), (false, &fq as &dyn Fn(f64) -> f64)] {
        let mut prev_w = 1e-6;
        let mut prev_v = f(prev_w);
        for j in 1..grid {
            let w = std::f64::consts::PI * j as f64 / grid as f64;
            let v = f(w);
            if prev_v == 0.0 || prev_v.signum() != v.signum() {
                let (mut lo, mut hi, mut flo) = (prev_w, w, prev_v);
                for _ in 0..30 {
                    let mid = 0.5 * (lo + hi);
                    let fm = f(mid);
                    if fm.signum() == flo.signum() {
                        lo = mid;
                        flo = fm;
                    } else {
                        hi = mid;
                    }
                }
                roots.push((0.5 * (lo + hi), is_p));
            }
            prev_w = w;
            prev_v = v;
        }
    }
    roots.sort_by(|x, y| x.0.total_cmp(&y.0));
    if roots.len() != d {
        return None;
    }
    for (i, r) in roots.iter().enumerate() {
        if r.1 != (i % 2 == 0) {
            return None;
        }
    }
    Some(roots.iter().map(|r| ((r.0 / std::f64::consts::PI) * 32768.0).round().clamp(1.0, 32767.0) as i32).collect())
}

/// The stage-2 dequantized step for index `i`.
fn dequant_step(i: i32, qstep: i32) -> i32 {
    ((((i << 10) - i.signum() * 102) * qstep) >> 16) as i32
}

/// Quantizes normalized LSFs; returns (I1, I2) whose decoded LSFs are
/// closest (IHMW-weighted, with a rate term).
fn quantize_nlsf(target: &[i32], wb: bool, voiced: bool) -> (usize, [i32; 16]) {
    let d = nlsf::order(wb);
    let qstep = nlsf::qstep(wb);
    // Weights of the input vector (Laroia), as floats.
    let mut w_in = vec![0.0f64; d];
    for k in 0..d {
        let prev = if k == 0 { 0 } else { target[k - 1] };
        let next = if k + 1 == d { 32768 } else { target[k + 1] };
        w_in[k] = 1.0 / f64::from((target[k] - prev).max(64)) + 1.0 / f64::from((next - target[k]).max(64));
    }
    let s1 = usize::from(wb) * 2 + usize::from(voiced);
    let mut stage1: Vec<(f64, usize)> = (0..32)
        .map(|i1| {
            let cb = nlsf::cb1(wb, i1);
            let e: f64 = (0..d).map(|k| w_in[k] * (f64::from(target[k] - (cb[k] << 7))).powi(2)).sum();
            (e, i1)
        })
        .collect();
    stage1.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut best = (f64::MAX, 0usize, [0i32; 16]);
    for &(_, i1) in stage1.iter().take(6) {
        let cb = nlsf::cb1(wb, i1);
        let wq9 = nlsf::weights_q9(cb, d);
        let mut i2 = [0i32; 16];
        let mut res_q = [0i32; 16];
        let mut rate = sym_bits(&NLSF_STAGE1_ICDF[s1], i1) as f64;
        for k in (0..d).rev() {
            let want = (f64::from(target[k] - (cb[k] << 7)) * f64::from(wq9[k]) / 16384.0) as i32;
            let pred = if k + 1 < d { (res_q[k + 1] * nlsf::pred_q8(wb, i1, k)) >> 8 } else { 0 };
            let icdf = nlsf::stage2_icdf(wb, i1, k);
            let mut bi = (f64::MAX, 0i32);
            for i in -10..=10 {
                let v = pred + dequant_step(i, qstep);
                let e = f64::from(want - v).powi(2);
                let mut bits = sym_bits(icdf, (i.clamp(-4, 4) + 4) as usize);
                if i.abs() >= 4 {
                    bits += sym_bits(&NLSF_EXT_ICDF, (i.abs() - 4) as usize);
                }
                let cost = e + 400.0 * f64::from(bits);
                if cost < bi.0 {
                    bi = (cost, i);
                }
            }
            i2[k] = bi.1;
            res_q[k] = pred + dequant_step(bi.1, qstep);
            let mut bits = sym_bits(icdf, (bi.1.clamp(-4, 4) + 4) as usize);
            if bi.1.abs() >= 4 {
                bits += sym_bits(&NLSF_EXT_ICDF, (bi.1.abs() - 4) as usize);
            }
            rate += f64::from(bits);
        }
        let mut rec = nlsf::reconstruct(wb, i1, &i2);
        nlsf::stabilize(&mut rec[..d], wb);
        let e: f64 = (0..d).map(|k| w_in[k] * f64::from(target[k] - rec[k]).powi(2)).sum();
        let cost = e + 0.02 * rate;
        if cost < best.0 {
            best = (cost, i1, i2);
        }
    }
    (best.1, best.2)
}

/// The gain index whose Q16 gain is closest (in log) to `want`.
fn gain_index(want: f64) -> i32 {
    let mut best = (f64::MAX, 0);
    for i in 0..64 {
        let g = f64::from(gain_q16(i));
        let e = (g.ln() - want.max(1.0).ln()).abs();
        if e < best.0 {
            best = (e, i);
        }
    }
    best.1
}

/// Coded gain indices (absolute or delta) that land the decoder's gain
/// indices on `targets`, from the decoder's previous index `prev`.
fn code_gains(targets: &[i32], prev: i32, cond: CondCoding) -> [i32; 4] {
    let mut out = [0i32; 4];
    let mut prev = prev;
    for (k, &t) in targets.iter().enumerate() {
        if k == 0 && cond != CondCoding::Conditional {
            out[0] = t;
            prev = t.max(prev - 16);
        } else {
            let mut best = (i32::MAX, 0);
            for d in 0..=40 {
                let v = (2 * d - 16).max(prev + d - 4).clamp(0, 63);
                let e = (v - t).abs();
                if e < best.0 {
                    best = (e, d);
                }
            }
            out[k] = best.1;
            prev = (2 * best.1 - 16).max(prev + best.1 - 4).clamp(0, 63);
        }
    }
    out
}

/// The analysis of one frame, independent of the rate.
#[derive(Clone)]
struct Analysis {
    active: bool,
    voiced: bool,
    nlsf_i1: usize,
    nlsf_i2: [i32; 16],
    lag: i32,
    contour: usize,
    periodicity: usize,
    ltp: [usize; 4],
    /// Residual RMS per subframe (output units), after LPC and LTP.
    sigma: [f64; 4],
    /// The 16 input samples before the frame, then the frame.
    ext: Vec<f32>,
    /// RMS of the frame.
    level: f64,
}

/// Normalized cross-correlation of `r[n]` and `r[n - lag]` over `range`.
fn ncorr(r: &[f32], start: usize, len: usize, lag: usize) -> f64 {
    let (mut xy, mut xx, mut yy) = (0.0f64, 0.0f64, 0.0f64);
    for n in start..start + len {
        let x = f64::from(r[n]);
        let y = f64::from(r[n - lag]);
        xy += x * y;
        xx += x * x;
        yy += y * y;
    }
    if xx <= 1e-12 || yy <= 1e-12 { 0.0 } else { xy / (xx * yy).sqrt() }
}

/// One channel of the SILK encoder.
#[derive(Clone)]
struct ChannelEncoder {
    state: ChannelState,
    /// Input history at the internal rate (for pitch and LPC analysis).
    hist: Vec<f32>,
    noise_floor: f64,
    seed: u32,
    prev_voiced: bool,
}

const HIST: usize = 640;

impl ChannelEncoder {
    fn new(fs_khz: usize) -> Self {
        Self { state: ChannelState::new(fs_khz), hist: vec![0.0; HIST], noise_floor: 1e-6, seed: 0, prev_voiced: false }
    }

    /// Analyses `frame` (with `ahead` samples of lookahead after it).
    fn analyse(&mut self, frame: &[f32], ahead: &[f32], nb_subfr: usize) -> Analysis {
        let fs = self.state.fs_khz;
        let d = self.state.lpc_order;
        let n = frame.len();
        let sub = n / nb_subfr;
        let wb = fs == 16;
        // Voice activity.
        let energy = frame.iter().map(|v| f64::from(*v) * f64::from(*v)).sum::<f64>() / n as f64;
        let active = energy > 3.0 * self.noise_floor && energy > 1e-8;
        self.noise_floor = if energy < self.noise_floor { 0.7 * self.noise_floor + 0.3 * energy } else { self.noise_floor * 1.02 };
        self.noise_floor = self.noise_floor.max(1e-9);
        // LPC on a block: 5 ms of history, the frame, the lookahead.
        let pre = 5 * fs;
        let mut block: Vec<f32> = self.hist[HIST - pre..].to_vec();
        block.extend_from_slice(frame);
        block.extend_from_slice(ahead);
        let a = lpc_analysis(&block, d);
        let nlsf_target = lpc_to_nlsf(&a).unwrap_or_else(|| {
            if self.state.first_frame_after_reset {
                (0..d).map(|k| ((k as i32 + 1) * 32768) / (d as i32 + 1)).collect()
            } else {
                self.state.prev_nlsf[..d].to_vec()
            }
        });
        let mut nlsf_target = nlsf_target;
        nlsf::stabilize(&mut nlsf_target, wb);
        // Pitch from the residual of history + frame.
        let mut sig: Vec<f32> = self.hist.clone();
        sig.extend_from_slice(frame);
        let mut res = vec![0.0f32; sig.len()];
        for i in d..sig.len() {
            let mut v = f64::from(sig[i]);
            for k in 0..d {
                v -= a[k] * f64::from(sig[i - k - 1]);
            }
            res[i] = v as f32;
        }
        let (_, min_lag, max_lag, _) = pitch_params(fs);
        let start = HIST;
        let mut best = (0.0f64, min_lag);
        for lag in min_lag..=max_lag {
            let c = ncorr(&res, start, n, lag as usize);
            // A slight preference for short lags avoids pitch multiples.
            let score = c * (1.0 - 0.15 * (f64::from(lag) / f64::from(min_lag)).log2() / 4.0);
            if score > best.0 {
                best = (score, lag);
            }
        }
        let threshold = if self.prev_voiced { 0.35 } else { 0.45 };
        let voiced = active && best.0 > threshold;
        let mut lag = best.1;
        let mut contour = 0;
        let mut ltp = [0usize; 4];
        let mut periodicity = 0;
        let mut lags = [0i32; 4];
        if voiced {
            // Refine: primary lag near the best one and a contour.
            let cb_count = contour_icdf(fs, nb_subfr).len();
            let mut bestc = (f64::MIN, lag, 0usize);
            for p in (lag - 2).max(min_lag)..=(lag + 2).min(max_lag) {
                for ci in 0..cb_count {
                    let off = contour_offsets(fs, nb_subfr, ci);
                    let mut s = 0.0;
                    for k in 0..nb_subfr {
                        let l = (p + off[k]).clamp(min_lag, max_lag) as usize;
                        s += ncorr(&res, start + k * sub, sub, l);
                    }
                    if s > bestc.0 {
                        bestc = (s, p, ci);
                    }
                }
            }
            lag = bestc.1;
            contour = bestc.2;
            let off = contour_offsets(fs, nb_subfr, contour);
            for k in 0..nb_subfr {
                lags[k] = (lag + off[k]).clamp(min_lag, max_lag);
            }
            // LTP: least squares per subframe on the residual.
            let mut fits = Vec::with_capacity(nb_subfr);
            for k in 0..nb_subfr {
                let l = lags[k] as usize;
                let s0 = start + k * sub;
                let mut wm = [[0.0f64; 5]; 5];
                let mut wv = [0.0f64; 5];
                let mut e0 = 1e-9;
                for nn in s0..s0 + sub {
                    let t = f64::from(res[nn]);
                    e0 += t * t;
                    let v: [f64; 5] = std::array::from_fn(|j| f64::from(res[nn + 2 - j - l]));
                    for i in 0..5 {
                        wv[i] += v[i] * t;
                        for j in 0..5 {
                            wm[i][j] += v[i] * v[j];
                        }
                    }
                }
                fits.push((wm, wv, e0));
            }
            let mut best_p = (f64::MAX, 0usize, [0usize; 4]);
            for p in 0..3 {
                let count = LTP_FILTER_ICDF[p].len();
                let mut total = sym_bits(&PERIODICITY_ICDF, p) as f64;
                let mut idx = [0usize; 4];
                for (k, (wm, wv, e0)) in fits.iter().enumerate() {
                    let mut bi = (f64::MAX, 0usize);
                    for i in 0..count {
                        let taps = ltp_taps(p, i);
                        let b: [f64; 5] = std::array::from_fn(|j| f64::from(taps[j]) / 128.0);
                        // Residual energy after prediction with b.
                        let mut e = *e0;
                        for x in 0..5 {
                            e -= 2.0 * b[x] * wv[x];
                            for y in 0..5 {
                                e += b[x] * b[y] * wm[x][y];
                            }
                        }
                        let e = e.max(1e-3 * e0);
                        let bits = f64::from(sym_bits(LTP_FILTER_ICDF[p], i));
                        let cost = 0.5 * (sub as f64) * (e / e0).log2() + bits;
                        if cost < bi.0 {
                            bi = (cost, i);
                        }
                    }
                    idx[k] = bi.1;
                    total += bi.0;
                }
                if total < best_p.0 {
                    best_p = (total, p, idx);
                }
            }
            periodicity = best_p.1;
            ltp = best_p.2;
        }
        // Quantize the LSFs and find the residual level with the quantized
        // filter.
        let (i1, i2) = quantize_nlsf(&nlsf_target, wb, voiced);
        let mut rec = nlsf::reconstruct(wb, i1, &i2);
        nlsf::stabilize(&mut rec[..d], wb);
        let aq: Vec<f64> = nlsf::nlsf_to_lpc(&rec[..d], wb).iter().map(|&v| f64::from(v) / 4096.0).collect();
        let mut sigma = [0.0f64; 4];
        for k in 0..nb_subfr {
            let s0 = start + k * sub;
            let mut e = 0.0f64;
            for nn in s0..s0 + sub {
                let mut v = f64::from(sig[nn]);
                for j in 0..d {
                    v -= aq[j] * f64::from(sig[nn - j - 1]);
                }
                e += v * v;
            }
            // The LPC residual level; voiced frames get part of it from the
            // LTP, so their step can be finer for the same rate.
            let scale = if voiced { 0.6 } else { 1.0 };
            let level = (frame[k * sub..(k + 1) * sub].iter().map(|v| f64::from(*v) * f64::from(*v)).sum::<f64>() / sub as f64).sqrt();
            sigma[k] = scale * (e / sub as f64).sqrt().max(0.1 * level).max(1e-6);
        }
        // A subframe much quieter than its neighbours still has to absorb
        // their filters' ringing.
        let top = sigma[..nb_subfr].iter().cloned().fold(0.0f64, f64::max);
        for v in sigma[..nb_subfr].iter_mut() {
            *v = v.max(0.2 * top);
        }
        let mut ext = self.hist[HIST - 16..].to_vec();
        ext.extend_from_slice(frame);
        // History for the next frame.
        let mut h = self.hist[n.min(HIST)..].to_vec();
        h.extend_from_slice(&frame[frame.len().saturating_sub(HIST)..]);
        self.hist = h[h.len() - HIST..].to_vec();
        self.prev_voiced = voiced;
        let level = energy.sqrt();
        Analysis { active, voiced, nlsf_i1: i1, nlsf_i2: i2, lag, contour, periodicity, ltp, sigma, ext, level }
    }

    /// One trial encoding of a frame at gain multiplier `mult`, on clones
    /// of the state and coder; returns them with the bits used.
    #[allow(clippy::too_many_arguments)]
    fn trial(
        &self,
        an: &Analysis,
        frame: &[f32],
        nb_subfr: usize,
        mult: f64,
        cond: CondCoding,
        active: bool,
        enc: &RangeEncoder,
        ltp_scale: usize,
    ) -> (ChannelState, RangeEncoder, i32, u32) {
        // A subframe whose pulses hit the representable limit had a step
        // too fine for its (closed-loop) residual: coarsen it and retry.
        let mut boost = [1.0f64; 4];
        for _ in 0..4 {
            let (st, e2, bits, seed, sat) = self.trial_once(an, frame, nb_subfr, mult, &boost, cond, active, enc, ltp_scale);
            if !sat.iter().any(|&s| s) {
                return (st, e2, bits, seed);
            }
            for k in 0..nb_subfr {
                if sat[k] {
                    boost[k] *= 16.0;
                }
            }
        }
        let (st, e2, bits, seed, _) = self.trial_once(an, frame, nb_subfr, mult, &boost, cond, active, enc, ltp_scale);
        (st, e2, bits, seed)
    }

    #[allow(clippy::too_many_arguments)]
    fn trial_once(
        &self,
        an: &Analysis,
        frame: &[f32],
        nb_subfr: usize,
        mult: f64,
        boost: &[f64; 4],
        cond: CondCoding,
        active: bool,
        enc: &RangeEncoder,
        ltp_scale: usize,
    ) -> (ChannelState, RangeEncoder, i32, u32, [bool; 4]) {
        let mut st = self.state.clone();
        let mut e2 = enc.clone();
        let start_bits = e2.tell_frac();
        let fs = st.fs_khz;
        let signal_type = if !active {
            SignalType::Inactive
        } else if an.voiced {
            SignalType::Voiced
        } else {
            SignalType::Unvoiced
        };
        // Gain targets: the step size tracks the residual level.
        // A step much above the signal itself only adds noise.
        let cap = if mult > 64.0 { 63 } else { gain_index(an.level.max(1e-5) * 2.0 * 2_147_483_648.0) };
        let targets: Vec<i32> =
            (0..nb_subfr).map(|k| gain_index(an.sigma[k] * mult * boost[k] * 2_147_483_648.0).min(cap)).collect();
        let mut sat = [false; 4];
        let gains = code_gains(&targets, st.last_gain_index, cond);
        let ix = FrameIndices {
            signal_type,
            qoff: 0,
            gains,
            nlsf_i1: an.nlsf_i1,
            nlsf_i2: an.nlsf_i2,
            interp_q2: 4,
            lag: an.lag,
            contour: an.contour,
            periodicity: an.periodicity,
            ltp: an.ltp,
            ltp_scale: if signal_type == SignalType::Voiced { ltp_scale } else { 0 },
            seed: self.seed & 3,
        };
        let prev_type = st.prev_signal_type;
        let prev_lag = st.prev_lag;
        let params: FrameParams = st.dequantize(&ix, nb_subfr, cond);
        // The LPC residual of the input with the decoder's coefficients:
        // quantizing it (rather than the closed-loop error) shapes the noise
        // like the spectrum and keeps the LPC loop from feeding the noise
        // back into the target.
        let n = frame.len();
        let d = st.lpc_order;
        let sub = n / nb_subfr;
        let mut r_in = vec![0.0f32; n];
        for (i, r) in r_in.iter_mut().enumerate() {
            let a = &params.a_q12[usize::from(!(i / sub < 2 && params.interp))];
            let mut v = an.ext[16 + i];
            for k in 0..d {
                v -= an.ext[16 + i - k - 1] * a[k] as f32 / 4096.0;
            }
            *r = v;
        }
        // Noise feedback: with output error e(n) = x(n) - y(n), quantizing
        // r(n) + sum_k a_k (1 - g^k) e(n - k) shapes the output noise by
        // 1 / A(z / g): white (plain closed loop) at g = 0, the spectral
        // envelope (open loop) at g = 1. Coarse steps shape more.
        let g = (0.7 + 0.08 * mult.log2().max(0.0)).min(0.94) as f32;
        let coef: Vec<[f32; 16]> = (0..2)
            .map(|h| {
                let mut c = [0.0f32; 16];
                let mut gk = 1.0f32;
                for k in 0..d {
                    gk *= g;
                    c[k] = params.a_q12[h][k] as f32 / 4096.0 * (1.0 - gk);
                }
                c
            })
            .collect();
        let mut err_hist = vec![0.0f32; n + 16];
        let offset = QUANT_OFFSETS_Q23[signal_type as usize][0];
        let mut seed = ix.seed;
        let mut q = vec![0i32; n.div_ceil(16) * 16];
        // Coarser steps also get a wider dead zone (fewer pulses).
        let dead_zone = (0.1 + 0.15 * mult.log2().max(0.0)).min(0.6) as f32;
        {
            let mut choose = |i: usize, ltp: f32, pred: f32, gain: f32| -> i32 {
                let c = &coef[usize::from(!(i / sub < 2 && params.interp))];
                let mut fb = 0.0f32;
                for k in 0..d {
                    fb += c[k] * err_hist[16 + i - k - 1];
                }
                let want_res = (r_in[i] + fb) * 65536.0 / gain - ltp;
                let want = want_res * 8_388_608.0;
                seed = seed.wrapping_mul(196_314_165).wrapping_add(907_633_515);
                let flip = seed & 0x8000_0000 != 0;
                let u = if flip { -want } else { want };
                // Find q with (q << 8) - sign(q) * 20 + offset closest to u.
                let x = (u - offset as f32) / 256.0;
                let mut qi = (x.abs() - dead_zone + 0.5).floor().max(0.0) as i32 * if x < 0.0 { -1 } else { 1 };
                let rec = |qq: i32| (qq << 8) - qq.signum() * 20 + offset;
                for cand in [qi - 1, qi + 1] {
                    let cost = |qq: i32| (rec(qq) as f32 - u).abs() + 64.0 * (qq.abs() as f32);
                    let better = cost(cand) < cost(qi);
                    if better {
                        qi = cand;
                    }
                }
                // Any block of such values fits 10 LSB levels.
                qi = qi.clamp(-1023, 1023);
                if qi.abs() >= 1023 {
                    sat[(i * nb_subfr / n).min(3)] = true;
                }
                q[i] = qi;
                let mut e = rec(qi);
                if flip {
                    e = -e;
                }
                seed = seed.wrapping_add(qi as u32);
                let y = gain / 65536.0 * (e as f32 / 8_388_608.0 + ltp) + pred;
                err_hist[16 + i] = frame[i] - y;
                e
            };
            st.synthesize_with(&params, nb_subfr, &mut choose);
        }
        encode_indices(&mut e2, &ix, fs, nb_subfr, active, cond, prev_type, prev_lag);
        encode_pulses(&mut e2, &q, signal_type, 0);
        let bits = e2.tell_frac() - start_bits;
        (st, e2, bits, ix.seed, sat)
    }
}

/// The SILK encoder of one Opus stream.
pub struct SilkEncoder {
    channels: usize,
    fs_khz: usize,
    ch: Vec<ChannelEncoder>,
    /// Stereo: previous quantized weights (Q13) and mid/side state.
    prev_w: [i32; 2],
    side_prev_uncoded: bool,
    /// The last sample of mid of the previous frame (for the low-pass).
    mid_tail: [f32; 2],
    /// LBRR data waiting for the next packet: per frame, per channel.
    lbrr: Vec<[Option<LbrrFrame>; 2]>,
    lbrr_w: Vec<[i32; 2]>,
    lbrr_mid_only: Vec<bool>,
    /// Hybrid frames must not exceed their SILK budget (the CELT layer
    /// follows); SILK-only frames may grow rather than turn to noise.
    pub hard_limit: bool,
}

#[derive(Clone)]
struct LbrrFrame {
    ix: FrameIndices,
    q: Vec<i32>,
    cond: CondCoding,
    prev_type: SignalType,
    prev_lag: i32,
}

/// Quantizes a stereo weight to the codebook: (index pair, quantized Q13).
fn quantize_weight(w: f64) -> (usize, i32, i32) {
    let t = &STEREO_WEIGHTS_Q13;
    let mut best = (f64::MAX, 0usize, 0i32, 0i32);
    for wi in 0..15 {
        let step = ((t[wi + 1] - t[wi]) * 6554) >> 16;
        for s in 0..5 {
            let v = t[wi] + step * (2 * s + 1);
            let e = (f64::from(v) / 8192.0 - w).abs();
            if e < best.0 {
                best = (e, wi, s, v);
            }
        }
    }
    (best.1, best.2, best.3)
}

impl SilkEncoder {
    pub fn new(channels: usize) -> Self {
        Self {
            channels,
            fs_khz: 0,
            ch: Vec::new(),
            prev_w: [0, 0],
            side_prev_uncoded: false,
            mid_tail: [0.0; 2],
            lbrr: Vec::new(),
            lbrr_w: Vec::new(),
            lbrr_mid_only: Vec::new(),
            hard_limit: false,
        }
    }

    /// Back to the initial state.
    pub fn reset(&mut self) {
        *self = Self::new(self.channels);
    }

    /// The internal rate in kHz (0 before the first frame).
    pub fn fs_khz(&self) -> usize {
        self.fs_khz
    }

    fn prepare(&mut self, fs_khz: usize) {
        if fs_khz != self.fs_khz || self.ch.is_empty() {
            self.fs_khz = fs_khz;
            self.ch = (0..self.channels).map(|_| ChannelEncoder::new(fs_khz)).collect();
            self.prev_w = [0, 0];
            self.side_prev_uncoded = false;
            self.mid_tail = [0.0; 2];
            self.lbrr.clear();
        }
    }

    /// Encodes the SILK layer of one Opus frame into `enc`. `input` holds
    /// the frame's samples at the internal rate (interleaved, the stream's
    /// channels) followed by `ahead` samples of lookahead; `budget_bits`
    /// is what the layer may use; `fec` codes LBRR data for the next packet.
    #[allow(clippy::too_many_arguments)]
    pub fn encode(
        &mut self,
        enc: &mut RangeEncoder,
        input: &[f32],
        frame_len: usize,
        fs_khz: usize,
        frame_ms: usize,
        budget_bits: i32,
        fec: bool,
    ) {
        self.prepare(fs_khz);
        let c = self.channels;
        let total = input.len() / c;
        let ahead = total - frame_len;
        let nb_subfr = if frame_ms == 10 { 2 } else { 4 };
        let nf = (frame_ms / 20).max(1);
        let sub_len = frame_len / nf;
        // Mid/side (stereo) or the mono signal.
        let mut chans: Vec<Vec<f32>> = vec![Vec::with_capacity(total); c];
        if c == 1 {
            chans[0] = input.to_vec();
        } else {
            for i in 0..total {
                let (l, r) = (input[2 * i], input[2 * i + 1]);
                chans[0].push(0.5 * (l + r));
                chans[1].push(0.5 * (l - r));
            }
        }
        // Stereo prediction weights and the side residual per frame.
        let mut weights = vec![[0i32; 2]; nf];
        let mut w_idx = vec![(0usize, 0i32, 0usize, 0i32); nf];
        let mut side_active = vec![false; nf];
        if c == 2 {
            let n1 = 8 * fs_khz;
            for f in 0..nf {
                let s0 = f * sub_len;
                // Least squares for side ~ w0 * LP(mid) + w1 * mid.
                let (mut a11, mut a12, mut a22, mut b1, mut b2, mut es, mut em) = (1e-9f64, 0.0, 1e-9, 0.0, 0.0, 0.0, 1e-12);
                for n in s0..s0 + sub_len {
                    let mm = |k: isize| -> f64 {
                        if k < 0 { f64::from(self.mid_tail[(2 + k) as usize]) } else { f64::from(chans[0][k as usize]) }
                    };
                    let m = mm(n as isize);
                    let lp = 0.25 * (mm(n as isize - 1) + 2.0 * m + mm((n + 1).min(total - 1) as isize));
                    let s = f64::from(chans[1][n]);
                    a11 += lp * lp;
                    a12 += lp * m;
                    a22 += m * m;
                    b1 += lp * s;
                    b2 += m * s;
                    es += s * s;
                    em += m * m;
                }
                let det = a11 * a22 - a12 * a12;
                let (w0, w1) = if det.abs() > 1e-12 { ((b1 * a22 - b2 * a12) / det, (b2 * a11 - b1 * a12) / det) } else { (0.0, 0.0) };
                // Quantize w1 and w0 + w1 (Table 7 interpolation).
                let (wi1, i3, q1) = quantize_weight(w1.clamp(-1.6, 1.6));
                let (wi0, i1, q0s) = quantize_weight((w0 + f64::from(q1) / 8192.0).clamp(-1.6, 1.6));
                weights[f] = [q0s - q1, q1];
                w_idx[f] = (wi0, i1, wi1, i3);
                // Side residual with the decoder's interpolated weights.
                let prev = if f == 0 { self.prev_w } else { weights[f - 1] };
                let mut res_e = 0.0f64;
                for j in 0..sub_len {
                    let n = s0 + j;
                    let pos = (j + 1).min(n1) as f64;
                    let w0i = (f64::from(prev[0]) + pos * f64::from(weights[f][0] - prev[0]) / n1 as f64) / 8192.0;
                    let w1i = (f64::from(prev[1]) + pos * f64::from(weights[f][1] - prev[1]) / n1 as f64) / 8192.0;
                    let mm = |k: isize| -> f64 {
                        if k < 0 { f64::from(self.mid_tail[(2 + k) as usize]) } else { f64::from(chans[0][k as usize]) }
                    };
                    let m = mm(n as isize);
                    let lp = 0.25 * (mm(n as isize - 1) + 2.0 * m + mm((n + 1).min(total - 1) as isize));
                    let r = f64::from(chans[1][n]) - w1i * m - w0i * lp;
                    chans[1][n] = r as f32;
                    res_e += r * r;
                }
                let _ = es;
                side_active[f] = res_e > 1e-3 * em && res_e / sub_len as f64 > 1e-8;
            }
            self.mid_tail = [chans[0][frame_len - 2], chans[0][frame_len - 1]];
        }
        // Analysis.
        let mut analyses: Vec<Vec<Analysis>> = Vec::with_capacity(c);
        let mut frames: Vec<Vec<Vec<f32>>> = Vec::with_capacity(c);
        for ch in 0..c {
            let mut an = Vec::with_capacity(nf);
            let mut fr = Vec::with_capacity(nf);
            for f in 0..nf {
                let s0 = f * sub_len;
                let frame = chans[ch][s0..s0 + sub_len].to_vec();
                let la_end = (s0 + sub_len + ahead.min(5 * fs_khz)).min(total);
                let la = chans[ch][s0 + sub_len..la_end].to_vec();
                let a = self.ch[ch].analyse(&frame, &la, nb_subfr);
                an.push(a);
                fr.push(frame);
            }
            analyses.push(an);
            frames.push(fr);
        }
        let mut vad = [[false; 3]; 2];
        for ch in 0..c {
            for f in 0..nf {
                vad[ch][f] = analyses[ch][f].active && (ch == 0 || side_active[f]);
            }
        }
        // Header: VAD flags and LBRR flags.
        let have_lbrr = !self.lbrr.is_empty() && self.lbrr.len() == nf;
        let mut lbrr_flags = [[false; 3]; 2];
        if have_lbrr {
            for f in 0..nf {
                for ch in 0..c {
                    lbrr_flags[ch][f] = self.lbrr[f][ch].is_some();
                }
            }
        }
        for ch in 0..c {
            for f in 0..nf {
                enc.bit_logp(vad[ch][f], 1);
            }
            enc.bit_logp(lbrr_flags[ch].iter().any(|&x| x), 1);
        }
        for ch in 0..c {
            if lbrr_flags[ch].iter().any(|&x| x) && nf > 1 {
                let v: usize = (0..nf).map(|f| usize::from(lbrr_flags[ch][f]) << f).sum();
                let table: &[u8] = if nf == 2 { &LBRR_FLAGS_2_ICDF } else { &LBRR_FLAGS_3_ICDF };
                enc.icdf(v - 1, table, 8);
            }
        }
        // LBRR frames of the previous packet.
        if have_lbrr {
            let lbrr = std::mem::take(&mut self.lbrr);
            for f in 0..nf {
                for ch in 0..c {
                    let Some(lf) = &lbrr[f][ch] else { continue };
                    if ch == 0 && c == 2 {
                        Self::encode_weights(enc, self.lbrr_w[f]);
                        if !lbrr_flags[1][f] {
                            enc.icdf(usize::from(self.lbrr_mid_only[f]), &MID_ONLY_ICDF, 8);
                        }
                    }
                    encode_indices(enc, &lf.ix, fs_khz, nb_subfr, true, lf.cond, lf.prev_type, lf.prev_lag);
                    encode_pulses(enc, &lf.q, lf.ix.signal_type, lf.ix.qoff);
                }
            }
        }
        // Budget per frame and channel.
        let side_share = 0.35f64;
        let mut new_lbrr: Vec<[Option<LbrrFrame>; 2]> = Vec::new();
        let mut new_lbrr_w = Vec::new();
        let mut new_lbrr_mid_only = Vec::new();
        for f in 0..nf {
            let mut mid_only = false;
            let mut wq = [0i32; 2];
            if c == 2 {
                let (wi0, i1, wi1, i3) = w_idx[f];
                let n = (wi0 / 3) * 5 + wi1 / 3;
                let idx = [n as i32, (wi0 % 3) as i32, i1, (wi1 % 3) as i32, i3];
                Self::encode_weight_indices(enc, idx);
                wq = weights[f];
                if !vad[1][f] {
                    mid_only = true;
                    enc.icdf(1, &MID_ONLY_ICDF, 8);
                }
                self.prev_w = weights[f];
            }
            let frames_left = (nf - f) as i32;
            let avail = budget_bits - enc.tell();
            let frame_budget = (avail - 8) / frames_left;
            for ch in 0..c {
                if ch == 1 && mid_only {
                    continue;
                }
                let share = if c == 1 {
                    1.0
                } else if ch == 0 {
                    if mid_only { 1.0 } else { 1.0 - side_share }
                } else {
                    1.0
                };
                let target = if ch == 1 { budget_bits - enc.tell() - (frames_left - 1) * frame_budget - 8 } else { (f64::from(frame_budget) * share) as i32 };
                let mut cond = if f == 0 { CondCoding::Independent } else { CondCoding::Conditional };
                if ch == 1 && self.side_prev_uncoded {
                    self.ch[1].state.reset(fs_khz);
                    if f > 0 {
                        cond = CondCoding::IndependentNoLtpScaling;
                    }
                }
                let snapshot = self.ch[ch].state.clone();
                let an = analyses[ch][f].clone();
                let frame = &frames[ch][f];
                let (st, e2, mult) = self.rate_loop(ch, &an, frame, nb_subfr, cond, vad[ch][f], enc, target.max(16) * 8);
                *enc = e2;
                self.ch[ch].state = st;
                self.ch[ch].seed = self.ch[ch].seed.wrapping_add(1);
                if ch == 1 {
                    self.side_prev_uncoded = false;
                }
                // LBRR copy at a coarser step, from the state before it.
                if fec && vad[ch][f] {
                    let lcond = if f > 0 && new_lbrr.get(f - 1).is_some_and(|x: &[Option<LbrrFrame>; 2]| x[ch].is_some()) {
                        CondCoding::Conditional
                    } else {
                        CondCoding::Independent
                    };
                    let mut tmp = self.ch[ch].clone();
                    tmp.state = snapshot;
                    let lf = Self::capture(&tmp, &an, frame, nb_subfr, mult * 2.0, lcond);
                    while new_lbrr.len() <= f {
                        new_lbrr.push([None, None]);
                        new_lbrr_w.push([0, 0]);
                        new_lbrr_mid_only.push(false);
                    }
                    new_lbrr[f][ch] = Some(lf);
                    if ch == 0 {
                        new_lbrr_w[f] = wq;
                        new_lbrr_mid_only[f] = mid_only;
                    }
                }
            }
            if c == 2 && mid_only {
                self.side_prev_uncoded = true;
            }
        }
        if fec && !new_lbrr.is_empty() {
            while new_lbrr.len() < nf {
                new_lbrr.push([None, None]);
                new_lbrr_w.push([0, 0]);
                new_lbrr_mid_only.push(false);
            }
            self.lbrr = new_lbrr;
            self.lbrr_w = new_lbrr_w;
            self.lbrr_mid_only = new_lbrr_mid_only;
        } else {
            self.lbrr.clear();
        }
    }

    /// The indices and pulses of a frame coded at `mult` (for LBRR).
    fn capture(ce: &ChannelEncoder, an: &Analysis, frame: &[f32], nb_subfr: usize, mult: f64, cond: CondCoding) -> LbrrFrame {
        let prev_type = ce.state.prev_signal_type;
        let prev_lag = ce.state.prev_lag;
        let (ix, q) = ce.indices_and_pulses(an, frame, nb_subfr, mult, cond, true);
        LbrrFrame { ix, q, cond, prev_type, prev_lag }
    }

    fn encode_weight_indices(enc: &mut RangeEncoder, idx: [i32; 5]) {
        enc.icdf(idx[0] as usize, &STEREO_STAGE1_ICDF, 8);
        enc.icdf(idx[1] as usize, &STEREO_STAGE2_ICDF, 8);
        enc.icdf(idx[2] as usize, &STEREO_STAGE3_ICDF, 8);
        enc.icdf(idx[3] as usize, &STEREO_STAGE2_ICDF, 8);
        enc.icdf(idx[4] as usize, &STEREO_STAGE3_ICDF, 8);
    }

    /// Encodes quantized weights (Q13) by finding their indices again.
    fn encode_weights(enc: &mut RangeEncoder, w: [i32; 2]) {
        let (wi1, i3, q1) = quantize_weight(f64::from(w[1]) / 8192.0);
        let (wi0, i1, _) = quantize_weight(f64::from(w[0] + q1) / 8192.0);
        let n = (wi0 / 3) * 5 + wi1 / 3;
        Self::encode_weight_indices(enc, [n as i32, (wi0 % 3) as i32, i1, (wi1 % 3) as i32, i3]);
    }

    /// Searches the gain multiplier so the frame uses about `target` 1/8
    /// bits; returns the committed state, coder and multiplier.
    #[allow(clippy::too_many_arguments)]
    fn rate_loop(
        &self,
        ch: usize,
        an: &Analysis,
        frame: &[f32],
        nb_subfr: usize,
        cond: CondCoding,
        active: bool,
        enc: &RangeEncoder,
        target: i32,
    ) -> (ChannelState, RangeEncoder, f64) {
        let ce = &self.ch[ch];
        let mut lo = -4.0f64;
        let mut hi = 3.0f64;
        let mut best: Option<(ChannelState, RangeEncoder, f64)> = None;
        for _ in 0..7 {
            let mid = 0.5 * (lo + hi);
            let mult = mid.exp2();
            let (st, e2, bits, _) = ce.trial(an, frame, nb_subfr, mult, cond, active, enc, 0);
            if bits <= target {
                hi = mid;
                best = Some((st, e2, mult));
            } else {
                lo = mid;
            }
        }
        if let Some(b) = best {
            return b;
        }
        // Nothing in the normal range fits (an onset, a tiny budget):
        // coarsen until it does.
        let mut last = None;
        let top = if self.hard_limit { 12 } else { 6 };
        for k in 4..=top {
            let mult = f64::from(k).exp2();
            let (st, e2, bits, _) = ce.trial(an, frame, nb_subfr, mult, cond, active, enc, 0);
            if bits <= target {
                return (st, e2, mult);
            }
            last = Some((st, e2, mult));
        }
        last.expect("at least one trial")
    }
}

impl ChannelEncoder {
    /// The indices and pulses [`Self::trial`] would code (without coding).
    fn indices_and_pulses(
        &self,
        an: &Analysis,
        frame: &[f32],
        nb_subfr: usize,
        mult: f64,
        cond: CondCoding,
        active: bool,
    ) -> (FrameIndices, Vec<i32>) {
        let scratch = RangeEncoder::new(1275);
        let (_, e2, _, _) = self.trial(an, frame, nb_subfr, mult, cond, active, &scratch, 0);
        // Decode the scratch frame to recover exactly what was coded.
        let bytes = e2.finish();
        let mut dec = crate::range::RangeDecoder::new(&bytes);
        let fs = self.state.fs_khz;
        let ix = super::decoder::decode_indices(&mut dec, fs, nb_subfr, active, cond, self.state.prev_signal_type, self.state.prev_lag);
        let raw = super::decoder::decode_pulses(&mut dec, ix.signal_type, ix.qoff, frame.len());
        (ix, raw)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packet::Bandwidth;
    use crate::silk::SilkDecoder;

    fn speechish(n: usize, fs: f32) -> Vec<f32> {
        let mut out = vec![0.0f32; n];
        let (mut y1, mut y2) = (0.0f32, 0.0f32);
        let mut phase = 0.0f32;
        for (i, o) in out.iter_mut().enumerate() {
            let t = i as f32 / fs;
            phase += 130.0 / fs;
            let mut e = 0.0;
            if phase >= 1.0 {
                phase -= 1.0;
                e = 1.0;
            }
            let r = 0.97f32;
            let th = 2.0 * std::f32::consts::PI * 600.0 / fs;
            let y = e + 2.0 * r * th.cos() * y1 - r * r * y2;
            y2 = y1;
            y1 = y;
            *o = 0.05 * y * (1.0 + 0.3 * (t * 3.0).sin());
        }
        out
    }

    /// SILK frames decode to what the encoder's closed loop synthesized:
    /// encoder and decoder share one synthesis, so the decoded signal tracks
    /// the input at the coded SNR.
    #[test]
    fn silk_round_trip_internal_rate() {
        for fs in [8usize, 12, 16] {
            let bw = match fs {
                8 => Bandwidth::Narrow,
                12 => Bandwidth::Medium,
                _ => Bandwidth::Wide,
            };
            let n = 20 * fs;
            let x = speechish(n * 60 + 5 * fs, fs as f32 * 1000.0);
            let mut enc = SilkEncoder::new(1);
            let mut dec = SilkDecoder::new();
            let mut out = Vec::new();
            let mut bits = 0;
            for f in 0..60 {
                let input = &x[f * n..f * n + n + 5 * fs];
                let mut re = RangeEncoder::new(1275);
                enc.encode(&mut re, input, n, fs, 20, 600, false);
                bits += re.tell();
                let used = ((re.tell() + 7) >> 3) as usize;
                re.shrink(used);
                let bytes = re.finish();
                let mut rd = crate::range::RangeDecoder::new(&bytes);
                let y = dec.decode(Some(&mut rd), false, bw, 20, 1, false);
                out.extend_from_slice(&y[0]);
            }
            let (mut s, mut e) = (0.0f64, 0.0f64);
            for i in 10 * n..60 * n - 1 {
                let a = f64::from(x[i]);
                let b = f64::from(out[i + 1]);
                s += a * a;
                e += (a - b) * (a - b);
            }
            let snr = 10.0 * (s / e).log10();
            eprintln!("fs {fs}: SNR {snr:.2} dB at {} b/s", bits as f64 / 1.2);
        }
    }
}
