//! The CELT encoder (RFC 6716 §5.3): pre-emphasis, MDCT analysis with
//! transient detection, band energies, the coarse/fine energy quantizers,
//! the explicit allocation adjustments (band boost, trim, skipping), stereo
//! decisions, PVQ search, and the same bit allocation as the decoder.
//!
//! The pitch pre-filter is not used (the post-filter flag is always coded
//! as off), and the TF analysis is the simple one: long-block resolution
//! for steady frames, short-block resolution for transients.

use super::Synth;
use super::bands::{self, FrameParams, SPREAD_AGGRESSIVE, SPREAD_LIGHT, SPREAD_NONE, SPREAD_NORMAL};
use super::energy;
use super::mode::{BITRES, init_caps};
use super::rate::{self, EncoderChoices};
use super::tables::*;
use crate::range::RangeEncoder;

/// A CELT encoder for one Opus stream.
pub struct CeltEncoder {
    channels: usize,
    /// Pre-emphasized input of the last `OVERLAP` samples, per channel.
    in_mem: Vec<Vec<f32>>,
    preemph_mem: [f32; 2],
    old_band_e: [f32; 2 * NB_EBANDS],
    old_log_e: [f32; 2 * NB_EBANDS],
    old_log_e2: [f32; 2 * NB_EBANDS],
    prev_coded: usize,
    force_intra: bool,
    /// Energy of the last frame's sub-blocks, for transient detection.
    last_block_energy: f32,
    synth: Synth,
    prev_intensity: usize,
}

/// Per-frame choices the Opus layer makes for the CELT layer.
#[derive(Clone, Copy, Debug)]
pub struct FrameConfig {
    /// First coded band (17 in hybrid mode).
    pub start: usize,
    /// One past the last coded band (from the bandwidth).
    pub end: usize,
    /// Equivalent bit rate of the whole stream in bit/s (drives the
    /// intensity-stereo threshold).
    pub bitrate: i32,
}

impl CeltEncoder {
    pub fn new(channels: usize) -> Self {
        Self {
            channels,
            in_mem: vec![vec![0.0; OVERLAP]; channels],
            preemph_mem: [0.0; 2],
            old_band_e: [0.0; 2 * NB_EBANDS],
            old_log_e: [-28.0; 2 * NB_EBANDS],
            old_log_e2: [-28.0; 2 * NB_EBANDS],
            prev_coded: NB_EBANDS,
            force_intra: true,
            last_block_energy: 0.0,
            synth: Synth::new(),
            prev_intensity: NB_EBANDS,
        }
    }

    /// Back to the initial state (the decoder resets at the same points).
    pub fn reset(&mut self) {
        *self = Self::new(self.channels);
    }

    /// Encodes one frame of `n` samples per channel (interleaved, ±1.0) into
    /// `enc`, whose whole frame size bounds the bits used. Returns whether
    /// the frame was coded as silence.
    pub fn encode(&mut self, pcm: &[f32], n: usize, enc: &mut RangeEncoder, cfg: FrameConfig) -> bool {
        let c = self.channels;
        let lm = match n {
            120 => 0,
            240 => 1,
            480 => 2,
            _ => 3,
        };
        let mm = 1usize << lm;
        let (start, end) = (cfg.start, cfg.end);
        let total_bits = (enc.storage() * 8) as i32;
        // Pre-emphasis into an N + OVERLAP block per channel.
        let mut x: Vec<Vec<f32>> = Vec::with_capacity(c);
        let mut peak = 0.0f32;
        for ch in 0..c {
            let mut b = self.in_mem[ch].clone();
            let mut m = self.preemph_mem[ch];
            for i in 0..n {
                let s = pcm[i * c + ch] * 32768.0;
                peak = peak.max(s.abs());
                b.push(s - DEEMPHASIS * m);
                m = s;
            }
            self.preemph_mem[ch] = m;
            self.in_mem[ch].copy_from_slice(&b[n..n + OVERLAP]);
            x.push(b);
        }
        if c == 1 {
            for i in 0..NB_EBANDS {
                self.old_band_e[i] = self.old_band_e[i].max(self.old_band_e[NB_EBANDS + i]);
            }
        }
        let mut tell = enc.tell();
        let silence = peak < 1e-4;
        // Hybrid frames (tell > 1) cannot signal silence; such a frame is
        // coded as a quiet one.
        if tell == 1 {
            enc.bit_logp(silence, 15);
        }
        let silence = silence && tell == 1;
        if silence {
            self.old_band_e = [-28.0; 2 * NB_EBANDS];
            self.finish_frame(c, start, end, false);
            // The decoder needs the history of a silent frame too.
            self.last_block_energy = 0.0;
            return true;
        }
        tell = enc.tell();
        if start == 0 && tell + 16 <= total_bits {
            enc.bit_logp(false, 1);
            tell = enc.tell();
        }
        // Transient detection on the new samples.
        let mut transient = false;
        if lm > 0 && tell + 3 <= total_bits {
            transient = self.detect_transient(&x, n);
            enc.bit_logp(transient, 3);
            tell = enc.tell();
        } else {
            self.detect_transient(&x, n);
        }
        let blocks = if transient { mm } else { 1 };
        let coefs: Vec<Vec<f32>> = x.iter().map(|xc| self.synth.mdct_blocks(xc, n, blocks)).collect();
        // Band energies.
        let mut band_amp = [0.0f32; 2 * NB_EBANDS];
        let mut log_e = [0.0f32; 2 * NB_EBANDS];
        for ch in 0..c {
            for i in 0..end {
                let e: f32 = coefs[ch][mm * EBANDS[i]..mm * EBANDS[i + 1]].iter().map(|v| v * v).sum::<f32>();
                let a = (e + 1e-27).sqrt();
                band_amp[ch * NB_EBANDS + i] = a;
                log_e[ch * NB_EBANDS + i] = a.log2() - E_MEANS[i];
            }
        }
        // Intra (no time prediction) for the first frame or on a big
        // change.
        let mut dist = 0.0f32;
        for ch in 0..c {
            for i in start..end {
                let d = log_e[ch * NB_EBANDS + i] - self.old_band_e[ch * NB_EBANDS + i];
                dist += d * d;
            }
        }
        let want_intra = self.force_intra || dist / ((end - start) * c) as f32 > 9.0;
        let intra = if tell + 3 <= total_bits {
            enc.bit_logp(want_intra, 3);
            want_intra
        } else {
            false
        };
        let mut err = [0.0f32; 2 * NB_EBANDS];
        energy::quant_coarse(enc, &log_e, &mut self.old_band_e, &mut err, start, end, intra, c, lm, total_bits);
        // TF resolution.
        let mut tf_choice = [0i32; NB_EBANDS];
        if transient {
            tf_choice[start..end].fill(1);
        }
        let tf_res = tf_encode(enc, start, end, transient, &mut tf_choice, lm, total_bits);
        tell = enc.tell();
        // Spreading.
        let mut xn: Vec<Vec<f32>> = coefs.clone();
        for ch in 0..c {
            for i in 0..NB_EBANDS {
                let a = band_amp[ch * NB_EBANDS + i];
                for v in &mut xn[ch][mm * EBANDS[i]..mm * EBANDS[i + 1]] {
                    if i < end {
                        *v /= a;
                    } else {
                        *v = 0.0;
                    }
                }
            }
            for v in &mut xn[ch][mm * EBANDS[end]..] {
                *v = 0.0;
            }
        }
        let spread = spreading_decision(&xn, start, end, mm, transient);
        if tell + 4 <= total_bits {
            enc.icdf(spread as usize, &SPREAD_ICDF, 5);
        }
        // Band boosts (§5.3.4.1).
        let cap = init_caps(lm, c);
        let mut want_boost = [0i32; NB_EBANDS];
        let (t1, t2) = if lm >= 1 { (2.0, 4.0) } else { (3.0, 5.0) };
        for i in start + 1..end.saturating_sub(1) {
            let mut d = 0.0f32;
            for ch in 0..c {
                let e = |j: usize| log_e[ch * NB_EBANDS + j] + E_MEANS[j];
                d = d.max(2.0 * e(i) - e(i - 1) - e(i + 1));
            }
            want_boost[i] = if d > t2 { 2 } else if d > t1 { 1 } else { 0 };
        }
        let mut offsets = [0i32; NB_EBANDS];
        let mut dynalloc_logp = 6u32;
        let mut total_bits8 = total_bits << BITRES;
        let mut tell8 = enc.tell_frac();
        // Only boost when the frame has room for it.
        let rich = total_bits > (n as i32 / 4) * c as i32;
        for i in start..end {
            let width = (c * (EBANDS[i + 1] - EBANDS[i]) << lm) as i32;
            let quanta = (width << BITRES).min((6 << BITRES).max(width));
            let mut loop_logp = dynalloc_logp;
            let mut boost = 0;
            let mut j = 0;
            while tell8 + ((loop_logp as i32) << BITRES) < total_bits8 && boost < cap[i] {
                let flag = rich && j < want_boost[i];
                enc.bit_logp(flag, loop_logp);
                tell8 = enc.tell_frac();
                if !flag {
                    break;
                }
                boost += quanta;
                total_bits8 -= quanta;
                loop_logp = 1;
                j += 1;
            }
            offsets[i] = boost;
            if boost > 0 {
                dynalloc_logp = 2.max(dynalloc_logp - 1);
            }
        }
        // Allocation trim (§5.3.4.2).
        let trim = if tell8 + (6 << BITRES) <= total_bits8 {
            let t = alloc_trim(&log_e, &xn, start, end, c, mm);
            enc.icdf(t as usize, &TRIM_ICDF, 7);
            t
        } else {
            5
        };
        let mut bits = (total_bits << BITRES) - enc.tell_frac() - 1;
        let anti_collapse_rsv = if transient && lm >= 2 && bits >= (lm as i32 + 2) << BITRES { 1 << BITRES } else { 0 };
        bits -= anti_collapse_rsv;
        // Stereo decisions (§5.3.5).
        let (intensity, dual) = if c == 2 {
            let equiv_kbps = cfg.bitrate / 1000 - (80 * 48000 / n as i32) / 1000;
            let thresh = match equiv_kbps {
                ..35 => 8,
                35..50 => 12,
                50..68 => 16,
                68..84 => 18,
                84..102 => 19,
                102..130 => 20,
                _ => NB_EBANDS,
            };
            let thresh = thresh.clamp(start, end);
            self.prev_intensity = thresh;
            (thresh, use_dual_stereo(&xn, lm, mm))
        } else {
            (0, false)
        };
        let choices = EncoderChoices { intensity, dual_stereo: dual, prev_coded: self.prev_coded, signal_bandwidth: end - 1 };
        let alloc = rate::compute_allocation(start, end, &offsets, &cap, trim, bits, c, lm, enc, choices);
        energy::code_fine(enc, &mut self.old_band_e, &mut err, &alloc.fine_quant, start, end, c);
        let mut xs = xn[0].clone();
        let mut ys = if c == 2 { xn[1].clone() } else { Vec::new() };
        let params = FrameParams {
            start,
            end,
            lm,
            short_blocks: blocks,
            spread,
            dual_stereo: alloc.dual_stereo,
            intensity: alloc.intensity,
            tf_res: &tf_res,
            total_bits: (total_bits << BITRES) - anti_collapse_rsv,
            balance: alloc.balance,
            pulses: &alloc.pulses,
            coded_bands: alloc.coded_bands,
            disable_inv: false,
            resynth: false,
        };
        bands::quant_all_bands(enc, &params, &mut xs, &mut ys, &band_amp, 0);
        if anti_collapse_rsv > 0 {
            enc.bits(u32::from(transient), 1);
        }
        let left = total_bits - enc.tell();
        energy::code_finalise(enc, &mut self.old_band_e, &mut err, &alloc.fine_quant, &alloc.fine_priority, left, start, end, c);
        self.prev_coded = alloc.coded_bands;
        self.finish_frame(c, start, end, transient);
        false
    }

    /// The decoder's end-of-frame energy bookkeeping, mirrored.
    fn finish_frame(&mut self, c: usize, start: usize, end: usize, transient: bool) {
        if c == 1 {
            let (a, b) = self.old_band_e.split_at_mut(NB_EBANDS);
            b.copy_from_slice(a);
        }
        if !transient {
            self.old_log_e2 = self.old_log_e;
            self.old_log_e = self.old_band_e;
        } else {
            for i in 0..2 * NB_EBANDS {
                self.old_log_e[i] = self.old_log_e[i].min(self.old_band_e[i]);
            }
        }
        for ch in 0..2 {
            for i in (0..start).chain(end..NB_EBANDS) {
                self.old_band_e[ch * NB_EBANDS + i] = 0.0;
                self.old_log_e[ch * NB_EBANDS + i] = -28.0;
                self.old_log_e2[ch * NB_EBANDS + i] = -28.0;
            }
        }
        self.force_intra = false;
    }

    /// A transient is a sudden rise in the energy of the high-passed
    /// signal from one eighth of the frame to the next.
    fn detect_transient(&mut self, x: &[Vec<f32>], n: usize) -> bool {
        let sub = 8;
        let len = n / sub;
        let mut energies = vec![0.0f32; sub];
        for xc in x {
            let new = &xc[OVERLAP..];
            for (k, e) in energies.iter_mut().enumerate() {
                for i in k * len..(k + 1) * len {
                    let prev = if i == 0 { xc[OVERLAP - 1] } else { new[i - 1] };
                    let h = new[i] - prev;
                    *e += h * h;
                }
            }
        }
        let mut prev = self.last_block_energy;
        let mut transient = false;
        let floor = 1e3 * len as f32;
        for &e in &energies {
            if e > 8.0 * (prev + floor) {
                transient = true;
            }
            prev = 0.7 * prev + 0.3 * e;
            prev = prev.max(e * 0.5);
        }
        self.last_block_energy = *energies.last().unwrap_or(&0.0);
        transient
    }
}

/// §4.3.1 encoder side; returns the TF changes the decoder will derive.
fn tf_encode(
    enc: &mut RangeEncoder,
    start: usize,
    end: usize,
    transient: bool,
    choice: &mut [i32; NB_EBANDS],
    lm: usize,
    total_bits: i32,
) -> [i32; NB_EBANDS] {
    let mut budget = total_bits;
    let mut tell = enc.tell();
    let mut logp = if transient { 2 } else { 4 };
    let tf_select_rsv = lm > 0 && tell + logp + 1 <= budget;
    budget -= i32::from(tf_select_rsv);
    let mut curr = 0;
    let mut tf_changed = 0;
    for r in choice.iter_mut().take(end).skip(start) {
        if tell + logp <= budget {
            enc.bit_logp((*r ^ curr) != 0, logp as u32);
            tell = enc.tell();
            curr = *r;
            tf_changed |= curr;
        } else {
            *r = curr;
        }
        logp = if transient { 4 } else { 5 };
    }
    let ti = 4 * usize::from(transient);
    if tf_select_rsv && TF_SELECT[lm][ti + tf_changed as usize] != TF_SELECT[lm][ti + 2 + tf_changed as usize] {
        enc.bit_logp(false, 1);
    }
    let mut res = [0i32; NB_EBANDS];
    for i in start..end {
        res[i] = i32::from(TF_SELECT[lm][ti + choice[i] as usize]);
    }
    res
}

/// §5.3.7: how tonal the normalized spectrum is decides the spreading.
fn spreading_decision(x: &[Vec<f32>], start: usize, end: usize, mm: usize, transient: bool) -> u32 {
    if transient {
        return SPREAD_NORMAL;
    }
    let mut sum = 0.0f32;
    let mut count = 0;
    for xc in x {
        for i in start..end {
            let n = mm * (EBANDS[i + 1] - EBANDS[i]);
            if n <= 8 {
                continue;
            }
            let band = &xc[mm * EBANDS[i]..mm * EBANDS[i + 1]];
            let nf = n as f32;
            let mut t = [0usize; 3];
            for &v in band {
                let x2n = v * v * nf;
                if x2n < 0.25 {
                    t[0] += 1;
                }
                if x2n < 0.0625 {
                    t[1] += 1;
                }
                if x2n < 0.015625 {
                    t[2] += 1;
                }
            }
            // Fraction of small values: high for peaky (tonal) bands.
            sum += (t[0] + t[1] + t[2]) as f32 / (3.0 * nf);
            count += 1;
        }
    }
    if count == 0 {
        return SPREAD_NORMAL;
    }
    let avg = sum / count as f32;
    if avg > 0.7 {
        SPREAD_NONE
    } else if avg > 0.55 {
        SPREAD_LIGHT
    } else if avg > 0.3 {
        SPREAD_NORMAL
    } else {
        SPREAD_AGGRESSIVE
    }
}

/// §5.3.4.2: the allocation trim from the spectral tilt and, for stereo,
/// the inter-channel correlation at low frequencies.
fn alloc_trim(log_e: &[f32], x: &[Vec<f32>], start: usize, end: usize, c: usize, mm: usize) -> i32 {
    let mut trim = 5.0f32;
    if c == 2 {
        let mut corr = 0.0f32;
        for i in 0..8.min(end) {
            let (lo, hi) = (mm * EBANDS[i], mm * EBANDS[i + 1]);
            let d: f32 = x[0][lo..hi].iter().zip(&x[1][lo..hi]).map(|(a, b)| a * b).sum();
            corr += d;
        }
        let corr = (corr / 8.0).abs().min(1.0);
        let logxc = (1.001 - corr * corr).log2();
        trim += (0.75 * logxc).max(-4.0);
    }
    // Spectral tilt: positive when the high bands are relatively louder.
    let mut diff = 0.0f32;
    let span = (end - start).max(2) as f32;
    for ch in 0..c {
        for i in start..end.saturating_sub(1) {
            diff += (log_e[ch * NB_EBANDS + i] + E_MEANS[i]) * (2.0 + 2.0 * (i - start) as f32 - span);
        }
    }
    diff /= c as f32 * (span - 1.0).max(1.0) * span;
    trim -= (diff / 6.0).clamp(-2.0, 2.0);
    (trim + 0.5).floor().clamp(0.0, 10.0) as i32
}

/// §5.3.5: dual (L/R) stereo when the L1 norms favour it over the first 13
/// bands.
fn use_dual_stereo(x: &[Vec<f32>], lm: usize, mm: usize) -> bool {
    let mut l1_lr = 0.0f32;
    let mut l1_ms = 0.0f32;
    let mut bins = 0usize;
    for i in 0..13 {
        for j in mm * EBANDS[i]..mm * EBANDS[i + 1] {
            let (l, r) = (x[0][j], x[1][j]);
            l1_lr += l.abs() + r.abs();
            l1_ms += (l + r).abs() * std::f32::consts::FRAC_1_SQRT_2 + (l - r).abs() * std::f32::consts::FRAC_1_SQRT_2;
            bins += 1;
        }
    }
    let e = if lm > 1 { 13.0 } else { 5.0 };
    // Mid/side iff L1_ms / (bins + E) < L1_lr / bins.
    !(l1_ms / (bins as f32 + e) < l1_lr / bins as f32)
}
