//! Band energy coding (RFC 6716 §4.3.2): the Laplace-coded coarse energy
//! with time and frequency prediction, fine energy, and the final bits.
//!
//! Energies are log2 amplitudes relative to [`E_MEANS`](super::tables::E_MEANS)
//! (1.0 = 6.02 dB), stored per band and channel as `e[c * 21 + band]`.

use super::rate::MAX_FINE_BITS;
use super::tables::{BETA_COEF, BETA_INTRA, E_PROB_MODEL, NB_EBANDS, PRED_COEF, SMALL_ENERGY_ICDF};
use crate::range::{Coder, RangeDecoder, RangeEncoder};

const LAPLACE_MINP: u32 = 1;
const LAPLACE_NMIN: u32 = 16;

fn laplace_freq1(fs0: u32, decay: u32) -> u32 {
    let ft = 32768 - LAPLACE_MINP * (2 * LAPLACE_NMIN) - fs0;
    (ft * (16384 - decay)) >> 15
}

/// Decodes a Laplace-distributed integer with probability of zero
/// `fs/32768` and decay `decay/16384`.
pub fn laplace_decode(dec: &mut RangeDecoder, mut fs: u32, decay: u32) -> i32 {
    let mut val = 0i32;
    let fm = dec.decode_bin(15);
    let mut fl = 0;
    if fm >= fs {
        val += 1;
        fl = fs;
        fs = laplace_freq1(fs, decay) + LAPLACE_MINP;
        while fs > LAPLACE_MINP && fm >= fl + 2 * fs {
            fs *= 2;
            fl += fs;
            fs = (((fs - 2 * LAPLACE_MINP) * decay) >> 15) + LAPLACE_MINP;
            val += 1;
        }
        if fs <= LAPLACE_MINP {
            let di = (fm - fl) >> 1;
            val += di as i32;
            fl += 2 * di * LAPLACE_MINP;
        }
        if fm < fl + fs {
            val = -val;
        } else {
            fl += fs;
        }
    }
    dec.update(fl, (fl + fs).min(32768), 32768);
    val
}

/// Encodes `val` with the model of [`laplace_decode`], clamping it to the
/// largest magnitude the model can represent; returns the value coded.
pub fn laplace_encode(enc: &mut RangeEncoder, val: i32, fs0: u32, decay: u32) -> i32 {
    // Walk the decoder's intervals: magnitude m has [fl, fl+fs) for -m and
    // [fl+fs, fl+2fs) for +m.
    if val == 0 {
        enc.encode_bin(0, fs0, 15);
        return 0;
    }
    let mag = val.unsigned_abs();
    let mut fl = fs0;
    let mut fs = laplace_freq1(fs0, decay) + LAPLACE_MINP;
    let mut m = 1u32;
    while m < mag && fs > LAPLACE_MINP {
        fl += 2 * fs;
        fs = (((2 * fs - 2 * LAPLACE_MINP) * decay) >> 15) + LAPLACE_MINP;
        m += 1;
    }
    if m < mag {
        // fs == 1: each further magnitude takes 2 slots.
        let max_extra = (32768 - fl) / 2 - 1;
        let extra = (mag - m).min(max_extra);
        fl += 2 * extra;
        m += extra;
    }
    if fl + 2 * fs > 32768 {
        // Cannot happen with the codec's models, but stay inside the range.
        fs = (32768 - fl) / 2;
    }
    let neg = val < 0;
    let lo = if neg { fl } else { fl + fs };
    enc.encode_bin(lo, lo + fs, 15);
    if neg { -(m as i32) } else { m as i32 }
}

/// §4.3.2.1: decodes the coarse energies of bands `start..end` into
/// `old_e` (which holds the previous frame's final energies).
pub fn unquant_coarse(
    dec: &mut RangeDecoder,
    old_e: &mut [f32],
    start: usize,
    end: usize,
    intra: bool,
    c: usize,
    lm: usize,
) {
    let model = &E_PROB_MODEL[lm][usize::from(intra)];
    let (coef, beta) = if intra { (0.0, BETA_INTRA) } else { (PRED_COEF[lm], BETA_COEF[lm]) };
    let budget = (dec.storage() * 8) as i32;
    let mut prev = [0.0f32; 2];
    for i in start..end {
        for ch in 0..c {
            let tell = dec.tell();
            let qi = if budget - tell >= 15 {
                let pi = 2 * i.min(20);
                laplace_decode(dec, u32::from(model[pi]) << 7, u32::from(model[pi + 1]) << 6)
            } else if budget - tell >= 2 {
                let q = dec.icdf(&SMALL_ENERGY_ICDF, 2) as i32;
                (q >> 1) ^ -(q & 1)
            } else if budget - tell >= 1 {
                -i32::from(dec.bit_logp(1))
            } else {
                -1
            };
            let q = qi as f32;
            let e = &mut old_e[i + ch * NB_EBANDS];
            *e = e.max(-9.0);
            let tmp = coef * *e + prev[ch] + q;
            *e = tmp;
            prev[ch] = prev[ch] + q - beta * q;
        }
    }
}

/// §4.3.2.1 encoder side: quantizes `target` energies of bands
/// `start..end` against the prediction from `old_e`, writing the coded
/// energies into `old_e` and the quantization error into `error`.
#[allow(clippy::too_many_arguments)]
pub fn quant_coarse(
    enc: &mut RangeEncoder,
    target: &[f32],
    old_e: &mut [f32],
    error: &mut [f32],
    start: usize,
    end: usize,
    intra: bool,
    c: usize,
    lm: usize,
    budget: i32,
) {
    let model = &E_PROB_MODEL[lm][usize::from(intra)];
    let (coef, beta) = if intra { (0.0, BETA_INTRA) } else { (PRED_COEF[lm], BETA_COEF[lm]) };
    let mut prev = [0.0f32; 2];
    for i in start..end {
        for ch in 0..c {
            let idx = i + ch * NB_EBANDS;
            let x = target[idx];
            let old = old_e[idx].max(-9.0);
            let f = x - coef * old - prev[ch];
            let mut qi = (f + 0.5).floor() as i32;
            // Do not let the energy drop much below the previous frame's
            // (keeps the decoder's energies from collapsing).
            let decay_bound = old_e[idx].max(-28.0) - 28.0;
            if qi < 0 && x < decay_bound {
                qi += (decay_bound - x) as i32;
                if qi > 0 {
                    qi = 0;
                }
            }
            let tell = enc.tell();
            let left = budget - tell;
            if left >= 15 {
                let pi = 2 * i.min(20);
                qi = laplace_encode(enc, qi, u32::from(model[pi]) << 7, u32::from(model[pi + 1]) << 6);
            } else if left >= 2 {
                qi = qi.clamp(-1, 1);
                let s = (2 * qi) ^ -i32::from(qi < 0);
                enc.icdf(s as usize, &SMALL_ENERGY_ICDF, 2);
            } else if left >= 1 {
                qi = qi.clamp(-1, 0);
                enc.bit_logp(qi != 0, 1);
            } else {
                qi = -1;
            }
            let q = qi as f32;
            error[idx] = f - q;
            old_e[idx] = coef * old + prev[ch] + q;
            prev[ch] = prev[ch] + q - beta * q;
        }
    }
}

/// §4.3.2.2: the fine energy refinement.
pub fn code_fine<E: Coder>(
    ec: &mut E,
    old_e: &mut [f32],
    error: &mut [f32],
    fine_quant: &[i32; NB_EBANDS],
    start: usize,
    end: usize,
    c: usize,
) {
    for i in start..end {
        let fq = fine_quant[i];
        if fq <= 0 {
            continue;
        }
        let levels = 1i32 << fq;
        for ch in 0..c {
            let idx = i + ch * NB_EBANDS;
            let want = if E::ENCODE {
                (((error[idx] + 0.5) * levels as f32).floor() as i32).clamp(0, levels - 1) as u32
            } else {
                0
            };
            let q2 = ec.bits(want, fq as u32);
            let offset = (q2 as f32 + 0.5) / levels as f32 - 0.5;
            old_e[idx] += offset;
            error[idx] -= offset;
        }
    }
}

/// §4.3.2.2: the final fine bits, given `bits_left` whole bits.
#[allow(clippy::too_many_arguments)]
pub fn code_finalise<E: Coder>(
    ec: &mut E,
    old_e: &mut [f32],
    error: &mut [f32],
    fine_quant: &[i32; NB_EBANDS],
    fine_priority: &[i32; NB_EBANDS],
    mut bits_left: i32,
    start: usize,
    end: usize,
    c: usize,
) {
    for prio in 0..2 {
        let mut i = start;
        while i < end && bits_left >= c as i32 {
            if fine_quant[i] >= MAX_FINE_BITS || fine_priority[i] != prio {
                i += 1;
                continue;
            }
            for ch in 0..c {
                let idx = i + ch * NB_EBANDS;
                let want = u32::from(error[idx] >= 0.0);
                let q2 = ec.bits(want, 1);
                let offset = (q2 as f32 - 0.5) / (1u32 << (fine_quant[i] + 1)) as f32;
                old_e[idx] += offset;
                error[idx] -= offset;
                bits_left -= 1;
            }
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn laplace_round_trip() {
        for &(fs, decay) in &[(72u32 << 7, 127u32 << 6), (24 << 7, 179 << 6), (177 << 7, 11 << 6), (42 << 7, 121 << 6)] {
            let vals: Vec<i32> = (-40..=40).chain([-200, 300, 1000, -1000]).collect();
            let mut enc = RangeEncoder::new(2000);
            let mut coded = Vec::new();
            for &v in &vals {
                coded.push(laplace_encode(&mut enc, v, fs, decay));
            }
            let bytes = enc.finish();
            let mut dec = RangeDecoder::new(&bytes);
            for (i, &c) in coded.iter().enumerate() {
                assert_eq!(laplace_decode(&mut dec, fs, decay), c, "fs {fs} value {}", vals[i]);
                // Small values are always representable; large ones saturate
                // where the model's tail runs out.
                if vals[i].abs() <= 10 {
                    assert_eq!(c, vals[i]);
                }
                assert_eq!(c.signum(), vals[i].signum());
            }
        }
    }
}
