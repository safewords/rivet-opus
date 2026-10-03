//! Band energies: the Laplace-coded coarse energy with its time/frequency
//! predictor, the fine energy bits and the final fine bits
//! (CELT_SPEC §2, §3; RFC 6716 §4.3.2, §5.3.2).
//!
//! Each procedure serves both sides through [`Coder`]: the decoder reads
//! the values, the encoder chooses them from its target energies and
//! writes them. Energies are log2 amplitudes relative to the band means
//! `E_MEANS`, kept per channel as `[channel][band]`.

use super::tables::{BETA_COEF, BETA_INTRA, E_PROB_MODEL, MAX_FINE_BITS, NB_EBANDS, PRED_COEF, SMALL_ENERGY_ICDF};
use crate::range::Coder;

/// Per-channel band energies.
pub(crate) type BandEnergies = [[f32; NB_EBANDS]; 2];

/// Smallest probability of a Laplace value (CELT_SPEC §2.1).
const LAPLACE_MINP: u32 = 1;
/// Values guaranteed the minimum probability in the tail.
const LAPLACE_NMIN: u32 = 16;

/// Probability of the first nonzero magnitude, from the probability of 0
/// (CELT_SPEC §2.1 step 2.1).
fn laplace_first(fs0: u32, decay: u32) -> u32 {
    (((32768 - 2 * LAPLACE_NMIN * LAPLACE_MINP - fs0) * (16384 - decay)) >> 15) + LAPLACE_MINP
}

/// Probability of the next magnitude (CELT_SPEC §2.1 step 2.2).
fn laplace_next(fs: u32, decay: u32) -> u32 {
    (((2 * fs - 2 * LAPLACE_MINP) * decay) >> 15) + LAPLACE_MINP
}

/// The Laplace-like symbol of CELT_SPEC §2.1, `fs0` being the probability
/// of 0 in 1/32768 and `decay` the Q14 ratio of successive magnitudes.
/// The encoder codes `val` clamped to what the distribution can represent
/// and returns the coded value; the decoder returns the decoded one.
fn laplace<C: Coder>(ec: &mut C, val: i32, fs0: u32, decay: u32) -> i32 {
    if C::ENCODE {
        let (fl, fs, coded) = laplace_interval(val, fs0, decay);
        ec.code(fl, (fl + fs).min(32768), 32768);
        coded
    } else {
        let fm = ec.decode_fs(32768);
        let mut fl = 0;
        let mut fs = fs0;
        let mut v = 0i32;
        if fm >= fs {
            v = 1;
            fl = fs;
            fs = laplace_first(fs, decay);
            while fs > LAPLACE_MINP && fm >= fl + 2 * fs {
                fl += 2 * fs;
                fs = laplace_next(fs, decay);
                v += 1;
            }
            if fs <= LAPLACE_MINP {
                let di = (fm - fl) >> 1;
                v += di as i32;
                fl += 2 * di;
            }
            if fm < fl + fs {
                v = -v;
            } else {
                fl += fs;
            }
        }
        ec.code(fl, (fl + fs).min(32768), 32768);
        v
    }
}

/// The interval `(fl, fs)` of `val` under the distribution of
/// [`laplace`], with the value clamped into the representable range.
fn laplace_interval(val: i32, fs0: u32, decay: u32) -> (u32, u32, i32) {
    if val == 0 {
        return (0, fs0, 0);
    }
    let mag = val.unsigned_abs();
    let mut fl = fs0;
    let mut fs = laplace_first(fs0, decay);
    let mut k = 1u32;
    while k < mag && fs > LAPLACE_MINP {
        fl += 2 * fs;
        fs = laplace_next(fs, decay);
        k += 1;
    }
    if k < mag {
        // The tail: every further magnitude has the minimum probability,
        // up to the end of the range (the positive value must still fit).
        let di = (mag - k).min((32766 - fl) >> 1);
        k += di;
        fl += 2 * di;
    }
    if val < 0 { (fl, fs, -(k as i32)) } else { (fl + fs, fs, k as i32) }
}

/// Coarse energy (CELT_SPEC §2.3; RFC 6716 §4.3.2.1, §5.3.2): for bands
/// `start..end` and the `c` coded channels, the integer residual of the
/// time/frequency prediction. `old` holds the previous frame's energies on
/// entry and the coarse energies on return. The encoder quantizes
/// `target` and leaves the remaining error in `err`; the decoder ignores
/// both.
pub(crate) fn code_coarse<C: Coder>(
    ec: &mut C,
    budget: i32,
    old: &mut BandEnergies,
    target: &BandEnergies,
    err: &mut BandEnergies,
    start: usize,
    end: usize,
    intra: bool,
    c: usize,
    lm: usize,
) {
    let (alpha, beta) = if intra { (0.0, BETA_INTRA) } else { (PRED_COEF[lm], BETA_COEF[lm]) };
    let model = &E_PROB_MODEL[lm][usize::from(intra)];
    let mut prev = [0.0f32; 2];
    for i in start..end {
        for ch in 0..c {
            let t = tell(ec);
            let old_e = old[ch][i].max(-9.0);
            let pred = alpha * old_e + prev[ch];
            let mut qi = 0;
            if C::ENCODE {
                let f = target[ch][i] - pred;
                qi = (f + 0.5).floor() as i32;
                // Encoder policy (CELT_SPEC §12): keep the last bands
                // cheap when the frame is nearly spent.
                let bits_left = budget - t - 3 * (c * (end - i)) as i32;
                if i != start && bits_left < 24 {
                    qi = qi.min(1);
                    if bits_left < 16 {
                        qi = qi.max(-1);
                    }
                }
            }
            let left = budget - t;
            let qi = if left >= 15 {
                let k = i.min(20);
                laplace(ec, qi, u32::from(model[2 * k]) << 7, u32::from(model[2 * k + 1]) << 6)
            } else if left >= 2 {
                // PDF {2, 1, 1}/4 over the values 0, −1, +1.
                let qi = qi.clamp(-1, 1);
                let s = if C::ENCODE {
                    match qi {
                        0 => 0,
                        -1 => 1,
                        _ => 2,
                    }
                } else {
                    let fs = ec.decode_fs(4);
                    (0..2).find(|&k| fs < 4 - u32::from(SMALL_ENERGY_ICDF[k])).unwrap_or(2)
                };
                let fl = if s == 0 { 0 } else { 4 - u32::from(SMALL_ENERGY_ICDF[s - 1]) };
                ec.code(fl, 4 - u32::from(SMALL_ENERGY_ICDF[s]), 4);
                (s >> 1) as i32 ^ -((s & 1) as i32)
            } else if left >= 1 {
                -i32::from(ec.bit_logp(qi < 0, 1))
            } else {
                -1
            };
            let q = qi as f32;
            old[ch][i] = pred + q;
            prev[ch] = prev[ch] + q - beta * q;
            if C::ENCODE {
                err[ch][i] = target[ch][i] - old[ch][i];
            }
        }
    }
}

/// Whole bits used so far, from the eighth-bit count: `ec_tell` equals
/// `ec_tell_frac` rounded up to whole bits (RFC 6716 §4.1.6.2 computes the
/// latter as 8·`nbits_total` minus a 3-bit refinement of `8·ilog(rng)`).
pub(crate) fn tell<C: Coder>(ec: &C) -> i32 {
    (ec.tell_frac() + 7) >> 3
}

/// Fine energy (CELT_SPEC §3.1; RFC 6716 §4.3.2.2): `fine_quant[i]` raw
/// bits per channel refine each band. The encoder codes the remaining
/// error `err` and reduces it.
pub(crate) fn code_fine<C: Coder>(
    ec: &mut C,
    old: &mut BandEnergies,
    err: &mut BandEnergies,
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
        let steps = 1i32 << fq;
        for ch in 0..c {
            let value = if C::ENCODE { (((err[ch][i] + 0.5) * steps as f32).floor() as i32).clamp(0, steps - 1) } else { 0 };
            let q2 = ec.bits(value as u32, fq as u32);
            let offset = (q2 as f32 + 0.5) * (1 << (14 - fq)) as f32 * (1.0 / 16384.0) - 0.5;
            old[ch][i] += offset;
            err[ch][i] -= offset;
        }
    }
}

/// The final fine bits (CELT_SPEC §3.2; RFC 6716 §4.3.2.2): the bits left
/// at the end of the frame give one more bit per channel to the bands with
/// fewer than eight fine bits, priority 0 first.
pub(crate) fn code_final<C: Coder>(
    ec: &mut C,
    old: &mut BandEnergies,
    err: &mut BandEnergies,
    fine_quant: &[i32; NB_EBANDS],
    fine_priority: &[bool; NB_EBANDS],
    mut bits_left: i32,
    start: usize,
    end: usize,
    c: usize,
) {
    for prio in [false, true] {
        for i in start..end {
            if bits_left < c as i32 {
                break;
            }
            if fine_quant[i] >= MAX_FINE_BITS || fine_priority[i] != prio {
                continue;
            }
            for ch in 0..c {
                let q2 = ec.bits(u32::from(C::ENCODE && err[ch][i] >= 0.0), 1);
                let offset = (q2 as f32 - 0.5) * (1 << (14 - fine_quant[i] - 1)) as f32 * (1.0 / 16384.0);
                old[ch][i] += offset;
                err[ch][i] -= offset;
                bits_left -= 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::range::{RangeDecoder, RangeEncoder};

    /// Every value, including clamped tail values, decodes as the encoder
    /// reports it.
    #[test]
    fn laplace_round_trips() {
        for &(fs0, decay) in &[(72u32 << 7, 127u32 << 6), (24 << 7, 179 << 6), (177 << 7, 11 << 6), (21 << 7, 178 << 6)] {
            let vals: Vec<i32> = (-40..=40).chain([-1000, 1000, 300, -300]).collect();
            let mut enc = RangeEncoder::new(4000);
            let coded: Vec<i32> = vals.iter().map(|&v| laplace(&mut enc, v, fs0, decay)).collect();
            let buf = enc.finish();
            let mut dec = RangeDecoder::new(&buf);
            for (&v, &c) in vals.iter().zip(&coded) {
                assert!(c.signum() == v.signum() && c.abs() <= v.abs());
                assert_eq!(laplace(&mut dec, 0, fs0, decay), c, "fs {fs0} decay {decay} value {v}");
            }
        }
    }
}
