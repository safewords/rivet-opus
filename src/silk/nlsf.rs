//! Normalized LSF reconstruction and conversion to LPC coefficients (RFC
//! 6716 §4.2.7.5), in the RFC's fixed-point arithmetic, which decoders
//! SHOULD match exactly (the encoder uses the same code so that its
//! prediction tracks the decoder's).

use super::tables::*;
use crate::range::ilog;

/// Whether a bandwidth uses the order-16 (WB) tables.
#[inline]
pub fn order(wb: bool) -> usize {
    if wb { 16 } else { 10 }
}

/// The stage-1 codebook vector `i1` (Q8).
pub fn cb1(wb: bool, i1: usize) -> &'static [i32] {
    if wb {
        &NLSF_CB1_WB_Q8[i1]
    } else {
        &NLSF_CB1_NB_Q8[i1]
    }
}

/// The stage-2 inverse CDF for coefficient `k` given `i1`.
pub fn stage2_icdf(wb: bool, i1: usize, k: usize) -> &'static [u8] {
    if wb {
        &NLSF_STAGE2_WB_ICDF[usize::from(NLSF_STAGE2_WB_SELECT[i1][k])]
    } else {
        &NLSF_STAGE2_NB_ICDF[usize::from(NLSF_STAGE2_NB_SELECT[i1][k])]
    }
}

/// The backward prediction weight (Q8) of coefficient `k < order - 1`.
pub fn pred_q8(wb: bool, i1: usize, k: usize) -> i32 {
    if wb {
        NLSF_PRED_WB_Q8[usize::from(NLSF_PRED_WB_SELECT[i1][k])][k]
    } else {
        NLSF_PRED_NB_Q8[usize::from(NLSF_PRED_NB_SELECT[i1][k])][k]
    }
}

/// The stage-2 step size (Q16).
pub fn qstep(wb: bool) -> i32 {
    if wb { 9830 } else { 11796 }
}

/// §4.2.7.5.2: undoes the backward prediction of the stage-2 indices.
pub fn residual_q10(wb: bool, i1: usize, i2: &[i32]) -> [i32; 16] {
    let d = order(wb);
    let q = qstep(wb);
    let mut res = [0i32; 16];
    for k in (0..d).rev() {
        let pred = if k + 1 < d {
            (res[k + 1] * pred_q8(wb, i1, k)) >> 8
        } else {
            0
        };
        let i = i2[k];
        res[k] = pred + ((((i << 10) - i.signum() * 102) * q) >> 16);
    }
    res
}

/// §4.2.7.5.3: the IHMW weights (Q9) of a stage-1 codebook vector.
pub fn weights_q9(cb: &[i32], d: usize) -> [i32; 16] {
    let mut w = [0i32; 16];
    for k in 0..d {
        let prev = if k == 0 { 0 } else { cb[k - 1] };
        let next = if k + 1 == d { 256 } else { cb[k + 1] };
        let w2 = (1024 / (cb[k] - prev) + 1024 / (next - cb[k])) << 16;
        let i = ilog(w2 as u32);
        let f = (w2 >> (i - 8)) & 127;
        let y = (if i & 1 != 0 { 32768 } else { 46214 }) >> ((32 - i) >> 1);
        w[k] = y + ((213 * f * y) >> 16);
    }
    w
}

/// §4.2.7.5.3: the normalized LSFs (Q15) from both stages, before
/// stabilization.
pub fn reconstruct(wb: bool, i1: usize, i2: &[i32]) -> [i32; 16] {
    let d = order(wb);
    let cb = cb1(wb, i1);
    let res = residual_q10(wb, i1, i2);
    let w = weights_q9(cb, d);
    let mut nlsf = [0i32; 16];
    for k in 0..d {
        nlsf[k] = ((cb[k] << 7) + (res[k] << 14) / w[k]).clamp(0, 32767);
    }
    nlsf
}

/// §4.2.7.5.4 (with RFC 8251 §7): enforces the minimum spacing.
pub fn stabilize(nlsf: &mut [i32], wb: bool) {
    let d = order(wb);
    let dmin: &[i32] = if wb {
        &NLSF_MIN_SPACING_WB
    } else {
        &NLSF_MIN_SPACING_NB
    };
    for _ in 0..20 {
        // Find the worst violation (lowest index on ties).
        let mut min_diff = i32::MAX;
        let mut idx = 0;
        for i in 0..=d {
            let lo = if i == 0 { 0 } else { nlsf[i - 1] };
            let hi = if i == d { 32768 } else { nlsf[i] };
            let diff = hi - lo - dmin[i];
            if diff < min_diff {
                min_diff = diff;
                idx = i;
            }
        }
        if min_diff >= 0 {
            return;
        }
        if idx == 0 {
            nlsf[0] = dmin[0];
        } else if idx == d {
            nlsf[d - 1] = 32768 - dmin[d];
        } else {
            let min_center: i32 = (dmin[idx] >> 1) + dmin[..idx].iter().sum::<i32>();
            let max_center: i32 = 32768 - (dmin[idx] >> 1) - dmin[idx + 1..=d].iter().sum::<i32>();
            let center = ((nlsf[idx - 1] + nlsf[idx] + 1) >> 1).clamp(min_center, max_center);
            nlsf[idx - 1] = center - (dmin[idx] >> 1);
            nlsf[idx] = nlsf[idx - 1] + dmin[idx];
        }
    }
    // Fallback.
    nlsf[..d].sort_unstable();
    nlsf[0] = nlsf[0].max(dmin[0]);
    for k in 1..d {
        nlsf[k] = nlsf[k].max((nlsf[k - 1] + dmin[k]).min(32767));
    }
    nlsf[d - 1] = nlsf[d - 1].min(32768 - dmin[d]);
    for k in (0..d - 1).rev() {
        nlsf[k] = nlsf[k].min(nlsf[k + 1] - dmin[k + 1]);
    }
}

/// `a[k] = (a[k] * sc[k]) >> 16` with `sc[k+1] = (sc[0] sc[k] + 32768) >> 16`
/// (silk_bwexpander_32 of §4.2.7.5.7).
pub fn bwexpand32(a: &mut [i64], sc0: i64) {
    let mut chirp = sc0;
    let d = a.len();
    for v in a.iter_mut().take(d - 1) {
        *v = (*v * chirp) >> 16;
        chirp = (sc0 * chirp + 32768) >> 16;
    }
    a[d - 1] = (a[d - 1] * chirp) >> 16;
}

/// §4.2.7.5.8: whether the Q12 filter passes the stability test.
pub fn is_stable(a_q12: &[i32]) -> bool {
    let d = a_q12.len();
    let dc: i32 = a_q12.iter().sum();
    if dc > 4096 {
        return false;
    }
    let mut a: Vec<i64> = a_q12.iter().map(|&v| i64::from(v) << 12).collect();
    let mut inv_gain: i64 = 1 << 30;
    for k in (0..d).rev() {
        if a[k].abs() > 16_773_022 {
            return false;
        }
        let rc = -a[k] << 7;
        let div = (1i64 << 30) - ((rc * rc) >> 32);
        inv_gain = ((inv_gain * div) >> 32) << 2;
        if inv_gain < 107_374 {
            return false;
        }
        if k > 0 {
            let b1 = ilog(div as u32);
            let b2 = b1 - 16;
            let inv_qb2 = ((1i64 << 29) - 1) / (div >> (b2 + 1));
            let err_q29 = (1i64 << 29) - (((div << (15 - b2)) * inv_qb2) >> 16);
            let gain_qb1 = (inv_qb2 << 16) + ((err_q29 * inv_qb2) >> 13);
            let mut next = vec![0i64; k];
            for (n, nv) in next.iter_mut().enumerate() {
                let num = a[n] - ((a[k - n - 1] * rc + (1 << 30)) >> 31);
                // RFC 8251 §6: saturate the numerator, and call a row that
                // does not fit in 32 bits unstable.
                let num = num.clamp(i64::from(i32::MIN), i64::from(i32::MAX));
                let v = (num * gain_qb1 + (1i64 << (b1 - 1))) >> b1;
                if v > i64::from(i32::MAX) || v < i64::from(i32::MIN) {
                    return false;
                }
                *nv = v;
            }
            a.truncate(k);
            a.copy_from_slice(&next);
        }
    }
    true
}

/// §4.2.7.5.6–§4.2.7.5.8: normalized LSFs (Q15) to stable Q12 LPC
/// coefficients.
pub fn nlsf_to_lpc(nlsf: &[i32], wb: bool) -> Vec<i32> {
    let d = order(wb);
    let ordering: &[usize] = if wb { &LSF_ORDER_WB } else { &LSF_ORDER_NB };
    let mut c = [0i64; 16];
    for k in 0..d {
        let i = (nlsf[k] >> 8) as usize;
        let f = i64::from(nlsf[k] & 255);
        let c0 = i64::from(COS_Q12[i]);
        let c1 = i64::from(COS_Q12[i + 1]);
        c[ordering[k]] = (c0 * 256 + (c1 - c0) * f + 4) >> 3;
    }
    let d2 = d / 2;
    // p[j], q[j] for the current row, j = 0..=d2.
    let mut p = vec![0i64; d2 + 2];
    let mut q = vec![0i64; d2 + 2];
    p[0] = 1 << 16;
    q[0] = 1 << 16;
    p[1] = -c[0];
    q[1] = -c[1];
    for k in 1..d2 {
        let mut np = vec![0i64; d2 + 2];
        let mut nq = vec![0i64; d2 + 2];
        // Row k-1 holds j = 0..=k; extend it by symmetry: x[k+1] = x[k-1].
        let get = |row: &[i64], j: isize| -> i64 {
            if j < 0 {
                0
            } else if j as usize == k + 1 {
                row[k - 1]
            } else {
                row[j as usize]
            }
        };
        for j in 0..=k + 1 {
            let ji = j as isize;
            np[j] = get(&p, ji) + get(&p, ji - 2) - ((c[2 * k] * get(&p, ji - 1) + 32768) >> 16);
            nq[j] =
                get(&q, ji) + get(&q, ji - 2) - ((c[2 * k + 1] * get(&q, ji - 1) + 32768) >> 16);
        }
        p = np;
        q = nq;
    }
    let mut a32 = vec![0i64; d];
    for k in 0..d2 {
        let qd = q[k + 1] - q[k];
        let ps = p[k + 1] + p[k];
        a32[k] = -qd - ps;
        a32[d - k - 1] = qd - ps;
    }
    // §4.2.7.5.7: limit the range.
    let mut saturate = true;
    for _ in 0..10 {
        let (mut kmax, mut maxabs) = (0usize, 0i64);
        for (k, &v) in a32.iter().enumerate() {
            if v.abs() > maxabs {
                maxabs = v.abs();
                kmax = k;
            }
        }
        let maxabs_q12 = ((maxabs + 16) >> 5).min(163_838);
        if maxabs_q12 > 32767 {
            let sc = 65470 - ((maxabs_q12 - 32767) << 14) / ((maxabs_q12 * (kmax as i64 + 1)) >> 2);
            bwexpand32(&mut a32, sc);
        } else {
            saturate = false;
            break;
        }
    }
    if saturate {
        for v in a32.iter_mut() {
            *v = (((*v + 16) >> 5).clamp(-32768, 32767)) << 5;
        }
    }
    // §4.2.7.5.8: limit the prediction gain.
    for i in 0..16 {
        let a_q12: Vec<i32> = a32.iter().map(|&v| ((v + 16) >> 5) as i32).collect();
        if is_stable(&a_q12) {
            return a_q12;
        }
        bwexpand32(&mut a32, 65536 - (2i64 << i));
    }
    a32.iter().map(|&v| ((v + 16) >> 5) as i32).collect()
}

/// §4.2.7.5.5: LSF interpolation for the first half of a 20 ms frame.
pub fn interpolate(n0: &[i32], n2: &[i32], w_q2: i32, d: usize) -> [i32; 16] {
    let mut n1 = [0i32; 16];
    for k in 0..d {
        n1[k] = n0[k] + ((w_q2 * (n2[k] - n0[k])) >> 2);
    }
    n1
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The stage-1 codebook vectors are increasing and inside (0, 256), and
    /// their IHMW weights lie in the RFC's stated range 1819..=5227.
    #[test]
    fn codebooks_and_weights() {
        for wb in [false, true] {
            let d = order(wb);
            for i1 in 0..32 {
                let cb = cb1(wb, i1);
                assert!(cb[0] > 0 && cb[d - 1] < 256);
                for k in 1..d {
                    assert!(cb[k] > cb[k - 1], "wb {wb} i1 {i1}");
                }
                for &w in &weights_q9(cb, d)[..d] {
                    assert!((1819..=5227).contains(&w), "{w}");
                }
            }
        }
    }

    /// Every PDF of the SILK tables sums to 256 (its inverse CDF ends in 0
    /// and never rises).
    #[test]
    fn icdfs_are_well_formed() {
        let mut all: Vec<&[u8]> = vec![
            &STEREO_STAGE1_ICDF,
            &STEREO_STAGE2_ICDF,
            &STEREO_STAGE3_ICDF,
            &MID_ONLY_ICDF,
            &FRAME_TYPE_INACTIVE_ICDF,
            &FRAME_TYPE_ACTIVE_ICDF,
            &GAIN_LSB_ICDF,
            &GAIN_DELTA_ICDF,
            &NLSF_EXT_ICDF,
            &NLSF_INTERP_ICDF,
            &PITCH_HIGH_ICDF,
            &PITCH_LOW_NB_ICDF,
            &PITCH_LOW_MB_ICDF,
            &PITCH_LOW_WB_ICDF,
            &PITCH_DELTA_ICDF,
            &CONTOUR_NB10_ICDF,
            &CONTOUR_NB20_ICDF,
            &CONTOUR_WB10_ICDF,
            &CONTOUR_WB20_ICDF,
            &PERIODICITY_ICDF,
            &LTP_SCALE_ICDF,
            &SEED_ICDF,
            &LSB_ICDF,
            &LBRR_FLAGS_2_ICDF,
            &LBRR_FLAGS_3_ICDF,
        ];
        all.extend(GAIN_MSB_ICDF.iter().map(|t| &t[..]));
        all.extend(NLSF_STAGE1_ICDF.iter().map(|t| &t[..]));
        all.extend(NLSF_STAGE2_NB_ICDF.iter().map(|t| &t[..]));
        all.extend(NLSF_STAGE2_WB_ICDF.iter().map(|t| &t[..]));
        all.extend(LTP_FILTER_ICDF.iter().copied());
        all.extend(RATE_LEVEL_ICDF.iter().map(|t| &t[..]));
        all.extend(PULSE_COUNT_ICDF.iter().map(|t| &t[..]));
        all.extend(SHELL16_ICDF.iter().copied());
        all.extend(SHELL8_ICDF.iter().copied());
        all.extend(SHELL4_ICDF.iter().copied());
        all.extend(SHELL2_ICDF.iter().copied());
        for t in all {
            assert_eq!(*t.last().unwrap(), 0);
            for w in t.windows(2) {
                assert!(w[0] >= w[1]);
            }
        }
        for (k, t) in SHELL16_ICDF.iter().enumerate() {
            assert_eq!(t.len(), k + 2, "a split of k+1 pulses has k+2 outcomes");
        }
    }

    /// A flat spectrum (equally spaced LSFs) gives a filter close to zero,
    /// and stabilized LSFs always give a stable filter.
    #[test]
    fn lsf_to_lpc_is_stable() {
        for wb in [false, true] {
            let d = order(wb);
            let flat: Vec<i32> = (0..d)
                .map(|k| ((k as i32 + 1) * 32768) / (d as i32 + 1))
                .collect();
            let a = nlsf_to_lpc(&flat, wb);
            assert!(a.iter().all(|v| v.abs() < 200), "{a:?}");
            let mut seed = 99u32;
            for _ in 0..500 {
                let mut n: Vec<i32> = (0..d)
                    .map(|_| {
                        seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                        (seed >> 17) as i32
                    })
                    .collect();
                stabilize(&mut n, wb);
                let dmin: &[i32] = if wb {
                    &NLSF_MIN_SPACING_WB
                } else {
                    &NLSF_MIN_SPACING_NB
                };
                assert!(n[0] >= dmin[0]);
                for k in 1..d {
                    assert!(n[k] - n[k - 1] >= dmin[k], "{n:?}");
                }
                let a = nlsf_to_lpc(&n, wb);
                assert!(is_stable(&a));
            }
        }
    }
}
