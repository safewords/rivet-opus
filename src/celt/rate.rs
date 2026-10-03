//! Bit allocation (CELT_SPEC §5.2, §6, §7.4; RFC 6716 §4.3.3, §5.3.4):
//! band boosts, caps, the reservations, the search over the static
//! allocation table, band skipping, intensity/dual stereo, and the split
//! of each band's budget into fine energy and shape bits. All integer and
//! shared by encoder and decoder through [`Coder`].

use super::mode::mode;
use super::tables::{BAND_ALLOCATION, BITRES, EBANDS, FINE_OFFSET, MAX_FINE_BITS, NB_EBANDS};
use crate::range::Coder;

/// Band width at LM = 0.
fn width(j: usize) -> i32 {
    (EBANDS[j + 1] - EBANDS[j]) as i32
}

/// The caps of CELT_SPEC §6.1: the most eighth bits each band may get.
pub(crate) fn init_caps(lm: usize, c: usize) -> [i32; NB_EBANDS] {
    let row = &mode().cache_caps[2 * lm + c - 1];
    let mut cap = [0i32; NB_EBANDS];
    for (i, v) in cap.iter_mut().enumerate() {
        *v = ((row[i] + 64) * c as i32 * (width(i) << lm)) >> 2;
    }
    cap
}

/// The pseudo-pulse index whose cost is nearest `b` eighth bits, ties to
/// the smaller (CELT_SPEC §7.4; RFC 6716 §4.3.4.1).
pub(crate) fn bits2pulses(band: usize, lm: i32, b: i32) -> usize {
    let cache = mode().cache(band, lm);
    let b = b - 1;
    let mut lo = 0usize;
    let mut hi = usize::from(cache[0]);
    for _ in 0..6 {
        let mid = (lo + hi + 1) >> 1;
        if i32::from(cache[mid]) >= b {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let lo_cost = if lo == 0 { -1 } else { i32::from(cache[lo]) };
    if b - lo_cost <= i32::from(cache[hi]) - b { lo } else { hi }
}

/// The cost in eighth bits of pseudo-pulse index `q` (CELT_SPEC §7.4).
pub(crate) fn pulses2bits(band: usize, lm: i32, q: usize) -> i32 {
    if q == 0 { 0 } else { i32::from(mode().cache(band, lm)[q]) + 1 }
}

/// Band boosts (CELT_SPEC §5.2; RFC 6716 §4.3.3). `want[i]` is how many
/// boost quanta the encoder asks for band `i` (ignored by the decoder).
/// Returns the boosts in eighth bits and the frame budget reduced by them
/// (`total` of §5.2, which also gates the trim).
pub(crate) fn code_boosts<C: Coder>(
    ec: &mut C,
    start: usize,
    end: usize,
    c: usize,
    lm: usize,
    cap: &[i32; NB_EBANDS],
    want: &[i32; NB_EBANDS],
    total_bytes: usize,
) -> ([i32; NB_EBANDS], i32) {
    let mut offsets = [0i32; NB_EBANDS];
    let mut dynalloc_logp = 6u32;
    let mut total = (total_bytes as i32 * 8) << BITRES;
    let mut t = ec.tell_frac();
    for i in start..end {
        let w = c as i32 * (width(i) << lm);
        let quanta = (w << BITRES).min((6 << BITRES).max(w));
        let mut loop_logp = dynalloc_logp;
        let mut boost = 0;
        let mut k = 0;
        while t + ((loop_logp as i32) << BITRES) < total && boost < cap[i] {
            let flag = ec.bit_logp(k < want[i], loop_logp);
            t = ec.tell_frac();
            if !flag {
                break;
            }
            boost += quanta;
            total -= quanta;
            loop_logp = 1;
            k += 1;
        }
        offsets[i] = boost;
        if boost > 0 {
            dynalloc_logp = 2.max(dynalloc_logp - 1);
        }
    }
    (offsets, total)
}

/// The encoder's free choices in the allocation (CELT_SPEC §6.6, §6.7).
#[derive(Clone, Copy, Debug)]
pub(crate) struct EncoderChoices {
    /// The intensity stereo band it would like (clamped to the coded
    /// bands).
    pub intensity: usize,
    /// Whether it would like dual (L/R) stereo.
    pub dual_stereo: bool,
    /// The previous frame's number of coded bands (skip hysteresis).
    pub prev_coded: usize,
}

/// The result of the allocation (CELT_SPEC §6.9 outputs).
#[derive(Clone, Debug)]
pub(crate) struct Allocation {
    /// Shape budget per band in eighth bits.
    pub pulses: [i32; NB_EBANDS],
    /// Fine energy bits per band and channel.
    pub fine_quant: [i32; NB_EBANDS],
    /// Priority of each band for the final fine bits.
    pub fine_priority: [bool; NB_EBANDS],
    /// One past the last band given shape bits.
    pub coded_bands: usize,
    /// The bits left over, carried into the band coding.
    pub balance: i32,
    /// First band coded with intensity stereo.
    pub intensity: usize,
    /// Dual (L/R) stereo.
    pub dual_stereo: bool,
}

/// The allocation of CELT_SPEC §6.2–§6.9 (RFC 6716 §4.3.3) for `bits`
/// eighth bits (after the anti-collapse reservation), coding the skip
/// flags, intensity and dual stereo on the way.
pub(crate) fn compute_allocation<C: Coder>(
    ec: &mut C,
    start: usize,
    end: usize,
    offsets: &[i32; NB_EBANDS],
    cap: &[i32; NB_EBANDS],
    alloc_trim: i32,
    bits: i32,
    c: usize,
    lm: usize,
    choices: EncoderChoices,
) -> Allocation {
    let m = mode();
    let ci = c as i32;
    let lmi = lm as i32;
    // §6.2 reservations.
    let mut total = bits.max(0);
    let skip_rsv = if total >= 1 << BITRES { 1 << BITRES } else { 0 };
    total -= skip_rsv;
    let mut intensity_rsv = 0;
    let mut dual_stereo_rsv = 0;
    if c == 2 {
        intensity_rsv = m.log2_frac_table[end - start];
        if intensity_rsv > total {
            intensity_rsv = 0;
        } else {
            total -= intensity_rsv;
            dual_stereo_rsv = if total >= 1 << BITRES { 1 << BITRES } else { 0 };
            total -= dual_stereo_rsv;
        }
    }
    // §6.3 thresholds and trim offsets.
    let mut thresh = [0i32; NB_EBANDS];
    let mut trim_offset = [0i32; NB_EBANDS];
    for j in start..end {
        let n = width(j);
        thresh[j] = (ci << BITRES).max((3 * (n << lmi) << BITRES) >> 4);
        trim_offset[j] = (ci * n * (alloc_trim - 5 - lmi) * (end - j - 1) as i32 * (1 << (lmi + BITRES))) >> 6;
        if n << lmi == 1 {
            trim_offset[j] -= ci << BITRES;
        }
    }
    // §6.4 the bracketing rows of the static table.
    let adj = |q: usize, j: usize| {
        let raw = (ci * width(j) * i32::from(BAND_ALLOCATION[q][j]) << lmi) >> 2;
        if raw > 0 { (raw + trim_offset[j]).max(0) } else { raw }
    };
    let mut lo = 1usize;
    let mut hi = 10usize;
    while lo <= hi {
        let mid = (lo + hi) >> 1;
        let mut psum = 0;
        let mut done = false;
        for j in (start..end).rev() {
            let v = adj(mid, j) + offsets[j];
            if v >= thresh[j] || done {
                done = true;
                psum += v.min(cap[j]);
            } else if v >= ci << BITRES {
                psum += ci << BITRES;
            }
        }
        if psum > total {
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    hi = lo;
    lo -= 1;
    let mut bits1 = [0i32; NB_EBANDS];
    let mut bits2 = [0i32; NB_EBANDS];
    let mut skip_start = start;
    for j in start..end {
        let mut b1 = adj(lo, j);
        let mut b2 = if hi <= 10 { adj(hi, j) } else { (cap[j] + trim_offset[j]).max(0) };
        if lo > 0 {
            b1 += offsets[j];
        }
        b2 += offsets[j];
        if offsets[j] > 0 {
            skip_start = j;
        }
        bits1[j] = b1;
        bits2[j] = (b2 - b1).max(0);
    }
    // §6.5 interpolation in 1/64 steps.
    let alloc_floor = ci << BITRES;
    let (mut lo, mut hi) = (0i32, 64i32);
    for _ in 0..6 {
        let mid = (lo + hi) >> 1;
        let mut psum = 0;
        let mut done = false;
        for j in (start..end).rev() {
            let v = bits1[j] + ((mid * bits2[j]) >> 6);
            if v >= thresh[j] || done {
                done = true;
                psum += v.min(cap[j]);
            } else if v >= alloc_floor {
                psum += alloc_floor;
            }
        }
        if psum > total {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let mut psum = 0;
    let mut done = false;
    let mut alloc = [0i32; NB_EBANDS];
    for j in (start..end).rev() {
        let mut v = bits1[j] + ((lo * bits2[j]) >> 6);
        if v < thresh[j] && !done {
            v = if v >= alloc_floor { alloc_floor } else { 0 };
        } else {
            done = true;
        }
        v = v.min(cap[j]);
        alloc[j] = v;
        psum += v;
    }
    // §6.6 skipping.
    let mut coded = end;
    loop {
        let j = coded - 1;
        if j <= skip_start {
            total += skip_rsv;
            break;
        }
        let mut left = total - psum;
        let span = (EBANDS[coded] - EBANDS[start]) as i32;
        let percoeff = left / span;
        left -= span * percoeff;
        let rem = (left - (EBANDS[j] - EBANDS[start]) as i32).max(0);
        let band_width = (EBANDS[coded] - EBANDS[j]) as i32;
        let mut band_bits = alloc[j] + percoeff * band_width + rem;
        if band_bits >= thresh[j].max(alloc_floor + (1 << BITRES)) {
            let keep = C::ENCODE && {
                let k = if j < choices.prev_coded { 7 } else { 9 };
                band_bits > ((k * band_width) << lmi << BITRES) >> 4
            };
            if ec.bit_logp(keep, 1) {
                break;
            }
            psum += 1 << BITRES;
            band_bits -= 1 << BITRES;
        }
        psum -= alloc[j] + intensity_rsv;
        if intensity_rsv > 0 {
            intensity_rsv = m.log2_frac_table[j - start];
        }
        psum += intensity_rsv;
        if band_bits >= alloc_floor {
            psum += alloc_floor;
            alloc[j] = alloc_floor;
        } else {
            alloc[j] = 0;
        }
        coded -= 1;
    }
    // §6.7 intensity and dual stereo.
    let intensity = if intensity_rsv > 0 {
        let want = choices.intensity.clamp(start, coded);
        start + ec.uint((want - start) as u32, (coded + 1 - start) as u32) as usize
    } else {
        0
    };
    if intensity <= start {
        total += dual_stereo_rsv;
        dual_stereo_rsv = 0;
    }
    let dual_stereo = dual_stereo_rsv > 0 && ec.bit_logp(choices.dual_stereo, 1);
    // §6.8 the remainder.
    let mut left = total - psum;
    let span = (EBANDS[coded] - EBANDS[start]) as i32;
    let percoeff = left / span;
    left -= span * percoeff;
    for j in start..coded {
        alloc[j] += percoeff * width(j);
    }
    for j in start..coded {
        let t = left.min(width(j));
        alloc[j] += t;
        left -= t;
    }
    // §6.9 fine energy / shape split.
    let mut fine_quant = [0i32; NB_EBANDS];
    let mut fine_priority = [false; NB_EBANDS];
    let mut balance = 0;
    let log_m = lmi << BITRES;
    for j in start..coded {
        let n = width(j) << lmi;
        alloc[j] += balance;
        let mut excess;
        if n > 1 {
            excess = (alloc[j] - cap[j]).max(0);
            alloc[j] -= excess;
            let den = ci * n + i32::from(c == 2 && n > 2 && !dual_stereo && j < intensity);
            let nclogn = den * (m.log_n[j] + log_m);
            let mut offset = (nclogn >> 1) - den * FINE_OFFSET;
            if n == 2 {
                offset += den << BITRES >> 2;
            }
            if alloc[j] + offset < den * 2 << BITRES {
                offset += nclogn >> 2;
            } else if alloc[j] + offset < den * 3 << BITRES {
                offset += nclogn >> 3;
            }
            let mut eb = ((alloc[j] + offset + (den << (BITRES - 1))) / (den << BITRES)).max(0);
            if ci * eb > alloc[j] >> BITRES {
                eb = alloc[j] >> (c - 1) >> BITRES;
            }
            eb = eb.min(MAX_FINE_BITS);
            fine_priority[j] = eb * (den << BITRES) >= alloc[j] + offset;
            alloc[j] -= ci * eb << BITRES;
            fine_quant[j] = eb;
        } else {
            excess = (alloc[j] - (ci << BITRES)).max(0);
            alloc[j] -= excess;
            fine_quant[j] = 0;
            fine_priority[j] = true;
        }
        if excess > 0 {
            let extra_fine = (excess >> (c - 1 + BITRES as usize)).min(MAX_FINE_BITS - fine_quant[j]);
            fine_quant[j] += extra_fine;
            let extra_bits = extra_fine * ci << BITRES;
            fine_priority[j] = extra_bits >= excess - balance;
            excess -= extra_bits;
        }
        balance = excess;
    }
    for j in coded..end {
        fine_quant[j] = alloc[j] >> (c - 1) >> BITRES;
        alloc[j] = 0;
        fine_priority[j] = fine_quant[j] < 1;
    }
    Allocation { pulses: alloc, fine_quant, fine_priority, coded_bands: coded, balance, intensity, dual_stereo }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bits2pulses_picks_the_nearest_cost() {
        for lm in -1..=3i32 {
            for band in 0..NB_EBANDS {
                if mode().cache_index[(lm + 1) as usize][band] < 0 {
                    continue;
                }
                let kmax = usize::from(mode().cache(band, lm)[0]);
                for q in 0..=kmax {
                    let cost = pulses2bits(band, lm, q);
                    // Equal costs (N = 1) resolve to the smallest count.
                    let back = bits2pulses(band, lm, cost);
                    assert!(back <= q && pulses2bits(band, lm, back) == cost, "band {band} lm {lm} q {q}");
                }
            }
        }
    }
}
