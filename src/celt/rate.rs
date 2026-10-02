//! The CELT bit allocation (RFC 6716 §4.3.3): the search of the static
//! allocation table, the interpolation between its rows, band skipping,
//! the intensity and dual stereo parameters, and the split of each band's
//! budget between fine energy and shape. Both sides run this code; it must
//! give identical results in encoder and decoder.

use super::mode::{BITRES, mode};
use super::tables::{ALLOC, EBANDS, NB_EBANDS};
use crate::range::Coder;

const ALLOC_STEPS: i32 = 6;
/// Fine energy is capped at this many bits per band and channel.
pub const MAX_FINE_BITS: i32 = 8;
const FINE_OFFSET: i32 = 21;

/// `ceil(8 * log2(i))`, the conservative cost of coding one of `i` values
/// (`LOG2_FRAC_TABLE`), for `i <= 24`.
pub fn log2_frac_table(i: usize) -> i32 {
    if i <= 1 { 0 } else { super::mode::log2_frac(i as u32, BITRES) }
}

/// The allocation of one frame.
pub struct Allocation {
    /// Shape bits per band, 1/8 bit.
    pub pulses: [i32; NB_EBANDS],
    /// Fine energy bits per band and channel.
    pub fine_quant: [i32; NB_EBANDS],
    /// Which bands get the leftover bits first (0) or second (1).
    pub fine_priority: [i32; NB_EBANDS],
    /// The first band with no shape bits.
    pub coded_bands: usize,
    /// Bits left over, 1/8 bit.
    pub balance: i32,
    /// First intensity-stereo band.
    pub intensity: usize,
    /// Dual (L/R) stereo.
    pub dual_stereo: bool,
}

/// Inputs the encoder decides (ignored by the decoder).
#[derive(Clone, Copy, Default)]
pub struct EncoderChoices {
    pub intensity: usize,
    pub dual_stereo: bool,
    /// Band skipping hysteresis: the previous frame's coded band count.
    pub prev_coded: usize,
    /// The highest band with signal.
    pub signal_bandwidth: usize,
}

#[allow(clippy::too_many_arguments)]
pub fn compute_allocation<E: Coder>(
    start: usize,
    end: usize,
    offsets: &[i32; NB_EBANDS],
    cap: &[i32; NB_EBANDS],
    alloc_trim: i32,
    total: i32,
    c: usize,
    lm: usize,
    ec: &mut E,
    choices: EncoderChoices,
) -> Allocation {
    let ci = c as i32;
    let mut total = total.max(0);
    let mut skip_start = start;
    let skip_rsv = if total >= 1 << BITRES { 1 << BITRES } else { 0 };
    total -= skip_rsv;
    let mut intensity_rsv = 0;
    let mut dual_stereo_rsv = 0;
    if c == 2 {
        intensity_rsv = log2_frac_table(end - start);
        if intensity_rsv > total {
            intensity_rsv = 0;
        } else {
            total -= intensity_rsv;
            dual_stereo_rsv = if total >= 1 << BITRES { 1 << BITRES } else { 0 };
            total -= dual_stereo_rsv;
        }
    }
    let mut thresh = [0i32; NB_EBANDS];
    let mut trim_offset = [0i32; NB_EBANDS];
    let width = |j: usize| (EBANDS[j + 1] - EBANDS[j]) as i32;
    for j in start..end {
        thresh[j] = (ci << BITRES).max((3 * width(j)) << lm << BITRES >> 4);
        trim_offset[j] =
            ci * width(j) * (alloc_trim - 5 - lm as i32) * (end - j - 1) as i32 * (1 << (lm as i32 + BITRES)) >> 6;
        if width(j) << lm == 1 {
            trim_offset[j] -= ci << BITRES;
        }
    }
    let alloc_bits = |q: usize, j: usize| (ci * width(j) * i32::from(ALLOC[q][j])) << lm >> 2;
    let mut lo = 1usize;
    let mut hi = ALLOC.len() - 1;
    while lo <= hi {
        let mid = (lo + hi) >> 1;
        let mut done = false;
        let mut psum = 0;
        for j in (start..end).rev() {
            let mut bitsj = alloc_bits(mid, j);
            if bitsj > 0 {
                bitsj = (bitsj + trim_offset[j]).max(0);
            }
            bitsj += offsets[j];
            if bitsj >= thresh[j] || done {
                done = true;
                psum += bitsj.min(cap[j]);
            } else if bitsj >= ci << BITRES {
                psum += ci << BITRES;
            }
        }
        if psum > total {
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    let hi = lo;
    let lo = lo - 1;
    let mut bits1 = [0i32; NB_EBANDS];
    let mut bits2 = [0i32; NB_EBANDS];
    for j in start..end {
        let mut b1 = alloc_bits(lo, j);
        let mut b2 = if hi >= ALLOC.len() { cap[j] } else { alloc_bits(hi, j) };
        if b1 > 0 {
            b1 = (b1 + trim_offset[j]).max(0);
        }
        if b2 > 0 {
            b2 = (b2 + trim_offset[j]).max(0);
        }
        if lo > 0 {
            b1 += offsets[j];
        }
        b2 += offsets[j];
        if offsets[j] > 0 {
            skip_start = j;
        }
        b2 = (b2 - b1).max(0);
        bits1[j] = b1;
        bits2[j] = b2;
    }
    interp_bits2pulses(
        start,
        end,
        skip_start,
        &bits1,
        &bits2,
        &thresh,
        cap,
        total,
        skip_rsv,
        intensity_rsv,
        dual_stereo_rsv,
        c,
        lm,
        ec,
        choices,
    )
}

#[allow(clippy::too_many_arguments)]
fn interp_bits2pulses<E: Coder>(
    start: usize,
    end: usize,
    skip_start: usize,
    bits1: &[i32; NB_EBANDS],
    bits2: &[i32; NB_EBANDS],
    thresh: &[i32; NB_EBANDS],
    cap: &[i32; NB_EBANDS],
    mut total: i32,
    skip_rsv: i32,
    mut intensity_rsv: i32,
    mut dual_stereo_rsv: i32,
    c: usize,
    lm: usize,
    ec: &mut E,
    choices: EncoderChoices,
) -> Allocation {
    let ci = c as i32;
    let alloc_floor = ci << BITRES;
    let stereo = i32::from(c > 1);
    let log_m = (lm as i32) << BITRES;
    let mut lo = 0;
    let mut hi = 1 << ALLOC_STEPS;
    for _ in 0..ALLOC_STEPS {
        let mid = (lo + hi) >> 1;
        let mut psum = 0;
        let mut done = false;
        for j in (start..end).rev() {
            let tmp = bits1[j] + (mid * bits2[j] >> ALLOC_STEPS);
            if tmp >= thresh[j] || done {
                done = true;
                psum += tmp.min(cap[j]);
            } else if tmp >= alloc_floor {
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
    let mut bits = [0i32; NB_EBANDS];
    for j in (start..end).rev() {
        let mut tmp = bits1[j] + (lo * bits2[j] >> ALLOC_STEPS);
        if tmp < thresh[j] && !done {
            tmp = if tmp >= alloc_floor { alloc_floor } else { 0 };
        } else {
            done = true;
        }
        tmp = tmp.min(cap[j]);
        bits[j] = tmp;
        psum += tmp;
    }
    // Skip bands from the top while that frees bits.
    let eb = |j: usize| EBANDS[j] as i32;
    let mut coded_bands = end;
    loop {
        let j = coded_bands - 1;
        if j <= skip_start {
            total += skip_rsv;
            break;
        }
        let mut left = total - psum;
        let percoeff = left / (eb(coded_bands) - eb(start));
        left -= (eb(coded_bands) - eb(start)) * percoeff;
        let rem = (left - (eb(j) - eb(start))).max(0);
        let band_width = eb(coded_bands) - eb(j);
        let mut band_bits = bits[j] + percoeff * band_width + rem;
        if band_bits >= thresh[j].max(alloc_floor + (1 << BITRES)) {
            let keep = if E::ENCODE {
                let depth_ok = band_bits > ((if j < choices.prev_coded { 7 } else { 9 }) * band_width << lm << BITRES) >> 4;
                coded_bands <= start + 2 || (depth_ok && j <= choices.signal_bandwidth)
            } else {
                false
            };
            if ec.bit_logp(keep, 1) {
                break;
            }
            psum += 1 << BITRES;
            band_bits -= 1 << BITRES;
        }
        psum -= bits[j] + intensity_rsv;
        if intensity_rsv > 0 {
            intensity_rsv = log2_frac_table(j - start);
        }
        psum += intensity_rsv;
        if band_bits >= alloc_floor {
            psum += alloc_floor;
            bits[j] = alloc_floor;
        } else {
            bits[j] = 0;
        }
        coded_bands -= 1;
    }
    let mut intensity = 0;
    if intensity_rsv > 0 {
        let v = choices.intensity.min(coded_bands).max(start);
        intensity = start + ec.uint((v - start) as u32, (coded_bands + 1 - start) as u32) as usize;
    }
    if intensity <= start {
        total += dual_stereo_rsv;
        dual_stereo_rsv = 0;
    }
    let dual_stereo = if dual_stereo_rsv > 0 { ec.bit_logp(choices.dual_stereo, 1) } else { false };
    // Spread what is left over the coded bands.
    let mut left = total - psum;
    let percoeff = left / (eb(coded_bands) - eb(start));
    left -= (eb(coded_bands) - eb(start)) * percoeff;
    for j in start..coded_bands {
        bits[j] += percoeff * (eb(j + 1) - eb(j));
    }
    for j in start..coded_bands {
        let tmp = left.min(eb(j + 1) - eb(j));
        bits[j] += tmp;
        left -= tmp;
    }
    let m = mode();
    let mut fine_quant = [0i32; NB_EBANDS];
    let mut fine_priority = [0i32; NB_EBANDS];
    let mut balance = 0;
    for j in start..coded_bands {
        let n0 = eb(j + 1) - eb(j);
        let n = n0 << lm;
        let bit = bits[j] + balance;
        let mut excess;
        if n > 1 {
            excess = (bit - cap[j]).max(0);
            bits[j] = bit - excess;
            let den = ci * n + i32::from(c == 2 && n > 2 && !dual_stereo && j < intensity);
            let nclogn = den * (m.log_n[j] + log_m);
            let mut offset = (nclogn >> 1) - den * FINE_OFFSET;
            if n == 2 {
                offset += den << BITRES >> 2;
            }
            if bits[j] + offset < den * 2 << BITRES {
                offset += nclogn >> 2;
            } else if bits[j] + offset < den * 3 << BITRES {
                offset += nclogn >> 3;
            }
            let mut eb_j = (bits[j] + offset + (den << (BITRES - 1))).max(0);
            eb_j = (eb_j / den) >> BITRES;
            if ci * eb_j > (bits[j] >> BITRES) {
                eb_j = bits[j] >> stereo >> BITRES;
            }
            eb_j = eb_j.min(MAX_FINE_BITS);
            fine_priority[j] = i32::from(eb_j * (den << BITRES) >= bits[j] + offset);
            bits[j] -= ci * eb_j << BITRES;
            fine_quant[j] = eb_j;
        } else {
            excess = (bit - (ci << BITRES)).max(0);
            bits[j] = bit - excess;
            fine_quant[j] = 0;
            fine_priority[j] = 1;
        }
        if excess > 0 {
            let extra_fine = (excess >> (stereo + BITRES)).min(MAX_FINE_BITS - fine_quant[j]);
            fine_quant[j] += extra_fine;
            let extra_bits = extra_fine * ci << BITRES;
            fine_priority[j] = i32::from(extra_bits >= excess - balance);
            excess -= extra_bits;
        }
        balance = excess;
    }
    for j in coded_bands..end {
        fine_quant[j] = bits[j] >> stereo >> BITRES;
        bits[j] = 0;
        fine_priority[j] = i32::from(fine_quant[j] < 1);
    }
    Allocation { pulses: bits, fine_quant, fine_priority, coded_bands, balance, intensity, dual_stereo }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log2_frac_table_values() {
        let expect = [0, 0, 8, 13, 16, 19, 21, 23, 24, 26, 27, 28, 29, 30, 31, 32, 32, 33, 34, 34, 35, 36, 36, 37];
        for (i, &e) in expect.iter().enumerate() {
            assert_eq!(log2_frac_table(i), e, "{i}");
        }
    }
}
