//! The derived tables of the CELT mode: everything CELT_SPEC gives a
//! construction for is computed here once (pulse cache §7.2, caps §7.3,
//! `log2_frac` §7.1 and the tables built on it, the interleaving orders
//! §8.2.2, `exp2_table8` §8.2.4). The unit tests check each construction
//! against the data extracted from RFC 6716 Appendix A.

use std::sync::OnceLock;

use super::cwrs;
use super::tables::{EBANDS, MAX_PSEUDO, NB_EBANDS, WINDOW120};
use crate::range::ilog;

/// The derived tables of the 48 kHz CELT mode.
pub(crate) struct Mode {
    /// The rising half of the low-overlap window (CELT_SPEC §10.4).
    pub window: Vec<f32>,
    /// `LOG2_FRAC_TABLE[i] = log2_frac(i + 1, 3)` (CELT_SPEC §6.2).
    pub log2_frac_table: [i32; 24],
    /// `logN[j] = log2_frac(width_j, 3)` (CELT_SPEC §6.9).
    pub log_n: [i32; NB_EBANDS],
    /// Offsets of the pulse-cache entries, by row `LM + 1` and band
    /// (−1: no entry) (CELT_SPEC §7.2).
    pub cache_index: [[i32; NB_EBANDS]; 5],
    /// The pulse-cache entries: `Kmax`, then the cost of each pseudo-pulse
    /// count (CELT_SPEC §7.2).
    pub cache_bits: Vec<u8>,
    /// The band caps by row `2·LM + C − 1` (CELT_SPEC §7.3).
    pub cache_caps: [[i32; NB_EBANDS]; 8],
    /// The Hadamard ("sequency") block orders for strides 2, 4, 8 and 16
    /// concatenated (CELT_SPEC §8.2.2).
    pub ordery: [usize; 30],
    /// Collapse-mask bit maps of the TF recombination (CELT_SPEC §8.2.2).
    pub bit_interleave: [u32; 16],
    pub bit_deinterleave: [u32; 16],
    /// `⌊16384·2^(k/8)⌋` (CELT_SPEC §8.2.4).
    pub exp2_table8: [i32; 8],
}

/// The mode's tables, built on first use.
pub(crate) fn mode() -> &'static Mode {
    static MODE: OnceLock<Mode> = OnceLock::new();
    MODE.get_or_init(Mode::build)
}

/// The conservative fixed-point log2 of CELT_SPEC §7.1: `frac` fractional
/// bits, never below the true value. `val > 0`.
pub(crate) fn log2_frac(val: u32, frac: u32) -> i32 {
    let l = ilog(val);
    if val & (val - 1) == 0 {
        return (l - 1) << frac;
    }
    // Normalise to 16 fractional bits, rounding up.
    let mut v: u32 = if l > 16 {
        let sh = (l - 16) as u32;
        (val >> sh) + u32::from(val & ((1 << sh) - 1) != 0)
    } else {
        val << (16 - l)
    };
    let mut r = (l - 1) << frac;
    for f in (0..=frac).rev() {
        let b = v >> 16;
        r += (b << f) as i32;
        v = (v + b) >> b;
        v = (v * v + 0x7FFF) >> 15;
    }
    r + i32::from(v > 0x8000)
}

/// The pulse count of pseudo-pulse index `q` (CELT_SPEC §7.2).
pub(crate) fn get_pulses(q: usize) -> usize {
    if q < 8 {
        q
    } else {
        (8 + (q & 7)) << ((q >> 3) - 1)
    }
}

/// Band width at LM = 0.
fn width(j: usize) -> usize {
    EBANDS[j + 1] - EBANDS[j]
}

impl Mode {
    fn build() -> Self {
        let mut log2_frac_table = [0i32; 24];
        for (i, v) in log2_frac_table.iter_mut().enumerate() {
            *v = log2_frac(i as u32 + 1, 3);
        }
        let mut log_n = [0i32; NB_EBANDS];
        for (j, v) in log_n.iter_mut().enumerate() {
            *v = log2_frac(width(j) as u32, 3);
        }
        let (cache_index, cache_bits) = build_cache();
        let cache_caps = build_caps(&cache_index, &cache_bits, &log_n);
        let mut ordery = [0usize; 30];
        let mut pos = 0;
        for s in 1..=4u32 {
            let stride = 1usize << s;
            for i in 0..stride {
                let rev = (i as u32).reverse_bits() >> (32 - s);
                ordery[pos] = stride - 1 - inverse_gray(rev) as usize;
                pos += 1;
            }
        }
        let mut bit_interleave = [0u32; 16];
        let mut bit_deinterleave = [0u32; 16];
        for x in 0..16u32 {
            bit_interleave[x as usize] = u32::from(x & 3 != 0) | (u32::from(x & 12 != 0) << 1);
            bit_deinterleave[x as usize] = (0..4)
                .filter(|k| x >> k & 1 != 0)
                .map(|k| 3 << (2 * k))
                .sum();
        }
        let mut exp2_table8 = [0i32; 8];
        for (k, v) in exp2_table8.iter_mut().enumerate() {
            *v = (16384.0 * (k as f64 / 8.0).exp2()).floor() as i32;
        }
        Self {
            window: WINDOW120.to_vec(),
            log2_frac_table,
            log_n,
            cache_index,
            cache_bits,
            cache_caps,
            ordery,
            bit_interleave,
            bit_deinterleave,
            exp2_table8,
        }
    }

    /// The pulse-cache entry of band `band` at `lm` (−1 ..= 3), `entry[0]`
    /// being `Kmax` (CELT_SPEC §7.2).
    pub fn cache(&self, band: usize, lm: i32) -> &[u8] {
        let o = self.cache_index[(lm + 1) as usize][band];
        debug_assert!(o >= 0);
        let o = o as usize;
        &self.cache_bits[o..=o + usize::from(self.cache_bits[o])]
    }
}

/// `x XOR x>>1 XOR x>>2 XOR …`, the inverse of the Gray code.
fn inverse_gray(mut x: u32) -> u32 {
    let mut r = 0;
    while x != 0 {
        r ^= x;
        x >>= 1;
    }
    r
}

/// The pulse cache of CELT_SPEC §7.2: one entry per distinct band size.
fn build_cache() -> ([[i32; NB_EBANDS]; 5], Vec<u8>) {
    let mut index = [[-1i32; NB_EBANDS]; 5];
    let mut bits = Vec::new();
    let mut seen: Vec<(usize, i32)> = Vec::new();
    for (r, row) in index.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            let n = (width(j) << r) >> 1;
            if n == 0 {
                continue;
            }
            if let Some(&(_, o)) = seen.iter().find(|(m, _)| *m == n) {
                *slot = o;
                continue;
            }
            let o = bits.len() as i32;
            seen.push((n, o));
            *slot = o;
            let kmax = (0..=MAX_PSEUDO)
                .rev()
                .find(|&q| cwrs::v(n, get_pulses(q)) < 1 << 32)
                .unwrap_or(0);
            bits.push(kmax as u8);
            for q in 1..=kmax {
                let v = cwrs::v(n, get_pulses(q)) as u32;
                bits.push((log2_frac(v, 3) - 1) as u8);
            }
        }
    }
    (index, bits)
}

/// The caps of CELT_SPEC §7.3: the most eighth bits a band can use at each
/// LM and channel count, in the byte form the allocation reads.
fn build_caps(
    index: &[[i32; NB_EBANDS]; 5],
    bits: &[u8],
    log_n: &[i32; NB_EBANDS],
) -> [[i32; NB_EBANDS]; 8] {
    let mut caps = [[0i32; NB_EBANDS]; 8];
    for lm in 0..4i32 {
        for c in 1..=2i32 {
            for j in 0..NB_EBANDS {
                let w = width(j) as i32;
                let max_bits = if w << lm == 1 {
                    8 * c * (1 + 8)
                } else {
                    let mut n0 = w;
                    let mut lm0 = 0i32;
                    if n0 > 2 {
                        n0 /= 2;
                        lm0 = -1;
                    } else if n0 <= 1 {
                        lm0 = lm.min(1);
                        n0 <<= lm0;
                    }
                    let o = index[(lm0 + 1) as usize][j] as usize;
                    let mut max_bits = i32::from(bits[o + usize::from(bits[o])]) + 1;
                    let mut n = n0;
                    for k in 0..lm - lm0 {
                        max_bits *= 2;
                        let offset = ((log_n[j] + 8 * (lm0 + k)) >> 1) - 4;
                        let num = 459 * ((2 * n - 1) * offset + max_bits);
                        let den = ((2 * n - 1) << 9) - 459;
                        max_bits += ((num + den / 2) / den).min(57);
                        n *= 2;
                    }
                    if c == 2 {
                        max_bits *= 2;
                        let two = n == 2;
                        let offset = ((log_n[j] + 8 * lm) >> 1) - if two { 16 } else { 4 };
                        let ndof = 2 * n - 1 - i32::from(two);
                        let f = if two { 512 } else { 487 };
                        let num = f * (max_bits + ndof * offset);
                        let den = (ndof << 9) - f;
                        max_bits += ((num + den / 2) / den).min(if two { 64 } else { 61 });
                    }
                    let ndof = c * n + i32::from(c == 2 && n > 2);
                    let offset = ((log_n[j] + 8 * lm) >> 1) - 21 + if n == 2 { 2 } else { 0 };
                    let num = max_bits + ndof * offset;
                    let den = (ndof - 1) << 3;
                    max_bits + 8 * c * ((num + den / 2) / den).min(8)
                };
                caps[(2 * lm + c - 1) as usize][j] = 4 * max_bits / (c * (w << lm)) - 64;
            }
        }
    }
    caps
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::celt::tables::{BETA_COEF, BETA_INTRA, E_MEANS, PRED_COEF};

    /// The data RFC 6716 Appendix A holds for the derived tables, as printed
    /// by tools/appendix_a_tables.py, to check the derivations against.
    mod extracted {
        /// RFC 6716 Appendix A celt/rate.c, `LOG2_FRAC_TABLE`, extracted by tools/appendix_a_tables.py.
        pub const LOG2_FRAC_TABLE: [u8; 24] = [
            0, 8, 13, 16, 19, 21, 23, 24, 26, 27, 28, 29, 30, 31, 32, 32, 33, 34, 34, 35, 36, 36,
            37, 37,
        ];

        /// RFC 6716 Appendix A celt/bands.c, `ordery_table`, extracted by tools/appendix_a_tables.py.
        pub const ORDERY_TABLE: [i32; 30] = [
            1, 0, 3, 0, 2, 1, 7, 0, 4, 3, 6, 1, 5, 2, 15, 0, 8, 7, 12, 3, 11, 4, 14, 1, 9, 6, 13,
            2, 10, 5,
        ];

        /// RFC 6716 Appendix A celt/bands.c, `bit_interleave_table`, extracted by tools/appendix_a_tables.py.
        pub const BIT_INTERLEAVE_TABLE: [u8; 16] = [0, 1, 1, 1, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3];

        /// RFC 6716 Appendix A celt/bands.c, `bit_deinterleave_table`, extracted by tools/appendix_a_tables.py.
        pub const BIT_DEINTERLEAVE_TABLE: [u8; 16] = [
            0, 3, 12, 15, 48, 51, 60, 63, 192, 195, 204, 207, 240, 243, 252, 255,
        ];

        /// RFC 6716 Appendix A celt/bands.c, `exp2_table8`, extracted by tools/appendix_a_tables.py.
        pub const EXP2_TABLE8: [i16; 8] = [16384, 17866, 19483, 21247, 23170, 25267, 27554, 30048];

        /// RFC 6716 Appendix A celt/static_modes_float.h, `logN400`, extracted by tools/appendix_a_tables.py.
        pub const LOGN400: [i16; 21] = [
            0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 16, 16, 16, 21, 21, 24, 29, 34, 36,
        ];

        /// RFC 6716 Appendix A celt/static_modes_float.h, `cache_index50`, extracted by tools/appendix_a_tables.py.
        pub const CACHE_INDEX50: [[i16; 21]; 5] = [
            [
                -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 41, 41, 41, 82, 82, 123, 164, 200, 222,
            ],
            [
                0, 0, 0, 0, 0, 0, 0, 0, 41, 41, 41, 41, 123, 123, 123, 164, 164, 240, 266, 283, 295,
            ],
            [
                41, 41, 41, 41, 41, 41, 41, 41, 123, 123, 123, 123, 240, 240, 240, 266, 266, 305,
                318, 328, 336,
            ],
            [
                123, 123, 123, 123, 123, 123, 123, 123, 240, 240, 240, 240, 305, 305, 305, 318,
                318, 343, 351, 358, 364,
            ],
            [
                240, 240, 240, 240, 240, 240, 240, 240, 305, 305, 305, 305, 343, 343, 343, 351,
                351, 370, 376, 382, 387,
            ],
        ];

        /// RFC 6716 Appendix A celt/static_modes_float.h, `cache_bits50`, extracted by tools/appendix_a_tables.py.
        #[rustfmt::skip]
        pub const CACHE_BITS50: [u8; 392] = [
            40, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7, 40, 15, 23, 28, 31, 34, 36,
            38, 39, 41, 42, 43, 44, 45, 46, 47, 47, 49, 50,
            51, 52, 53, 54, 55, 55, 57, 58, 59, 60, 61, 62,
            63, 63, 65, 66, 67, 68, 69, 70, 71, 71, 40, 20,
            33, 41, 48, 53, 57, 61, 64, 66, 69, 71, 73, 75,
            76, 78, 80, 82, 85, 87, 89, 91, 92, 94, 96, 98,
            101, 103, 105, 107, 108, 110, 112, 114, 117, 119, 121, 123,
            124, 126, 128, 40, 23, 39, 51, 60, 67, 73, 79, 83,
            87, 91, 94, 97, 100, 102, 105, 107, 111, 115, 118, 121,
            124, 126, 129, 131, 135, 139, 142, 145, 148, 150, 153, 155,
            159, 163, 166, 169, 172, 174, 177, 179, 35, 28, 49, 65,
            78, 89, 99, 107, 114, 120, 126, 132, 136, 141, 145, 149,
            153, 159, 165, 171, 176, 180, 185, 189, 192, 199, 205, 211,
            216, 220, 225, 229, 232, 239, 245, 251, 21, 33, 58, 79,
            97, 112, 125, 137, 148, 157, 166, 174, 182, 189, 195, 201,
            207, 217, 227, 235, 243, 251, 17, 35, 63, 86, 106, 123,
            139, 152, 165, 177, 187, 197, 206, 214, 222, 230, 237, 250,
            25, 31, 55, 75, 91, 105, 117, 128, 138, 146, 154, 161,
            168, 174, 180, 185, 190, 200, 208, 215, 222, 229, 235, 240,
            245, 255, 16, 36, 65, 89, 110, 128, 144, 159, 173, 185,
            196, 207, 217, 226, 234, 242, 250, 11, 41, 74, 103, 128,
            151, 172, 191, 209, 225, 241, 255, 9, 43, 79, 110, 138,
            163, 186, 207, 227, 246, 12, 39, 71, 99, 123, 144, 164,
            182, 198, 214, 228, 241, 253, 9, 44, 81, 113, 142, 168,
            192, 214, 235, 255, 7, 49, 90, 127, 160, 191, 220, 247,
            6, 51, 95, 134, 170, 203, 234, 7, 47, 87, 123, 155,
            184, 212, 237, 6, 52, 97, 137, 174, 208, 240, 5, 57,
            106, 151, 192, 231, 5, 59, 111, 158, 202, 243, 5, 55,
            103, 147, 187, 224, 5, 60, 113, 161, 206, 248, 4, 65,
            122, 175, 224, 4, 67, 127, 182, 234,
        ];

        /// RFC 6716 Appendix A celt/static_modes_float.h, `cache_caps50`, extracted by tools/appendix_a_tables.py.
        pub const CACHE_CAPS50: [[u8; 21]; 8] = [
            [
                224, 224, 224, 224, 224, 224, 224, 224, 160, 160, 160, 160, 185, 185, 185, 178,
                178, 168, 134, 61, 37,
            ],
            [
                224, 224, 224, 224, 224, 224, 224, 224, 240, 240, 240, 240, 207, 207, 207, 198,
                198, 183, 144, 66, 40,
            ],
            [
                160, 160, 160, 160, 160, 160, 160, 160, 185, 185, 185, 185, 193, 193, 193, 183,
                183, 172, 138, 64, 38,
            ],
            [
                240, 240, 240, 240, 240, 240, 240, 240, 207, 207, 207, 207, 204, 204, 204, 193,
                193, 180, 143, 66, 40,
            ],
            [
                185, 185, 185, 185, 185, 185, 185, 185, 193, 193, 193, 193, 193, 193, 193, 183,
                183, 172, 138, 65, 39,
            ],
            [
                207, 207, 207, 207, 207, 207, 207, 207, 204, 204, 204, 204, 201, 201, 201, 188,
                188, 176, 141, 66, 40,
            ],
            [
                193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 194, 194, 194, 184,
                184, 173, 139, 65, 39,
            ],
            [
                204, 204, 204, 204, 204, 204, 204, 204, 201, 201, 201, 201, 198, 198, 198, 187,
                187, 175, 140, 66, 40,
            ],
        ];

        /// RFC 6716 Appendix A celt/quant_bands.c, `eMeans` (FIXED_POINT branch), extracted by tools/appendix_a_tables.py.
        pub const E_MEANS_Q4: [i8; 25] = [
            103, 100, 92, 85, 81, 77, 72, 70, 78, 75, 73, 71, 78, 74, 69, 72, 70, 74, 76, 71, 60,
            60, 60, 60, 60,
        ];

        /// RFC 6716 Appendix A celt/quant_bands.c, `pred_coef` (FIXED_POINT branch), extracted by tools/appendix_a_tables.py.
        pub const PRED_COEF_Q15: [i16; 4] = [29440, 26112, 21248, 16384];

        /// RFC 6716 Appendix A celt/quant_bands.c, `beta_coef` (FIXED_POINT branch), extracted by tools/appendix_a_tables.py.
        pub const BETA_COEF_Q15: [i16; 4] = [30147, 22282, 12124, 6554];

        /// RFC 6716 Appendix A celt/quant_bands.c, `beta_intra` (FIXED_POINT branch), extracted by tools/appendix_a_tables.py.
        pub const BETA_INTRA_Q15: i16 = 4915;
    }

    #[test]
    fn log2_frac_tables_match_appendix_a() {
        let m = mode();
        for i in 0..24 {
            assert_eq!(
                m.log2_frac_table[i],
                i32::from(extracted::LOG2_FRAC_TABLE[i]),
                "LOG2_FRAC_TABLE[{i}]"
            );
            // Also the ceiling of 8·log2 (CELT_SPEC §6.2).
            assert_eq!(
                m.log2_frac_table[i],
                (8.0 * ((i + 1) as f64).log2()).ceil() as i32
            );
        }
        for j in 0..NB_EBANDS {
            assert_eq!(m.log_n[j], i32::from(extracted::LOGN400[j]), "logN[{j}]");
        }
    }

    #[test]
    fn pulse_cache_matches_appendix_a() {
        let m = mode();
        for r in 0..5 {
            for j in 0..NB_EBANDS {
                assert_eq!(
                    m.cache_index[r][j],
                    i32::from(extracted::CACHE_INDEX50[r][j]),
                    "cache_index[{r}][{j}]"
                );
            }
        }
        assert_eq!(m.cache_bits, extracted::CACHE_BITS50.to_vec());
    }

    #[test]
    fn caps_match_appendix_a() {
        let m = mode();
        for r in 0..8 {
            for j in 0..NB_EBANDS {
                assert_eq!(
                    m.cache_caps[r][j],
                    i32::from(extracted::CACHE_CAPS50[r][j]),
                    "caps[{r}][{j}]"
                );
            }
        }
    }

    #[test]
    fn small_tables_match_appendix_a() {
        let m = mode();
        for i in 0..30 {
            assert_eq!(
                m.ordery[i] as i32,
                extracted::ORDERY_TABLE[i],
                "ordery[{i}]"
            );
        }
        for x in 0..16 {
            assert_eq!(
                m.bit_interleave[x],
                u32::from(extracted::BIT_INTERLEAVE_TABLE[x])
            );
            assert_eq!(
                m.bit_deinterleave[x],
                u32::from(extracted::BIT_DEINTERLEAVE_TABLE[x])
            );
        }
        for k in 0..8 {
            assert_eq!(m.exp2_table8[k], i32::from(extracted::EXP2_TABLE8[k]));
        }
    }

    /// The float coefficient tables are their fixed-point forms scaled
    /// exactly (CELT_SPEC §2.3, §10.1).
    #[test]
    fn float_tables_are_scaled_fixed_point() {
        for i in 0..25 {
            assert_eq!(E_MEANS[i], f32::from(extracted::E_MEANS_Q4[i]) / 16.0);
        }
        for k in 0..4 {
            assert_eq!(
                PRED_COEF[k],
                f32::from(extracted::PRED_COEF_Q15[k]) / 32768.0
            );
            assert_eq!(
                BETA_COEF[k],
                f32::from(extracted::BETA_COEF_Q15[k]) / 32768.0
            );
        }
        assert_eq!(BETA_INTRA, f32::from(extracted::BETA_INTRA_Q15) / 32768.0);
    }

    /// The window table is the RFC 6716 §4.3.7 formula to within one float
    /// ulp (CELT_SPEC §10.4 notes five entries one ulp off its rounding).
    #[test]
    fn window_is_the_rfc_formula() {
        let mut off = 0;
        for (i, &w) in WINDOW120.iter().enumerate() {
            let s = (std::f64::consts::FRAC_PI_2 * (i as f64 + 0.5) / 120.0).sin();
            let f = (std::f64::consts::FRAC_PI_2 * s * s).sin() as f32;
            let ulps = (f.to_bits() as i64 - w.to_bits() as i64).abs();
            assert!(ulps <= 1, "window[{i}]: {w} vs {f}");
            off += usize::from(ulps != 0);
        }
        assert!(off <= 5);
    }

    /// log2_frac never undershoots the true log2 and is within 1/8 bit
    /// above it.
    #[test]
    fn log2_frac_is_conservative() {
        for v in (1u32..5000).chain([65535, 65537, 1 << 20 | 1, u32::MAX]) {
            let exact = 8.0 * f64::from(v).log2();
            let l = f64::from(log2_frac(v, 3));
            assert!(
                l >= exact - 1e-9 && l <= exact + 1.0 + 1e-9,
                "{v}: {l} vs {exact}"
            );
        }
    }
}
