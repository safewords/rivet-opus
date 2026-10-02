//! Derived constants of the 48 kHz CELT mode: the per-band log2 widths, the
//! pulse cache that maps bit budgets to pulse counts (§4.3.4.1), the band
//! caps (§4.3.3) and the low-overlap window (§4.3.7). Everything here is
//! computed from its mathematical definition once, on first use.

use std::sync::OnceLock;

use super::cwrs;
use super::tables::{CACHE_CAPS, EBANDS, NB_EBANDS, OVERLAP};

/// Fractional bits of the allocation (1/8 bit).
pub const BITRES: i32 = 3;
/// Largest pseudo-pulse index in the cache.
const MAX_PSEUDO: usize = 40;

/// A conservative log2 of `val` in `1/2^frac` units, rounded up to the
/// next representable step for values that are not powers of two.
pub fn log2_frac(val: u32, frac: i32) -> i32 {
    let mut l = crate::range::ilog(val);
    if val & val.wrapping_sub(1) != 0 {
        // Normalise to a Q15 mantissa in [32768, 65536), rounding up.
        let mut v: u64 = if l > 16 { (u64::from(val - 1) >> (l - 16)) + 1 } else { u64::from(val) << (16 - l) };
        l = (l - 1) << frac;
        let mut f = frac;
        loop {
            let b = (v >> 16) as i32;
            l += b << f;
            v = (v + b as u64) >> b;
            v = (v * v + 0x7FFF) >> 15;
            if f == 0 {
                break;
            }
            f -= 1;
        }
        l + i32::from(v > 0x8000)
    } else {
        (l - 1) << frac
    }
}

/// The pseudo-pulse index `i` as a pulse count: 0–7 exactly, then eight
/// steps per octave.
pub fn get_pulses(i: usize) -> usize {
    if i < 8 { i } else { (8 + (i & 7)) << ((i >> 3) - 1) }
}

/// Mode data computed once.
pub struct Mode {
    /// log2 of each band's width at LM 0, in 1/8 bits.
    pub log_n: [i32; NB_EBANDS],
    /// For `(LM + 1) * NB_EBANDS + band`, the start of that band's entry
    /// in `cache_bits` (entries shared between equal widths).
    cache_index: Vec<usize>,
    /// Entry `[0]` is the largest pseudo-pulse index; `[k]` the cost in
    /// 1/8 bits, minus one, of `get_pulses(k)` pulses.
    cache_bits: Vec<u8>,
    /// The low-overlap window's rising half, `OVERLAP` samples.
    pub window: Vec<f32>,
}

impl Mode {
    fn new() -> Self {
        let mut log_n = [0; NB_EBANDS];
        for (i, l) in log_n.iter_mut().enumerate() {
            *l = log2_frac((EBANDS[i + 1] - EBANDS[i]) as u32, BITRES);
        }
        let mut cache_index = vec![usize::MAX; 5 * NB_EBANDS];
        let mut cache_bits = Vec::new();
        let mut entries: Vec<(usize, usize)> = Vec::new(); // (N, start)
        for lm1 in 0..5 {
            for j in 0..NB_EBANDS {
                let n = ((EBANDS[j + 1] - EBANDS[j]) << lm1) >> 1;
                if n == 0 {
                    continue;
                }
                if let Some(&(_, start)) = entries.iter().find(|e| e.0 == n) {
                    cache_index[lm1 * NB_EBANDS + j] = start;
                    continue;
                }
                let fits = |k: usize| cwrs::v(n, k) <= u64::from(u32::MAX);
                let mut kmax = 0;
                while kmax < MAX_PSEUDO && fits(get_pulses(kmax + 1)) {
                    kmax += 1;
                }
                let start = cache_bits.len();
                cache_bits.push(kmax as u8);
                for k in 1..=kmax {
                    let vk = cwrs::v(n, get_pulses(k)) as u32;
                    cache_bits.push((log2_frac(vk, BITRES) - 1) as u8);
                }
                entries.push((n, start));
                cache_index[lm1 * NB_EBANDS + j] = start;
            }
        }
        let window = (0..OVERLAP)
            .map(|i| {
                let x = (std::f64::consts::FRAC_PI_2 * (i as f64 + 0.5) / OVERLAP as f64).sin();
                (std::f64::consts::FRAC_PI_2 * x * x).sin() as f32
            })
            .collect();
        Self { log_n, cache_index, cache_bits, window }
    }

    /// The cache row for `band` at `lm` (which may be -1 after a split).
    pub fn cache(&self, band: usize, lm: i32) -> &[u8] {
        let start = self.cache_index[((lm + 1) as usize) * NB_EBANDS + band];
        let len = usize::from(self.cache_bits[start]) + 1;
        &self.cache_bits[start..start + len]
    }

    /// §4.3.4.1: the pseudo-pulse index whose cost is nearest `bits`
    /// (1/8 bits), rounding down at the halfway point.
    pub fn bits2pulses(&self, band: usize, lm: i32, bits: i32) -> usize {
        let cache = self.cache(band, lm);
        let mut lo = 0usize;
        let mut hi = usize::from(cache[0]);
        let bits = bits - 1;
        for _ in 0..6 {
            let mid = (lo + hi + 1) >> 1;
            if i32::from(cache[mid]) >= bits {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let lo_bits = if lo == 0 { -1 } else { i32::from(cache[lo]) };
        if bits - lo_bits <= i32::from(cache[hi]) - bits { lo } else { hi }
    }

    /// The cost in 1/8 bits of pseudo-pulse index `q`.
    pub fn pulses2bits(&self, band: usize, lm: i32, q: usize) -> i32 {
        if q == 0 { 0 } else { i32::from(self.cache(band, lm)[q]) + 1 }
    }
}

/// The shared mode.
pub fn mode() -> &'static Mode {
    static MODE: OnceLock<Mode> = OnceLock::new();
    MODE.get_or_init(Mode::new)
}

/// §4.3.3: the per-band maximum allocation `cap[]` in 1/8 bits.
pub fn init_caps(lm: usize, c: usize) -> [i32; NB_EBANDS] {
    let mut cap = [0; NB_EBANDS];
    for (i, cp) in cap.iter_mut().enumerate() {
        let n = ((EBANDS[i + 1] - EBANDS[i]) << lm) as i32;
        *cp = (i32::from(CACHE_CAPS[NB_EBANDS * (2 * lm + c - 1) + i]) + 64) * c as i32 * n >> 2;
    }
    cap
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log2_frac_is_a_conservative_eighth_bit_log() {
        for v in 1..5000u32 {
            let l = log2_frac(v, 3);
            let exact = (v as f64).log2() * 8.0;
            assert!(f64::from(l) >= exact - 1e-9, "{v}: {l} < {exact}");
            assert!(f64::from(l) < exact + 1.0 + 1e-9, "{v}: {l} vs {exact}");
        }
        // The band widths' logs (in 1/8 bit).
        let m = mode();
        assert_eq!(m.log_n, [0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 16, 16, 16, 21, 21, 24, 29, 34, 36]);
    }

    #[test]
    fn cache_costs_rise_and_round_trip() {
        let m = mode();
        for lm in -1..=3 {
            for band in 0..NB_EBANDS {
                let n = ((EBANDS[band + 1] - EBANDS[band]) << (lm + 1)) >> 1;
                if n <= 1 {
                    // N = 1 has only a sign: every pulse count costs a bit.
                    continue;
                }
                let c = m.cache(band, lm);
                for k in 2..c.len() {
                    assert!(c[k] >= c[k - 1], "band {band} lm {lm}");
                }
                for q in 0..c.len() {
                    let bits = m.pulses2bits(band, lm, q);
                    let back = m.bits2pulses(band, lm, bits);
                    assert_eq!(m.pulses2bits(band, lm, back), bits, "band {band} lm {lm} q {q}");
                }
            }
        }
    }

    /// The window is power complementary (Princen-Bradley).
    #[test]
    fn window_is_power_complementary() {
        let w = &mode().window;
        for i in 0..OVERLAP {
            let s = w[i] * w[i] + w[OVERLAP - 1 - i] * w[OVERLAP - 1 - i];
            assert!((s - 1.0).abs() < 1e-6);
        }
    }
}
