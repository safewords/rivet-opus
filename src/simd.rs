//! Vector kernels, chosen at run time (AVX2 on x86-64, NEON on AArch64)
//! with a portable fallback.
//!
//! Every path computes the same thing in the same order: eight lanes, each
//! a plain multiply then add (no fused multiply-add), reduced by the same
//! tree. Their results are therefore bit-identical to one another, so the
//! codec's output does not depend on the CPU it runs on. The `force-scalar`
//! feature compiles the vector paths out (CI tests the fallback with it).

/// Lanes of the accumulation, whatever the instruction set.
const LANES: usize = 8;

/// The fixed reduction of the eight lane sums:
/// `((l0 + l4) + (l2 + l6)) + ((l1 + l5) + (l3 + l7))`, the order the
/// vector paths reduce in (halves, then quarters, then the pair).
#[inline(always)]
fn reduce(l: [f32; LANES]) -> f32 {
    let h = [l[0] + l[4], l[1] + l[5], l[2] + l[6], l[3] + l[7]];
    (h[0] + h[2]) + (h[1] + h[3])
}

/// `Σ a[i]·b[i]` over `a.len()` (`b` at least as long), in the lane order
/// described above.
#[inline]
pub(crate) fn dot(a: &[f32], b: &[f32]) -> f32 {
    let b = &b[..a.len()];
    #[cfg(all(target_arch = "x86_64", not(feature = "force-scalar")))]
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 is available on this CPU (checked just above).
        return unsafe { x86::dot_avx2(a, b) };
    }
    #[cfg(all(target_arch = "aarch64", not(feature = "force-scalar")))]
    {
        // SAFETY: NEON is part of the AArch64 baseline.
        return unsafe { arm::dot_neon(a, b) };
    }
    #[allow(unreachable_code)]
    dot_scalar(a, b)
}

/// The portable definition of [`dot`].
pub(crate) fn dot_scalar(a: &[f32], b: &[f32]) -> f32 {
    let b = &b[..a.len()];
    let mut acc = [0.0f32; LANES];
    let (ca, ra) = a.as_chunks::<LANES>();
    let (cb, rb) = b.as_chunks::<LANES>();
    for (x, y) in ca.iter().zip(cb) {
        for j in 0..LANES {
            acc[j] += x[j] * y[j];
        }
    }
    for (j, (x, y)) in ra.iter().zip(rb).enumerate() {
        acc[j] += x * y;
    }
    reduce(acc)
}

#[cfg(all(target_arch = "x86_64", not(feature = "force-scalar")))]
#[allow(unsafe_code)]
mod x86 {
    use super::{LANES, reduce};
    use std::arch::x86_64::*;

    /// [`super::dot`] with 256-bit vectors: one vector holds the eight
    /// lanes.
    ///
    /// # Safety
    /// The CPU must support AVX2, and `b.len() >= a.len()`.
    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn dot_avx2(a: &[f32], b: &[f32]) -> f32 {
        let n = a.len() / LANES * LANES;
        let mut acc = _mm256_setzero_ps();
        let mut i = 0;
        while i < n {
            // SAFETY: i + 8 <= n <= a.len() <= b.len().
            let (x, y) = unsafe { (_mm256_loadu_ps(a.as_ptr().add(i)), _mm256_loadu_ps(b.as_ptr().add(i))) };
            acc = _mm256_add_ps(acc, _mm256_mul_ps(x, y));
            i += LANES;
        }
        let mut l = [0.0f32; LANES];
        // SAFETY: `l` holds eight floats.
        unsafe { _mm256_storeu_ps(l.as_mut_ptr(), acc) };
        for j in 0..a.len() - n {
            l[j] += a[n + j] * b[n + j];
        }
        reduce(l)
    }
}

#[cfg(all(target_arch = "aarch64", not(feature = "force-scalar")))]
#[allow(unsafe_code)]
mod arm {
    use super::{LANES, reduce};
    use std::arch::aarch64::*;

    /// [`super::dot`] with two 128-bit vectors holding lanes 0-3 and 4-7.
    ///
    /// # Safety
    /// `b.len() >= a.len()`.
    #[target_feature(enable = "neon")]
    pub(super) unsafe fn dot_neon(a: &[f32], b: &[f32]) -> f32 {
        let n = a.len() / LANES * LANES;
        let (mut lo, mut hi) = (vdupq_n_f32(0.0), vdupq_n_f32(0.0));
        let mut i = 0;
        while i < n {
            // SAFETY: i + 8 <= n <= a.len() <= b.len().
            unsafe {
                let (pa, pb) = (a.as_ptr().add(i), b.as_ptr().add(i));
                lo = vaddq_f32(lo, vmulq_f32(vld1q_f32(pa), vld1q_f32(pb)));
                hi = vaddq_f32(hi, vmulq_f32(vld1q_f32(pa.add(4)), vld1q_f32(pb.add(4))));
            }
            i += LANES;
        }
        let mut l = [0.0f32; LANES];
        // SAFETY: `l` holds eight floats.
        unsafe {
            vst1q_f32(l.as_mut_ptr(), lo);
            vst1q_f32(l.as_mut_ptr().add(4), hi);
        }
        for j in 0..a.len() - n {
            l[j] += a[n + j] * b[n + j];
        }
        reduce(l)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn noise(len: usize, seed: u32) -> Vec<f32> {
        let mut s = seed;
        (0..len)
            .map(|_| {
                s = s.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                (s >> 8) as f32 / (1u32 << 24) as f32 - 0.5
            })
            .collect()
    }

    /// The dispatched kernel is bit-identical to the portable one at every
    /// length (tails included), on random data and on extremes.
    #[test]
    fn dot_matches_scalar_bit_for_bit() {
        for len in 0..200 {
            let a = noise(len, len as u32 + 1);
            let b = noise(len + 3, len as u32 + 1000);
            assert_eq!(dot(&a, &b).to_bits(), dot_scalar(&a, &b).to_bits(), "len {len}");
        }
        let big: Vec<f32> = (0..37).map(|i| if i % 3 == 0 { 3.0e38 } else { -1.0e-38 }).collect();
        assert_eq!(dot(&big, &big).to_bits(), dot_scalar(&big, &big).to_bits());
        let mixed: Vec<f32> = (0..64).map(|i| [1.0e20, -1.0e20, 1.0, f32::MIN_POSITIVE][i % 4]).collect();
        assert_eq!(dot(&mixed[1..], &mixed).to_bits(), dot_scalar(&mixed[1..], &mixed).to_bits());
    }
}
