//! The PVQ codebook of RFC 6716 §4.3.4.2: V(N, K), the number of integer
//! vectors of dimension N whose absolute values sum to K, and the bijection
//! between such vectors and indices `0 .. V(N, K)`.

/// V(N, K) by the recurrence V(N,K) = V(N-1,K) + V(N,K-1) + V(N-1,K-1),
/// V(N,0) = 1, V(0,K) = 0 for K != 0. Saturates at `u64::MAX`.
pub fn v(n: usize, k: usize) -> u64 {
    if let Some(t) = Shared::covers(n, k) {
        let x = t.get(n, k);
        if x != SATURATED {
            return u64::from(x);
        }
    }
    v_direct(n, k)
}

fn v_direct(n: usize, k: usize) -> u64 {
    let mut row = vec![0u64; k + 1];
    row[0] = 1; // n = 0
    for _ in 0..n {
        // new[k] = old[k] + new[k-1] + old[k-1]
        let mut prev_old = row[0];
        row[0] = 1;
        for j in 1..=k {
            let old = row[j];
            row[j] = old.saturating_add(row[j - 1]).saturating_add(prev_old);
            prev_old = old;
        }
    }
    row[k]
}

/// The largest band (N) and pulse count (K) CELT codes: band 20's 22 bins
/// at LM 3 (8 × 22), and `get_pulses(MAX_PSEUDO)`. Larger codebooks fall
/// back to a table of their own.
const TABLE_N: usize = 176;
const TABLE_K: usize = 240;

/// `u32::MAX` marks a V(N, K) of 2^32 or more. No codebook has exactly
/// 2^32 - 1 entries: V(N, K) is even for K > 0 (each codeword's sign
/// flipped is another) and 1 for K = 0.
const SATURATED: u32 = u32::MAX;

/// V(N, K) for every `N <= TABLE_N`, `K <= TABLE_K`, built once. Codebooks
/// the range coder can carry have fewer than 2^32 entries, and every value
/// [`decode`] and [`encode`] read for them is at most V(N, K), so 32 bits
/// hold them exactly.
struct Shared {
    v: Vec<u32>,
}

impl Shared {
    const W: usize = TABLE_K + 1;

    fn covers(n: usize, k: usize) -> Option<&'static Self> {
        static T: std::sync::OnceLock<Shared> = std::sync::OnceLock::new();
        (n <= TABLE_N && k <= TABLE_K).then(|| {
            T.get_or_init(|| {
                let t = table(TABLE_N, TABLE_K);
                Shared {
                    v: t.iter()
                        .map(|&x| u32::try_from(x).unwrap_or(SATURATED))
                        .collect(),
                }
            })
        })
    }

    #[inline(always)]
    fn get(&self, n: usize, k: usize) -> u32 {
        self.v[n * Self::W + k]
    }
}

/// The table `V[n][k]` for `n <= N`, `k <= K`, flattened by rows of `K + 1`.
fn table(n: usize, k: usize) -> Vec<u64> {
    let w = k + 1;
    let mut t = vec![0u64; (n + 1) * w];
    t[0] = 1;
    for i in 1..=n {
        t[i * w] = 1;
        for j in 1..=k {
            t[i * w + j] = t[(i - 1) * w + j]
                .saturating_add(t[i * w + j - 1])
                .saturating_add(t[(i - 1) * w + j - 1]);
        }
    }
    t
}

/// The shared table if it holds V(n, k) exactly.
fn shared_exact(n: usize, k: usize) -> Option<&'static Shared> {
    Shared::covers(n, k).filter(|t| t.get(n, k) != SATURATED)
}

/// The vector with index `i` among the V(N, K) codewords (§4.3.4.2), into
/// `y[..n]`.
pub fn decode(i: u32, n: usize, k: usize, y: &mut [i32]) {
    if let Some(t) = shared_exact(n, k) {
        return decode_with(|nn, kk| u64::from(t.get(nn, kk)), i, n, k, y);
    }
    let t = table(n, k);
    decode_with(|nn, kk| t[nn * (k + 1) + kk], i, n, k, y);
}

#[inline(always)]
fn decode_with(vv: impl Fn(usize, usize) -> u64, i: u32, n: usize, k: usize, y: &mut [i32]) {
    let mut i = u64::from(i);
    let mut k = k;
    for j in 0..n {
        let rem = n - j;
        let mut p = (vv(rem - 1, k) + vv(rem, k)) / 2;
        let neg = i >= p;
        if neg {
            i -= p;
        }
        let k0 = k;
        p -= vv(rem - 1, k);
        while p > i {
            k -= 1;
            p -= vv(rem - 1, k);
        }
        let mag = (k0 - k) as i32;
        y[j] = if neg { -mag } else { mag };
        i -= p;
    }
}

/// The index of `y` (whose absolute values sum to K) — the inverse of
/// [`decode`].
pub fn encode(y: &[i32], k: usize) -> u32 {
    if let Some(t) = shared_exact(y.len(), k) {
        return encode_with(|nn, kk| u64::from(t.get(nn, kk)), y, k);
    }
    let t = table(y.len(), k);
    encode_with(|nn, kk| t[nn * (k + 1) + kk], y, k)
}

#[inline(always)]
fn encode_with(vv: impl Fn(usize, usize) -> u64, y: &[i32], k: usize) -> u32 {
    let n = y.len();
    let mut i = 0u64;
    let mut k = k;
    for (j, &yj) in y.iter().enumerate() {
        let rem = n - j;
        if yj < 0 {
            i += (vv(rem - 1, k) + vv(rem, k)) / 2;
        }
        let m = yj.unsigned_abs() as usize;
        // Magnitudes are laid out largest first: the block of magnitude m
        // starts after those of magnitudes k .. m+1.
        for r in 0..k - m {
            i += vv(rem - 1, r);
        }
        k -= m;
    }
    i as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_values() {
        assert_eq!(v(1, 5), 2);
        assert_eq!(v(2, 1), 4);
        assert_eq!(v(4, 1), 8);
        assert_eq!(v(4, 2), 32);
        assert_eq!(v(3, 0), 1);
        assert_eq!(v(0, 3), 0);
        // V(N,K) = sum over the number of non-zero positions m of
        // 2^m C(N,m) C(K-1,m-1).
        fn c(n: u64, k: u64) -> u64 {
            if k > n {
                return 0;
            }
            (0..k).fold(1, |acc, i| acc * (n - i) / (i + 1))
        }
        for n in 1..9u64 {
            for k in 1..9u64 {
                let direct: u64 = (1..=n.min(k))
                    .map(|m| (1 << m) * c(n, m) * c(k - 1, m - 1))
                    .sum();
                assert_eq!(v(n as usize, k as usize), direct, "V({n},{k})");
            }
        }
    }

    /// Every codeword of small codebooks is enumerated exactly once:
    /// decoding each index gives a distinct vector of L1 norm K, and encoding
    /// it gives the index back.
    /// The shared table agrees with the recurrence wherever it is exact,
    /// and marks exactly the codebooks of 2^32 entries or more.
    #[test]
    fn shared_table_matches_recurrence() {
        let t = Shared::covers(TABLE_N, TABLE_K).unwrap();
        let full = table(TABLE_N, TABLE_K);
        for n in 0..=TABLE_N {
            for kk in 0..=TABLE_K {
                let x = full[n * (TABLE_K + 1) + kk];
                assert_eq!(t.get(n, kk) == SATURATED, x >= 1 << 32, "V({n},{kk})");
                assert_eq!(v(n, kk), x, "V({n},{kk})");
                if x < 1 << 32 {
                    assert_eq!(u64::from(t.get(n, kk)), x);
                }
            }
        }
        assert_eq!(v(300, 2), v_direct(300, 2));
    }

    #[test]
    fn enumeration_is_a_bijection() {
        for n in 1..=6 {
            for k in 1..=6 {
                let total = v(n, k) as u32;
                let mut seen = std::collections::HashSet::new();
                for i in 0..total {
                    let mut y = vec![0; n];
                    decode(i, n, k, &mut y);
                    assert_eq!(
                        y.iter().map(|x| x.unsigned_abs() as usize).sum::<usize>(),
                        k
                    );
                    assert!(seen.insert(y.clone()), "duplicate {y:?}");
                    assert_eq!(encode(&y, k), i, "N={n} K={k} y={y:?}");
                }
                assert_eq!(seen.len() as u32, total);
            }
        }
    }

    #[test]
    fn large_codebooks_round_trip() {
        let mut seed = 12345u32;
        for &(n, k) in &[
            (176usize, 3usize),
            (16, 30),
            (8, 128),
            (2, 128),
            (40, 12),
            (3, 1000),
        ] {
            for _ in 0..50 {
                let mut y = vec![0i32; n];
                for _ in 0..k {
                    seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    let pos = (seed >> 8) as usize % n;
                    let s = if y[pos] != 0 {
                        y[pos].signum()
                    } else if seed & 1 == 0 {
                        1
                    } else {
                        -1
                    };
                    y[pos] += s;
                }
                if v(n, k) > u64::from(u32::MAX) {
                    continue;
                }
                let i = encode(&y, k);
                let mut back = vec![0; n];
                decode(i, n, k, &mut back);
                assert_eq!(back, y);
            }
        }
    }
}
