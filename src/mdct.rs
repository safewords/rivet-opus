//! A mixed-radix FFT and the MDCT built on it.
//!
//! The MDCT here is the textbook transform, unscaled:
//!
//! ```text
//!   X[k] = sum_{n<2N} x[n] cos(pi/N (n + 1/2 + N/2)(k + 1/2)),   k < N
//!   y[n] = sum_{k<N}  X[k] cos(pi/N (n + 1/2 + N/2)(k + 1/2)),   n < 2N
//! ```
//!
//! computed as a DCT-IV of a folded input through an N/2-point complex FFT.
//! Windowing, overlap-add and scaling belong to the callers.

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Complex {
    pub re: f32,
    pub im: f32,
}

impl Complex {
    #[inline]
    fn mul(self, o: Complex) -> Complex {
        Complex { re: self.re * o.re - self.im * o.im, im: self.re * o.im + self.im * o.re }
    }
    #[inline]
    fn add(self, o: Complex) -> Complex {
        Complex { re: self.re + o.re, im: self.im + o.im }
    }
}

/// A forward complex FFT of a fixed size whose prime factors are 2, 3 and 5.
pub(crate) struct Fft {
    n: usize,
    factors: Vec<usize>,
    /// `e^{-2 pi i j / n}`.
    twiddles: Vec<Complex>,
}

impl Fft {
    pub fn new(n: usize) -> Self {
        let mut factors = Vec::new();
        let mut m = n;
        for p in [4, 2, 3, 5] {
            while m % p == 0 {
                factors.push(p);
                m /= p;
            }
        }
        assert_eq!(m, 1, "FFT size {n} has a prime factor above 5");
        let twiddles = (0..n)
            .map(|j| {
                let a = -2.0 * std::f64::consts::PI * j as f64 / n as f64;
                Complex { re: a.cos() as f32, im: a.sin() as f32 }
            })
            .collect();
        Self { n, factors, twiddles }
    }

    /// `out = FFT(inp)`.
    pub fn process(&self, inp: &[Complex], out: &mut [Complex]) {
        self.work(out, inp, 0, 1, 0);
    }

    fn work(&self, out: &mut [Complex], inp: &[Complex], offset: usize, stride: usize, level: usize) {
        let n = out.len();
        if n == 1 {
            out[0] = inp[offset];
            return;
        }
        let p = self.factors[level];
        let m = n / p;
        for q in 0..p {
            self.work(&mut out[q * m..(q + 1) * m], inp, offset + q * stride, stride * p, level + 1);
        }
        let tw_step = self.n / n;
        let mut a = [Complex::default(); 5];
        for k in 0..m {
            for q in 0..p {
                let t = self.twiddles[(q * k * tw_step) % self.n];
                a[q] = out[q * m + k].mul(t);
            }
            for q2 in 0..p {
                let mut s = Complex::default();
                for (q, aq) in a.iter().enumerate().take(p) {
                    let idx = ((q * q2) % p) * (self.n / p);
                    s = s.add(aq.mul(self.twiddles[idx]));
                }
                out[q2 * m + k] = s;
            }
        }
    }
}

/// An MDCT with `n` coefficients (`2n` time samples).
pub(crate) struct Mdct {
    n: usize,
    fft: Fft,
    pre: Vec<Complex>,
    post: Vec<Complex>,
}

impl Mdct {
    pub fn new(n: usize) -> Self {
        assert!(n % 4 == 0);
        let pre = (0..n / 2)
            .map(|i| {
                let a = -std::f64::consts::PI * i as f64 / n as f64;
                Complex { re: a.cos() as f32, im: a.sin() as f32 }
            })
            .collect();
        let post = (0..n / 2)
            .map(|k| {
                let a = -std::f64::consts::PI * (k as f64 + 0.25) / n as f64;
                Complex { re: a.cos() as f32, im: a.sin() as f32 }
            })
            .collect();
        Self { n, fft: Fft::new(n / 2), pre, post }
    }

    /// The number of coefficients.
    pub fn len(&self) -> usize {
        self.n
    }

    /// `out[k] = sum_n v[n] cos(pi/N (n + 1/2)(k + 1/2))`.
    pub fn dct4(&self, v: &[f32], out: &mut [f32]) {
        let n = self.n;
        let h = n / 2;
        let z: Vec<Complex> =
            (0..h).map(|i| Complex { re: v[2 * i], im: v[n - 1 - 2 * i] }.mul(self.pre[i])).collect();
        let mut zf = vec![Complex::default(); h];
        self.fft.process(&z, &mut zf);
        for k in 0..h {
            let y = zf[k].mul(self.post[k]);
            out[2 * k] = y.re;
            out[n - 1 - 2 * k] = -y.im;
        }
    }

    /// The forward MDCT of `x[..2n]` into `out[..n]`.
    pub fn forward(&self, x: &[f32], out: &mut [f32]) {
        let n = self.n;
        let h = n / 2;
        let mut v = vec![0.0f32; n];
        for i in 0..h {
            v[i] = -x[3 * h - 1 - i] - x[3 * h + i];
            v[h + i] = x[i] - x[n - 1 - i];
        }
        self.dct4(&v, out);
    }

    /// The inverse MDCT of `coefs[..n]` into `y[..2n]`.
    pub fn inverse(&self, coefs: &[f32], y: &mut [f32]) {
        let n = self.n;
        let h = n / 2;
        let mut u = vec![0.0f32; n];
        self.dct4(coefs, &mut u);
        for j in 0..h {
            y[j] = u[j + h];
            y[h + j] = -u[n - 1 - j];
            y[n + j] = -u[h - 1 - j];
            y[3 * h + j] = -u[j];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn direct_mdct(x: &[f32], n: usize) -> Vec<f64> {
        (0..n)
            .map(|k| {
                (0..2 * n)
                    .map(|i| {
                        f64::from(x[i])
                            * (std::f64::consts::PI / n as f64 * (i as f64 + 0.5 + n as f64 / 2.0) * (k as f64 + 0.5))
                                .cos()
                    })
                    .sum()
            })
            .collect()
    }

    fn noise(len: usize, seed: u32) -> Vec<f32> {
        let mut s = seed;
        (0..len)
            .map(|_| {
                s = s.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                (s >> 8) as f32 / (1u32 << 24) as f32 - 0.5
            })
            .collect()
    }

    #[test]
    fn fft_matches_dft() {
        for n in [1usize, 2, 3, 4, 5, 6, 8, 15, 60, 120, 240, 480] {
            let re = noise(n, 7);
            let im = noise(n, 9);
            let x: Vec<Complex> = (0..n).map(|i| Complex { re: re[i], im: im[i] }).collect();
            let mut out = vec![Complex::default(); n];
            Fft::new(n).process(&x, &mut out);
            for k in 0..n {
                let (mut sr, mut si) = (0.0f64, 0.0f64);
                for (i, v) in x.iter().enumerate() {
                    let a = -2.0 * std::f64::consts::PI * (i * k) as f64 / n as f64;
                    sr += f64::from(v.re) * a.cos() - f64::from(v.im) * a.sin();
                    si += f64::from(v.re) * a.sin() + f64::from(v.im) * a.cos();
                }
                assert!((sr - f64::from(out[k].re)).abs() < 1e-3 && (si - f64::from(out[k].im)).abs() < 1e-3, "n {n} k {k}");
            }
        }
    }

    /// The fast MDCT agrees with the defining sum (a float reference).
    #[test]
    fn mdct_matches_definition() {
        for n in [8usize, 60, 120, 240, 480, 960] {
            let x = noise(2 * n, n as u32);
            let mut fast = vec![0.0; n];
            let m = Mdct::new(n);
            m.forward(&x, &mut fast);
            let slow = direct_mdct(&x, n);
            let scale = (n as f64).sqrt();
            for k in 0..n {
                assert!((slow[k] - f64::from(fast[k])).abs() < 2e-4 * scale, "n {n} k {k}: {} vs {}", slow[k], fast[k]);
            }
            // Inverse against its definition.
            let coefs = noise(n, 3 + n as u32);
            let mut y = vec![0.0; 2 * n];
            m.inverse(&coefs, &mut y);
            for (i, &yi) in y.iter().enumerate() {
                let d: f64 = (0..n)
                    .map(|k| {
                        f64::from(coefs[k])
                            * (std::f64::consts::PI / n as f64 * (i as f64 + 0.5 + n as f64 / 2.0) * (k as f64 + 0.5))
                                .cos()
                    })
                    .sum();
                assert!((d - f64::from(yi)).abs() < 2e-4 * scale, "inverse n {n} i {i}");
            }
        }
    }

    /// Windowed MDCT, inverse and overlap-add reconstruct the input (TDAC),
    /// scaled by N/2, with a sine window.
    #[test]
    fn tdac_reconstructs() {
        let n = 120;
        let m = Mdct::new(n);
        let w: Vec<f32> =
            (0..2 * n).map(|i| (std::f32::consts::PI * (i as f32 + 0.5) / (2 * n) as f32).sin()).collect();
        let x = noise(6 * n, 5);
        let mut acc = vec![0.0f32; 6 * n];
        for b in 0..5 {
            let seg: Vec<f32> = (0..2 * n).map(|i| x[b * n + i] * w[i]).collect();
            let mut c = vec![0.0; n];
            m.forward(&seg, &mut c);
            let mut y = vec![0.0; 2 * n];
            m.inverse(&c, &mut y);
            for i in 0..2 * n {
                acc[b * n + i] += y[i] * w[i] * 2.0 / n as f32;
            }
        }
        for i in n..5 * n {
            assert!((acc[i] - x[i]).abs() < 1e-4, "{i}: {} vs {}", acc[i], x[i]);
        }
    }
}
