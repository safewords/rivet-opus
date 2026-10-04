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
        Complex {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }
    #[inline]
    fn add(self, o: Complex) -> Complex {
        Complex {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }
}

/// One radix-`p` pass of the Stockham autosort FFT (decimation in
/// frequency): with the current sub-transform length `n = p·m` and stride
/// `s`, `y[k + s(p·q + t)] = w_n^{q·t} · Σ_r x[k + s(q + r·m)] · w_p^{r·t}`.
struct Pass {
    p: usize,
    m: usize,
    s: usize,
    /// `w_n^{q·t}` for `q < m`, `1 <= t < p`, by `q` then `t`.
    tw: Vec<Complex>,
}

/// A forward complex FFT of a fixed size whose prime factors are 2, 3 and 5.
pub(crate) struct Fft {
    n: usize,
    passes: Vec<Pass>,
}

fn unit(a: f64) -> Complex {
    Complex {
        re: a.cos() as f32,
        im: a.sin() as f32,
    }
}

/// `-i·z`.
#[inline(always)]
fn mul_neg_i(z: Complex) -> Complex {
    Complex {
        re: z.im,
        im: -z.re,
    }
}

#[inline(always)]
fn sub(a: Complex, b: Complex) -> Complex {
    Complex {
        re: a.re - b.re,
        im: a.im - b.im,
    }
}

#[inline(always)]
fn scale(a: Complex, f: f32) -> Complex {
    Complex {
        re: a.re * f,
        im: a.im * f,
    }
}

/// The largest FFT the codec runs (the 960-coefficient MDCT's), whose
/// scratch lives on the stack.
const MAX_STACK_FFT: usize = MAX_MDCT / 2;

impl Fft {
    pub fn new(n: usize) -> Self {
        let mut factors = Vec::new();
        let mut m = n;
        for p in [4, 2, 3, 5] {
            while m.is_multiple_of(p) {
                factors.push(p);
                m /= p;
            }
        }
        assert_eq!(m, 1, "FFT size {n} has a prime factor above 5");
        let mut passes = Vec::with_capacity(factors.len());
        let (mut len, mut s) = (n, 1);
        for &p in &factors {
            let m = len / p;
            let mut tw = Vec::with_capacity(m * (p - 1));
            for q in 0..m {
                for t in 1..p {
                    tw.push(unit(
                        -2.0 * std::f64::consts::PI * (q * t) as f64 / len as f64,
                    ));
                }
            }
            passes.push(Pass { p, m, s, tw });
            len = m;
            s *= p;
        }
        Self { n, passes }
    }

    /// `out = FFT(inp)`.
    pub fn process(&self, inp: &[Complex], out: &mut [Complex]) {
        let n = self.n;
        let (inp, out) = (&inp[..n], &mut out[..n]);
        if self.passes.is_empty() {
            out[0] = inp[0];
            return;
        }
        let mut stack = [Complex::default(); MAX_STACK_FFT];
        let mut heap = Vec::new();
        let tmp: &mut [Complex] = if n <= MAX_STACK_FFT {
            &mut stack[..n]
        } else {
            heap.resize(n, Complex::default());
            &mut heap
        };
        // The passes ping-pong between `out` and `tmp`, starting from `inp`,
        // so that the last one writes `out`.
        let mut to_out = self.passes.len() % 2 == 1;
        for (i, pass) in self.passes.iter().enumerate() {
            match (i == 0, to_out) {
                (true, true) => pass.run(inp, out),
                (true, false) => pass.run(inp, tmp),
                (false, true) => pass.run(tmp, out),
                (false, false) => pass.run(out, tmp),
            }
            to_out = !to_out;
        }
    }
}

impl Pass {
    fn run(&self, x: &[Complex], y: &mut [Complex]) {
        let (m, s) = (self.m, self.s);
        match self.p {
            2 => {
                for q in 0..m {
                    let w1 = self.tw[q];
                    for k in 0..s {
                        let a = x[k + s * q];
                        let b = x[k + s * (q + m)];
                        y[k + s * 2 * q] = a.add(b);
                        y[k + s * (2 * q + 1)] = sub(a, b).mul(w1);
                    }
                }
            }
            3 => {
                // w_3 = -1/2 - i·sqrt(3)/2.
                let h = 0.75f64.sqrt() as f32;
                for q in 0..m {
                    let (w1, w2) = (self.tw[2 * q], self.tw[2 * q + 1]);
                    for k in 0..s {
                        let a0 = x[k + s * q];
                        let a1 = x[k + s * (q + m)];
                        let a2 = x[k + s * (q + 2 * m)];
                        let t1 = a1.add(a2);
                        let t2 = sub(a0, scale(t1, 0.5));
                        // (a1 - a2)·(-i·sqrt(3)/2)
                        let t3 = mul_neg_i(scale(sub(a1, a2), h));
                        y[k + s * 3 * q] = a0.add(t1);
                        y[k + s * (3 * q + 1)] = t2.add(t3).mul(w1);
                        y[k + s * (3 * q + 2)] = sub(t2, t3).mul(w2);
                    }
                }
            }
            4 => {
                for q in 0..m {
                    let (w1, w2, w3) = (self.tw[3 * q], self.tw[3 * q + 1], self.tw[3 * q + 2]);
                    for k in 0..s {
                        let a0 = x[k + s * q];
                        let a1 = x[k + s * (q + m)];
                        let a2 = x[k + s * (q + 2 * m)];
                        let a3 = x[k + s * (q + 3 * m)];
                        let b0 = a0.add(a2);
                        let b1 = sub(a0, a2);
                        let b2 = a1.add(a3);
                        let b3 = mul_neg_i(sub(a1, a3));
                        y[k + s * 4 * q] = b0.add(b2);
                        y[k + s * (4 * q + 1)] = b1.add(b3).mul(w1);
                        y[k + s * (4 * q + 2)] = sub(b0, b2).mul(w2);
                        y[k + s * (4 * q + 3)] = sub(b1, b3).mul(w3);
                    }
                }
            }
            _ => {
                // Radix 5, with w_5^t = cos(2πt/5) - i·sin(2πt/5).
                let c1 = (0.4 * std::f64::consts::PI).cos() as f32;
                let c2 = (0.8 * std::f64::consts::PI).cos() as f32;
                let s1 = (0.4 * std::f64::consts::PI).sin() as f32;
                let s2 = (0.8 * std::f64::consts::PI).sin() as f32;
                for q in 0..m {
                    let tw = &self.tw[4 * q..4 * q + 4];
                    for k in 0..s {
                        let a0 = x[k + s * q];
                        let a1 = x[k + s * (q + m)];
                        let a2 = x[k + s * (q + 2 * m)];
                        let a3 = x[k + s * (q + 3 * m)];
                        let a4 = x[k + s * (q + 4 * m)];
                        let (p1, m1) = (a1.add(a4), sub(a1, a4));
                        let (p2, m2) = (a2.add(a3), sub(a2, a3));
                        let r1 = a0.add(scale(p1, c1)).add(scale(p2, c2));
                        let r2 = a0.add(scale(p1, c2)).add(scale(p2, c1));
                        let i1 = mul_neg_i(scale(m1, s1).add(scale(m2, s2)));
                        let i2 = mul_neg_i(sub(scale(m1, s2), scale(m2, s1)));
                        y[k + s * 5 * q] = a0.add(p1).add(p2);
                        y[k + s * (5 * q + 1)] = r1.add(i1).mul(tw[0]);
                        y[k + s * (5 * q + 2)] = r2.add(i2).mul(tw[1]);
                        y[k + s * (5 * q + 3)] = sub(r2, i2).mul(tw[2]);
                        y[k + s * (5 * q + 4)] = sub(r1, i1).mul(tw[3]);
                    }
                }
            }
        }
    }
}

/// The largest MDCT (coefficients): CELT's 20 ms frame. Working buffers of
/// this size live on the stack.
const MAX_MDCT: usize = 960;

/// An MDCT with `n` coefficients (`2n` time samples).
pub(crate) struct Mdct {
    n: usize,
    fft: Fft,
    pre: Vec<Complex>,
    post: Vec<Complex>,
}

impl Mdct {
    pub fn new(n: usize) -> Self {
        assert!(n.is_multiple_of(4) && n <= MAX_MDCT, "MDCT size {n}");
        let pre = (0..n / 2)
            .map(|i| {
                let a = -std::f64::consts::PI * i as f64 / n as f64;
                Complex {
                    re: a.cos() as f32,
                    im: a.sin() as f32,
                }
            })
            .collect();
        let post = (0..n / 2)
            .map(|k| {
                let a = -std::f64::consts::PI * (k as f64 + 0.25) / n as f64;
                Complex {
                    re: a.cos() as f32,
                    im: a.sin() as f32,
                }
            })
            .collect();
        Self {
            n,
            fft: Fft::new(n / 2),
            pre,
            post,
        }
    }

    /// `out[k] = sum_n v[n] cos(pi/N (n + 1/2)(k + 1/2))`.
    pub fn dct4(&self, v: &[f32], out: &mut [f32]) {
        let n = self.n;
        let h = n / 2;
        let mut z = [Complex::default(); MAX_MDCT / 2];
        let mut zf = [Complex::default(); MAX_MDCT / 2];
        let (z, zf) = (&mut z[..h], &mut zf[..h]);
        for (i, zi) in z.iter_mut().enumerate() {
            *zi = Complex {
                re: v[2 * i],
                im: v[n - 1 - 2 * i],
            }
            .mul(self.pre[i]);
        }
        self.fft.process(z, zf);
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
        let mut v = [0.0f32; MAX_MDCT];
        let v = &mut v[..n];
        for i in 0..h {
            v[i] = -x[3 * h - 1 - i] - x[3 * h + i];
            v[h + i] = x[i] - x[n - 1 - i];
        }
        self.dct4(v, out);
    }

    /// The inverse MDCT of `coefs[..n]` into `y[..2n]`.
    pub fn inverse(&self, coefs: &[f32], y: &mut [f32]) {
        let n = self.n;
        let h = n / 2;
        let mut u = [0.0f32; MAX_MDCT];
        let u = &mut u[..n];
        self.dct4(coefs, u);
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
                            * (std::f64::consts::PI / n as f64
                                * (i as f64 + 0.5 + n as f64 / 2.0)
                                * (k as f64 + 0.5))
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
            let x: Vec<Complex> = (0..n)
                .map(|i| Complex {
                    re: re[i],
                    im: im[i],
                })
                .collect();
            let mut out = vec![Complex::default(); n];
            Fft::new(n).process(&x, &mut out);
            for k in 0..n {
                let (mut sr, mut si) = (0.0f64, 0.0f64);
                for (i, v) in x.iter().enumerate() {
                    let a = -2.0 * std::f64::consts::PI * (i * k) as f64 / n as f64;
                    sr += f64::from(v.re) * a.cos() - f64::from(v.im) * a.sin();
                    si += f64::from(v.re) * a.sin() + f64::from(v.im) * a.cos();
                }
                assert!(
                    (sr - f64::from(out[k].re)).abs() < 1e-3
                        && (si - f64::from(out[k].im)).abs() < 1e-3,
                    "n {n} k {k}"
                );
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
                assert!(
                    (slow[k] - f64::from(fast[k])).abs() < 2e-4 * scale,
                    "n {n} k {k}: {} vs {}",
                    slow[k],
                    fast[k]
                );
            }
            // Inverse against its definition.
            let coefs = noise(n, 3 + n as u32);
            let mut y = vec![0.0; 2 * n];
            m.inverse(&coefs, &mut y);
            for (i, &yi) in y.iter().enumerate() {
                let d: f64 = (0..n)
                    .map(|k| {
                        f64::from(coefs[k])
                            * (std::f64::consts::PI / n as f64
                                * (i as f64 + 0.5 + n as f64 / 2.0)
                                * (k as f64 + 0.5))
                                .cos()
                    })
                    .sum();
                assert!(
                    (d - f64::from(yi)).abs() < 2e-4 * scale,
                    "inverse n {n} i {i}"
                );
            }
        }
    }

    /// Windowed MDCT, inverse and overlap-add reconstruct the input (TDAC),
    /// scaled by N/2, with a sine window.
    #[test]
    fn tdac_reconstructs() {
        let n = 120;
        let m = Mdct::new(n);
        let w: Vec<f32> = (0..2 * n)
            .map(|i| (std::f32::consts::PI * (i as f32 + 0.5) / (2 * n) as f32).sin())
            .collect();
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
