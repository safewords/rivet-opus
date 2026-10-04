//! Sample-rate conversion between the SILK rates (8, 12, 16 kHz) and the
//! Opus API rates (8, 12, 16, 24, 48 kHz).
//!
//! RFC 6716 §4.2.9 makes the resampler non-normative ("a decoder can use any
//! method it wants") but makes its delay normative: the SILK output must lag
//! by the Table 54 allocation (0.538, 0.692, 0.706 ms for NB, MB, WB), which
//! the encoder compensates on its MDCT path. This delay applies at every
//! output rate, including the SILK rate itself.
//!
//! # Design
//!
//! A rational converter `in → out` runs, conceptually, at the common rate
//! `F = lcm(in, out)`: zero-stuff by `L = F / in`, filter with a prototype
//! low-pass `h[0..N]`, keep every `M = F / out`-th sample. It is executed
//! as `L` polyphase branches so only the kept outputs are computed.
//!
//! The prototype is a **delay-constrained weighted least-squares FIR**,
//! designed here at construction (nothing is tabulated). With
//! `lo = min(in, out)`, `fc = lo / 2` and `Δ = 0.05 · lo` (400 Hz at 8 kHz),
//! the target response is `D(f) = L · A(f) · exp(−j2πfτ/F)`: a pure delay
//! of `τ` samples of `F` (the Table 54 delay, for the decoder) with
//! `A(f) = 1` up to `fc − Δ`, a raised-cosine roll-off
//! `A(f) = ½(1 + cos(π(f − fc + Δ) / 2Δ))` across the transition and
//! `A(f) = 0` from `fc + Δ`. `h` minimises the weighted error
//! `∫_0^{F/2} W(f) |H(f) − D(f)|² df` with `W = 1` in the passband, `0.1`
//! across the transition, `100` in the near stopband `[fc + Δ, 0.75 · lo]`
//! and `10⁴` beyond. The near stopband holds the images (or aliases) of the
//! top of the band, where speech has little energy; the far stopband those
//! of its strong low frequencies. `N = ⌈5 F / (2Δ)⌉`, 50 taps per branch.
//!
//! The objective is a quadratic form in `h`, so `h` solves the normal
//! equations `R h = p` with `R[m][n] = ∫ W(f) cos(2πf(m−n)/F) df` (closed
//! form; a symmetric Toeplitz matrix) and
//! `p[n] = ∫ W(f) L A(f) cos(2πf(n−τ)/F) df` (closed form in the passband,
//! Simpson's rule across the transition), by Cholesky decomposition. As
//! `τ` is far shorter than `N / 2` the result is not linear-phase: it is
//! the least-squares approximation, within the band, of a pure delay of
//! exactly the allowed amount, as §4.2.9 anticipates ("may add a delay that
//! is not an exact integer, or is not linear-phase"). Each polyphase branch
//! is finally scaled to unit DC gain so a constant input leaves no image
//! tones at multiples of `in`.
//!
//! The weights and the two band fractions are the only free choices. They
//! trade passband flatness against image rejection under the short NB
//! delay budget (4.3 samples at 8 kHz), and were chosen on a small grid for
//! passband flatness within ±0.8 dB with the conformance metric of
//! RFC 6716 §6 met at every output rate. The response of every converter
//! the codec builds is measured by the unit tests below and tabulated in
//! `docs/VALIDATION.md`.

use crate::simd;
use std::sync::{Arc, Mutex, OnceLock};

/// Transition half-width as a fraction of the lower rate.
const TRANSITION: f64 = 0.05;
/// Weight of the transition band (raised-cosine target).
const TRANSITION_WEIGHT: f64 = 0.1;
/// Weight of the near stopband, `[fc + Δ, FAR_STOP · lo]`.
const NEAR_STOP_WEIGHT: f64 = 100.0;
/// Start of the far stopband, as a fraction of the lower rate.
const FAR_STOP: f64 = 0.75;
/// Weight of the far stopband, `[FAR_STOP · lo, F/2]`.
const STOP_WEIGHT: f64 = 10_000.0;
/// Filter length in units of `F / (2Δ)`.
const LENGTH_FACTOR: f64 = 5.0;

/// A streaming resampler for one channel.
pub(crate) struct Resampler {
    in_rate: usize,
    /// Output samples per period of the rate ratio, and input samples.
    phases: usize,
    in_step: usize,
    /// Kernel taps per phase, starting at input offset `first[p]` relative
    /// to the period's first input sample (shared by every resampler of the
    /// same conversion).
    kernel: Arc<Kernel>,
    /// Input history, oldest first.
    hist: Vec<f32>,
}

/// A converter's polyphase kernel: taps per phase and the input offset
/// each phase starts at.
struct Kernel {
    taps: Vec<Vec<f32>>,
    first: Vec<isize>,
}

/// The kernel of `in_rate → out_rate` with a `delay_ms` delay, designed on
/// first use and then shared (the design solves a 150-tap least-squares
/// problem: a decoder or encoder would otherwise pay for it on creation and
/// at every SILK rate change).
fn kernel(in_rate: usize, out_rate: usize, delay_ms: f64) -> Arc<Kernel> {
    type Cache = Mutex<Vec<((usize, usize, u64), Arc<Kernel>)>>;
    static CACHE: OnceLock<Cache> = OnceLock::new();
    let key = (in_rate, out_rate, delay_ms.to_bits());
    let cache = CACHE.get_or_init(Default::default);
    if let Some((_, k)) = cache
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .iter()
        .find(|(k, _)| *k == key)
    {
        return k.clone();
    }
    let k = Arc::new(design_kernel(in_rate, out_rate, delay_ms));
    let mut c = cache.lock().unwrap_or_else(|e| e.into_inner());
    if let Some((_, k)) = c.iter().find(|(k, _)| *k == key) {
        return k.clone();
    }
    c.push((key, k.clone()));
    k
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// `∫_a^b cos(ωk) dω`.
fn int_cos(a: f64, b: f64, k: f64) -> f64 {
    if k.abs() < 1e-12 {
        b - a
    } else {
        ((b * k).sin() - (a * k).sin()) / k
    }
}

/// Solves `R x = p` for symmetric positive definite `R` (row-major, `n×n`)
/// by Cholesky decomposition.
fn solve_spd(mut r: Vec<f64>, mut p: Vec<f64>, n: usize) -> Vec<f64> {
    for j in 0..n {
        let mut d = r[j * n + j];
        for k in 0..j {
            d -= r[j * n + k] * r[j * n + k];
        }
        let d = d.max(1e-300).sqrt();
        r[j * n + j] = d;
        for i in j + 1..n {
            let mut s = r[i * n + j];
            for k in 0..j {
                s -= r[i * n + k] * r[j * n + k];
            }
            r[i * n + j] = s / d;
        }
    }
    for i in 0..n {
        let mut s = p[i];
        for k in 0..i {
            s -= r[i * n + k] * p[k];
        }
        p[i] = s / r[i * n + i];
    }
    for i in (0..n).rev() {
        let mut s = p[i];
        for k in i + 1..n {
            s -= r[k * n + i] * p[k];
        }
        p[i] = s / r[i * n + i];
    }
    p
}

/// The delay-constrained least-squares prototype (see the module
/// documentation): `n` taps, band edges `wp`, `ws` and far-stopband start
/// `wf` in radians per sample, passband gain `gain` and delay `tau`
/// samples.
pub(crate) fn design(n: usize, wp: f64, ws: f64, wf: f64, tau: f64, gain: f64) -> Vec<f64> {
    use std::f64::consts::PI;
    let (near_w, stop_w, tw) = (NEAR_STOP_WEIGHT, STOP_WEIGHT, TRANSITION_WEIGHT);
    let we = ws.min(PI);
    let wf = wf.clamp(we, PI);
    let row: Vec<f64> = (0..n)
        .map(|k| {
            let k = k as f64;
            int_cos(0.0, wp, k)
                + tw * int_cos(wp, we, k)
                + near_w * int_cos(we, wf, k)
                + stop_w * int_cos(wf, PI, k)
        })
        .collect();
    let mut r = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            r[i * n + j] = row[i.abs_diff(j)];
        }
    }
    // A negligible ridge keeps the factorisation well defined.
    let ridge = 1e-10 * row[0];
    for i in 0..n {
        r[i * n + i] += ridge;
    }
    // The target: a pure delay in the passband, tapering as a raised cosine
    // to zero across the transition band (Simpson's rule there).
    let steps = 512;
    let h = (we - wp) / steps as f64;
    let p: Vec<f64> = (0..n)
        .map(|k| {
            let k = k as f64 - tau;
            let mut t = 0.0;
            for i in 0..=steps {
                let w = wp + i as f64 * h;
                let a = 0.5 + 0.5 * (PI * (w - wp) / (ws - wp)).cos();
                let c = if i == 0 || i == steps {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                t += c * a * (w * k).cos();
            }
            gain * (int_cos(0.0, wp, k) + tw * t * h / 3.0)
        })
        .collect();
    solve_spd(r, p, n)
}

/// Designs the kernel of a converter (see the module documentation).
fn design_kernel(in_rate: usize, out_rate: usize, delay_ms: f64) -> Kernel {
    let g = gcd(in_rate, out_rate);
    if in_rate == out_rate {
        // No filtering is needed, only the delay, rounded to whole
        // samples as §4.2.9 allows ("it may not be possible to achieve
        // exactly these delays while using a whole number of input or
        // output samples"): a transparent delay line.
        let d = (delay_ms * in_rate as f64 / 1000.0).round() as usize;
        return Kernel {
            taps: vec![vec![1.0]],
            first: vec![-(d as isize)],
        };
    }
    let up = out_rate / g; // L
    let down = in_rate / g; // M
    let f = (in_rate * up) as f64; // F
    let lo = in_rate.min(out_rate) as f64;
    let delta = TRANSITION * lo;
    let fc = 0.5 * lo;
    let wp = 2.0 * std::f64::consts::PI * (fc - delta) / f;
    let ws = 2.0 * std::f64::consts::PI * (fc + delta) / f;
    let n = (LENGTH_FACTOR * f / (2.0 * delta)).ceil() as usize;
    let tau = delay_ms * f / 1000.0;
    let wf = 2.0 * std::f64::consts::PI * FAR_STOP * lo / f;
    let h = design(n, wp, ws, wf, tau, up as f64);

    // Phase p of a period (output time p·M on the F grid, the period
    // starting at input 0) takes input j with h index p·M − j·L.
    let mut taps = Vec::with_capacity(up);
    let mut first = Vec::with_capacity(up);
    for p in 0..up {
        let t = (p * down) as isize;
        let l = up as isize;
        let j_max = t.div_euclid(l);
        let j_min = (t - n as isize + 1 + l - 1).div_euclid(l);
        let mut kern: Vec<f64> = (j_min..=j_max).map(|j| h[(t - j * l) as usize]).collect();
        let dc: f64 = kern.iter().sum();
        kern.iter_mut().for_each(|v| *v /= dc);
        taps.push(kern.iter().map(|&v| v as f32).collect());
        first.push(j_min);
    }
    Kernel { taps, first }
}

impl Resampler {
    /// A resampler from `in_rate` to `out_rate` whose passband delay is
    /// `delay_ms`. Equal rates make a delay line of the nearest whole
    /// number of samples.
    pub fn new(in_rate: usize, out_rate: usize, delay_ms: f64) -> Self {
        let kernel = kernel(in_rate, out_rate, delay_ms);
        let g = gcd(in_rate, out_rate);
        let reach = kernel
            .first
            .iter()
            .map(|&f| (-f).max(0) as usize)
            .max()
            .unwrap_or(0)
            + 1;
        Self {
            in_rate,
            phases: out_rate / g,
            in_step: in_rate / g,
            kernel,
            hist: vec![0.0; reach],
        }
    }

    /// The input rate.
    pub fn in_rate(&self) -> usize {
        self.in_rate
    }

    /// Converts `input` (a whole number of rate periods) and appends the
    /// output to `out`.
    pub fn process(&mut self, input: &[f32], out: &mut Vec<f32>) {
        let base = self.hist.len();
        let buf = &mut self.hist;
        buf.extend_from_slice(input);
        let periods = input.len() / self.in_step;
        debug_assert_eq!(periods * self.in_step, input.len());
        let k = &*self.kernel;
        out.reserve(periods * self.phases);
        for per in 0..periods {
            let origin = base as isize + (per * self.in_step) as isize;
            for p in 0..self.phases {
                let start = (origin + k.first[p]) as usize;
                out.push(simd::dot(&k.taps[p], &buf[start..]));
            }
        }
        // Keep the last `base` samples as the next call's history.
        buf.drain(..input.len());
    }

    /// Clears the history (§4.2.9: "re-initialized with silence").
    pub fn reset(&mut self) {
        self.hist.iter_mut().for_each(|v| *v = 0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Every converter the codec builds: the decoder's (each SILK rate to
    /// each output rate with the Table 54 delay) and the encoder's (input
    /// to 48 kHz, 1 ms; 48 kHz to the SILK rate, the rest of its lookahead).
    pub(super) fn converters() -> Vec<(usize, usize, f64)> {
        let mut v = Vec::new();
        for (fin, d) in [(8000, 0.538), (12000, 0.692), (16000, 0.706)] {
            for fout in [8000, 12000, 16000, 24000, 48000] {
                v.push((fin, fout, d));
            }
        }
        for fin in [8000, 12000, 16000, 24000] {
            v.push((fin, 48000, 1.0));
        }
        for (fout, d) in [(8000, 1.337), (12000, 1.2247), (16000, 1.2315)] {
            v.push((48000, fout, d));
        }
        v
    }

    /// The measured response of one converter.
    pub(super) struct Response {
        pub taps: usize,
        pub taps_per_phase: usize,
        /// Largest passband gain deviation, dB, over `[0, fc − Δ]`.
        pub ripple_db: f64,
        /// Largest group-delay deviation from the target, ms, over
        /// `[0, 0.9 (fc − Δ)]`.
        pub delay_err_ms: f64,
        /// Gain at `fc` (the lower rate's Nyquist), dB.
        pub edge_db: f64,
        /// Smallest attenuation over the near stopband `[fc + Δ, 0.75 lo]`
        /// and over the far stopband `[0.75 lo, F/2]`, dB (none for equal
        /// rates).
        pub stop_db: Option<(f64, f64)>,
    }

    /// The response of the polyphase filter as a whole, from its
    /// prototype reassembled out of the branches.
    pub(super) fn measure(fin: usize, fout: usize, d: f64) -> Response {
        let r = Resampler::new(fin, fout, d);
        let up = r.phases;
        let f_hi = (fin * up) as f64;
        let mut pairs = Vec::new();
        for (p, (taps, &first)) in r.kernel.taps.iter().zip(&r.kernel.first).enumerate() {
            let t = (p * r.in_step) as isize;
            for (k, &v) in taps.iter().enumerate() {
                pairs.push((
                    (t - (first + k as isize) * up as isize) as usize,
                    f64::from(v) / up as f64,
                ));
            }
        }
        let mut h = vec![0.0f64; pairs.iter().map(|&(i, _)| i + 1).max().unwrap()];
        for (i, v) in pairs {
            h[i] = v;
        }
        let len = h.len();
        let resp = |f: f64| {
            let w = 2.0 * PI * f / f_hi;
            h.iter().enumerate().fold((0.0, 0.0), |(re, im), (k, &v)| {
                (re + v * (w * k as f64).cos(), im - v * (w * k as f64).sin())
            })
        };
        let gain_db = |f: f64| {
            let (re, im) = resp(f);
            10.0 * (re * re + im * im).log10()
        };
        let phase = |f: f64| {
            let (re, im) = resp(f);
            im.atan2(re)
        };
        let lo = fin.min(fout) as f64;
        let (fc, delta) = (0.5 * lo, TRANSITION * lo);
        let fp = fc - delta;
        let mut ripple: f64 = 0.0;
        let mut delay_err: f64 = 0.0;
        let mut f = 10.0;
        while f <= fp {
            ripple = ripple.max(gain_db(f).abs());
            if f <= 0.9 * fp {
                let df = 1.0;
                let mut dp = phase(f + df) - phase(f);
                if dp > PI {
                    dp -= 2.0 * PI;
                } else if dp < -PI {
                    dp += 2.0 * PI;
                }
                let gd = -dp / (2.0 * PI * df) * 1e3;
                delay_err = delay_err.max((gd - d).abs());
            }
            f += 10.0;
        }
        let stop = (fin != fout).then(|| {
            let (mut near, mut far) = (f64::MAX, f64::MAX);
            let mut f = fc + delta;
            while f < f_hi / 2.0 {
                if f < FAR_STOP * lo {
                    near = near.min(-gain_db(f));
                } else {
                    far = far.min(-gain_db(f));
                }
                f += 5.0;
            }
            (near, far)
        });
        Response {
            taps: len,
            taps_per_phase: r.kernel.taps.iter().map(Vec::len).max().unwrap(),
            ripple_db: ripple,
            delay_err_ms: delay_err,
            edge_db: gain_db(fc),
            stop_db: stop,
        }
    }

    /// Passband flat to ±0.8 dB, group delay within 0.04 ms of the target
    /// over 90 % of the passband, −6 dB ± 1.5 at the lower Nyquist
    /// frequency, image/alias rejection of 30 dB next to the band and 60 dB
    /// beyond.
    #[test]
    fn frequency_response() {
        for (fin, fout, d) in converters() {
            let m = measure(fin, fout, d);
            if fin == fout {
                // A pure delay of the nearest whole number of samples.
                assert!(
                    m.ripple_db < 1e-4 && m.delay_err_ms <= 500.0 / fin as f64,
                    "{fin}: {} {}",
                    m.ripple_db,
                    m.delay_err_ms
                );
                continue;
            }
            assert!(
                m.ripple_db < 0.8,
                "{fin}->{fout}: ripple {:.3} dB",
                m.ripple_db
            );
            assert!(
                m.delay_err_ms < 0.04,
                "{fin}->{fout}: delay error {:.4} ms",
                m.delay_err_ms
            );
            assert!(
                (m.edge_db + 6.0).abs() < 1.5,
                "{fin}->{fout}: {:.2} dB at the edge",
                m.edge_db
            );
            if let Some((near, far)) = m.stop_db {
                assert!(
                    near > 30.0 && far > 60.0,
                    "{fin}->{fout}: stopband {near:.1} / {far:.1} dB"
                );
            }
        }
    }

    /// Prints the table of `docs/VALIDATION.md`.
    #[test]
    #[ignore]
    fn print_response_table() {
        println!(
            "| in → out | delay (ms) | taps (per branch) | passband ripple | delay error | at fc | near stopband | far stopband |"
        );
        println!("|---|---|---|---|---|---|---|---|");
        for (fin, fout, d) in converters() {
            let m = measure(fin, fout, d);
            let stop = m.stop_db.map_or("— | —".to_string(), |(n, f)| {
                format!("{n:.1} dB | {f:.1} dB")
            });
            println!(
                "| {} → {} kHz | {d} | {} ({}) | ±{:.2} dB | {:.3} ms | {:.1} dB | {stop} |",
                fin / 1000,
                fout / 1000,
                m.taps,
                m.taps_per_phase,
                m.ripple_db,
                m.delay_err_ms,
                m.edge_db
            );
        }
    }

    /// Streaming in 10 ms blocks gives the same output as one block, and a
    /// 500 Hz tone comes out delayed by the stated amount (amplitude within
    /// the passband ripple).
    #[test]
    fn streaming_tone() {
        for (fin, fout, d) in converters() {
            let mut r = Resampler::new(fin, fout, d);
            let f = 500.0;
            let n = fin / 5;
            let x: Vec<f32> = (0..n)
                .map(|i| (2.0 * PI * f * i as f64 / fin as f64).sin() as f32)
                .collect();
            let mut y = Vec::new();
            for chunk in x.chunks(fin / 100) {
                r.process(chunk, &mut y);
            }
            assert_eq!(y.len(), n * fout / fin);
            let mut whole = Vec::new();
            Resampler::new(fin, fout, d).process(&x, &mut whole);
            assert!(
                y.iter().zip(&whole).all(|(a, b)| (a - b).abs() < 1e-5),
                "{fin}->{fout}: blocks differ"
            );
            let d = if fin == fout {
                (d * fin as f64 / 1000.0).round() * 1000.0 / fin as f64
            } else {
                d
            };
            let mut err = 0.0f64;
            for (m, &v) in y.iter().enumerate().skip(fout / 25) {
                let want = (2.0 * PI * f * (m as f64 / fout as f64 - d / 1000.0)).sin();
                err = err.max((want - f64::from(v)).abs());
            }
            assert!(err < 0.1, "{fin}->{fout}: max error {err}");
        }
    }
}
