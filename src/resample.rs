//! Sample-rate conversion between the SILK rates (8, 12, 16 kHz) and the
//! Opus API rates (8, 12, 16, 24, 48 kHz).
//!
//! RFC 6716 §4.2.9 leaves the resampler non-normative but fixes the delay it
//! may add (Table 54). From the SILK rates to 48 kHz the filters are the
//! ones fitted to the reference decoder's output (see `resample_fit`); for
//! every other pair this is a linear-phase windowed-sinc interpolator whose
//! group delay is exactly the Table 54 allocation, so SILK output lines up
//! with the CELT layer as the encoder expects.

use crate::resample_fit::{FIT_8K, FIT_12K, FIT_16K};

/// A streaming resampler for one channel.
pub(crate) struct Resampler {
    in_rate: usize,
    out_rate: usize,
    /// Output samples per period of the rate ratio, and input samples.
    phases: usize,
    in_step: usize,
    /// Kernel taps per phase, starting at input offset `first[p]` relative
    /// to the period's first input sample.
    taps: Vec<Vec<f32>>,
    first: Vec<isize>,
    /// Input history, oldest first.
    hist: Vec<f32>,
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 { a } else { gcd(b, a % b) }
}

impl Resampler {
    /// A resampler from `in_rate` to `out_rate` delaying by `delay_ms`.
    pub fn new(in_rate: usize, out_rate: usize, delay_ms: f64) -> Self {
        let g = gcd(in_rate, out_rate);
        let phases = out_rate / g;
        let in_step = in_rate / g;
        if in_rate == out_rate {
            // A plain copy.
            return Self { in_rate, out_rate, phases: 1, in_step: 1, taps: vec![vec![1.0]], first: vec![0], hist: vec![0.0; 1] };
        }

        // Delay in input samples; the kernel's half-width equals it so the
        // filter is causal.
        let delay = delay_ms * in_rate as f64 / 1000.0;
        let half = delay;
        let cutoff = 0.5 * (out_rate.min(in_rate) as f64 / in_rate as f64) * 0.97;
        let mut taps = Vec::with_capacity(phases);
        let mut first = Vec::with_capacity(phases);
        for p in 0..phases {
            // Output p of a period sits at input time p * in/out - delay.
            let t = p as f64 * in_rate as f64 / out_rate as f64 - delay;
            let k0 = (t - half).floor() as isize + 1;
            let k1 = (t + half).ceil() as isize - 1;
            let mut kern = Vec::new();
            let mut sum = 0.0;
            for k in k0..=k1 {
                let u = t - k as f64;
                let x = 2.0 * cutoff * u;
                let sinc = if x.abs() < 1e-12 { 1.0 } else { (std::f64::consts::PI * x).sin() / (std::f64::consts::PI * x) };
                let wnd = {
                    let r = u / half;
                    if r.abs() >= 1.0 { 0.0 } else { 0.5 + 0.5 * (std::f64::consts::PI * r).cos() }
                };
                let v = 2.0 * cutoff * sinc * wnd;
                sum += v;
                kern.push(v);
            }
            // Unity gain at DC.
            let kern: Vec<f32> = kern.iter().map(|v| (v / sum) as f32).collect();
            taps.push(kern);
            first.push(k0);
        }
        let reach = first.iter().map(|&f| (-f).max(0) as usize).max().unwrap_or(0) + 1;
        Self { in_rate, out_rate, phases, in_step, taps, first, hist: vec![0.0; reach] }
    }

    /// The decoder's SILK output resampler: the fitted filters to 48 kHz,
    /// otherwise the windowed-sinc design with the Table 54 delay.
    pub fn silk_output(in_rate: usize, out_rate: usize, delay_ms: f64) -> Self {
        let fitted = match (in_rate, out_rate) {
            (8000, 48000) => Some(FIT_8K),
            (12000, 48000) => Some(FIT_12K),
            (16000, 48000) => Some(FIT_16K),
            _ => None,
        };
        match fitted {
            Some(fit) => {
                let k = fit[0].len();
                let taps: Vec<Vec<f32>> = fit.iter().map(|g| g.iter().rev().copied().collect()).collect();
                let first = vec![-(k as isize - 1); fit.len()];
                Self { in_rate, out_rate, phases: fit.len(), in_step: 1, taps, first, hist: vec![0.0; k] }
            }
            None => Self::new(in_rate, out_rate, delay_ms),
        }
    }

    /// The input rate.
    pub fn in_rate(&self) -> usize {
        self.in_rate
    }

    /// The output rate.
    pub fn out_rate(&self) -> usize {
        self.out_rate
    }

    /// Converts `input` (a whole number of rate periods) and appends the
    /// output to `out`.
    pub fn process(&mut self, input: &[f32], out: &mut Vec<f32>) {
        let base = self.hist.len();
        let mut buf = std::mem::take(&mut self.hist);
        buf.extend_from_slice(input);
        let periods = input.len() / self.in_step;
        debug_assert_eq!(periods * self.in_step, input.len());
        for per in 0..periods {
            let origin = base as isize + (per * self.in_step) as isize;
            for p in 0..self.phases {
                let start = origin + self.first[p];
                let mut acc = 0.0f32;
                for (k, &h) in self.taps[p].iter().enumerate() {
                    let idx = start + k as isize;
                    if idx >= 0 && (idx as usize) < buf.len() {
                        acc += buf[idx as usize] * h;
                    }
                }
                out.push(acc);
            }
        }
        let keep = base;
        self.hist = buf[buf.len() - keep..].to_vec();
    }

    /// Clears the history (§4.2.9: "re-initialized with silence").
    pub fn reset(&mut self) {
        self.hist.iter_mut().for_each(|v| *v = 0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tone keeps its amplitude and is delayed by the stated amount.
    #[test]
    fn tone_through_resampler() {
        for (fin, fout, d) in [(16000, 48000, 0.706), (8000, 48000, 0.538), (12000, 48000, 0.692), (16000, 8000, 0.706), (12000, 24000, 0.692)] {
            let mut r = Resampler::new(fin, fout, d);
            let f = 440.0;
            let n = fin / 10;
            let x: Vec<f32> = (0..n).map(|i| (2.0 * std::f64::consts::PI * f * i as f64 / fin as f64).sin() as f32).collect();
            let mut y = Vec::new();
            for chunk in x.chunks(fin / 100) {
                r.process(chunk, &mut y);
            }
            assert_eq!(y.len(), n * fout / fin);
            let delay = d / 1000.0;
            let mut err = 0.0f64;
            let mut cnt = 0;
            for (m, &v) in y.iter().enumerate().skip(fout / 50) {
                let t = m as f64 / fout as f64 - delay;
                let want = (2.0 * std::f64::consts::PI * f * t).sin();
                err = err.max((want - f64::from(v)).abs());
                cnt += 1;
            }
            assert!(cnt > 0 && err < 0.02, "{fin}->{fout}: max error {err}");
        }
    }
}
