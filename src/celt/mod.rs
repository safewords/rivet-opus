//! The CELT layer (RFC 6716 §4.3 and §5.3).

pub(crate) mod bands;
pub(crate) mod cwrs;
pub(crate) mod decoder;
pub(crate) mod encoder;
pub(crate) mod energy;
pub(crate) mod mode;
pub(crate) mod rate;
pub(crate) mod tables;

use crate::mdct::Mdct;
use tables::OVERLAP;

/// Gain of the synthesis (inverse) MDCT relative to the textbook transform
/// (the RFC describes the output as "scaled by 1/2" in the reference
/// implementation's own transform convention; measured against the test
/// vectors, the textbook inverse needs no extra factor).
const SYNTH_SCALE: f32 = 1.0;

/// The MDCTs of the four block sizes and the window, shared by synthesis
/// and analysis.
pub(crate) struct Synth {
    mdcts: [Mdct; 4],
}

impl Synth {
    pub fn new() -> Self {
        Self { mdcts: [Mdct::new(120), Mdct::new(240), Mdct::new(480), Mdct::new(960)] }
    }

    fn mdct(&self, n: usize) -> &Mdct {
        &self.mdcts[match n {
            120 => 0,
            240 => 1,
            480 => 2,
            _ => 3,
        }]
    }

    /// The window over a block of `nb` coefficients (2·nb samples): zero,
    /// rising over the overlap, one, falling, zero.
    fn window_at(nb: usize, t: usize) -> f32 {
        let w = &mode::mode().window;
        let rise = nb / 2 - OVERLAP / 2;
        let fall = 3 * nb / 2 - OVERLAP / 2;
        if t < rise || t >= fall + OVERLAP {
            0.0
        } else if t < rise + OVERLAP {
            w[t - rise]
        } else if t < fall {
            1.0
        } else {
            w[OVERLAP - 1 - (t - fall)]
        }
    }

    /// Inverse MDCTs of `blocks` interleaved blocks of `freq` (a frame of `n`
    /// coefficients), windowed and overlap-added into `out`, whose first
    /// `OVERLAP` samples hold the previous frame's tail; `out[OVERLAP ..
    /// n + OVERLAP]` is overwritten.
    pub fn imdct_ola(&self, freq: &[f32], out: &mut [f32], n: usize, blocks: usize) {
        let nb = n / blocks;
        let mdct = self.mdct(nb);
        out[OVERLAP..n + OVERLAP].fill(0.0);
        let rise = nb / 2 - OVERLAP / 2;
        let mut coefs = vec![0.0f32; nb];
        let mut y = vec![0.0f32; 2 * nb];
        for b in 0..blocks {
            for (k, c) in coefs.iter_mut().enumerate() {
                *c = freq[b + k * blocks];
            }
            mdct.inverse(&coefs, &mut y);
            for t in rise..rise + nb + OVERLAP {
                out[b * nb + t - rise] += SYNTH_SCALE * y[t] * Self::window_at(nb, t);
            }
        }
    }

    /// The forward MDCTs of `blocks` blocks covering `x[.. n + OVERLAP]`
    /// (block `b` takes `x[b·nb .. b·nb + nb + OVERLAP]`), interleaved, with
    /// the scaling that makes [`Self::imdct_ola`] its inverse.
    pub fn mdct_blocks(&self, x: &[f32], n: usize, blocks: usize) -> Vec<f32> {
        let nb = n / blocks;
        let mdct = self.mdct(nb);
        let rise = nb / 2 - OVERLAP / 2;
        let scale = 2.0 / nb as f32 / SYNTH_SCALE;
        let mut out = vec![0.0f32; n];
        let mut buf = vec![0.0f32; 2 * nb];
        let mut coefs = vec![0.0f32; nb];
        for b in 0..blocks {
            buf.fill(0.0);
            for t in rise..rise + nb + OVERLAP {
                buf[t] = x[b * nb + t - rise] * Self::window_at(nb, t);
            }
            mdct.forward(&buf, &mut coefs);
            for (k, c) in coefs.iter().enumerate() {
                out[b + k * blocks] = c * scale;
            }
        }
        out
    }

    /// One long-block MDCT of `x[.. n + OVERLAP]`.
    pub fn mdct_windowed(&self, x: &[f32], n: usize) -> Vec<f32> {
        self.mdct_blocks(x, n, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Analysis then synthesis with overlap-add gives the input back.
    #[test]
    fn analysis_synthesis_is_identity() {
        let s = Synth::new();
        for (n, blocks) in [(960, 1), (960, 8), (480, 4), (120, 1), (240, 2)] {
            let total = 6 * n + OVERLAP;
            let x: Vec<f32> = (0..total).map(|i| ((i * 7919) % 1000) as f32 / 500.0 - 1.0).collect();
            let mut out = vec![0.0f32; total + n];
            for f in 0..6 {
                let coefs = s.mdct_blocks(&x[f * n..f * n + n + OVERLAP], n, blocks);
                s.imdct_ola(&coefs, &mut out[f * n..], n, blocks);
            }
            for i in OVERLAP..5 * n {
                assert!((out[i] - x[i]).abs() < 1e-3, "n {n} b {blocks} i {i}: {} vs {}", out[i], x[i]);
            }
        }
    }
}
