//! The Opus conformance quality metric (RFC 6716 §6, §6.1).
//!
//! RFC 6716 §6 defines decoder compliance as reproducing the reference
//! decoder's final range state on every packet *and* producing output that
//! is "within the thresholds specified by the opus_compare.c tool" against
//! the reference output of each test vector, at every output rate and
//! channel count supported. §6.1 describes the result — a quality of 100
//! for identical output, 0 at the pass threshold, calibrated to additive
//! white noise at 48 dB SNR — but the metric itself is given only by that
//! tool, which is part of the normative Appendix A. This is a
//! re-implementation of the metric from its specification there, not a
//! translation of the tool; it is used only by the tests.
//!
//! The metric, as specified:
//!
//! * Both signals are 16-bit PCM in sample units. The reference is the
//!   48 kHz stereo `.dec` file; for a mono comparison it is downmixed to
//!   `(L + R) / 2`. The tested signal is at the output rate `fs`, with
//!   `D = 48000 / fs`; its length must be the reference's divided by `D`.
//! * Frames: Hann windows (`0.5 − 0.5·cos(2πk/(W−1))`) of `W = 480 / D`
//!   samples every `120 / D` samples, `(len48 − 480 + 120) / 120` of them;
//!   every frame's DFT bins are 100 Hz apart at every rate. The power of bin
//!   `j` is `|D · Σ w[k]·x[k]·e^{−2πijk/W}|² + 100000`.
//! * 21 bands on bin edges 0, 2, 4, …, 156, 200 (the CELT bands in 100 Hz
//!   units); at 8, 12, 16 and 24 kHz only the first 13, 15, 17 and 19.
//!   The reference's mean bin power per band is spread into a masking
//!   level: upwards across bands (+0.1 of the band below, cumulatively),
//!   downwards (+0.03 of the band above), forwards in time (+0.5 of the
//!   same band in the previous frame, once masked) and, in stereo, across
//!   channels (+0.01 of the other channel's level). A tenth of the mask is
//!   added to every bin of both spectra.
//! * Each bin's power is then summed with its own (unmasked-sum) value in
//!   the previous frame.
//! * Per bin and channel the distortion is `r − ln r − 1` with
//!   `r = P_test / P_ref`, weighted by 0.1 on bins 79–81 and by 0.01 on bin
//!   80 (the SILK/CELT cross-over at 8 kHz). Below 48 kHz the top 300 Hz
//!   are ignored (none at 12 kHz, whose last band already stops 400 Hz
//!   short). A band's distortion is the mean over its bins and channels;
//!   a frame's is the mean of the bands' squares over 21 (a fixed count),
//!   raised to the 4th power; the error is the 16th root of the mean of the
//!   frames' 4th powers, and `Q = 100·(1 − ln(1 + err) / (2·ln 1.13))`.
//!   The vector passes at `Q ≥ 0`.

#![allow(dead_code)]

/// Band edges in 100 Hz bins.
const BANDS: [usize; 22] = [
    0, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 68, 80, 96, 120, 156, 200,
];
const NBANDS: usize = 21;
/// Bins per frame of the 48 kHz analysis that are stored (W/2).
const NFREQS: usize = 240;
const WIN: usize = 480;
const STEP: usize = 120;

/// The spectral analysis of one signal: per frame, per bin, per channel
/// power, and per band mean power.
pub struct Spectrum {
    pub frames: usize,
    pub channels: usize,
    pub freqs: usize,
    /// `[frame][bin][channel]`, `freqs` bins per frame.
    pub power: Vec<f32>,
    /// `[frame][band][channel]`, `NBANDS` per frame (only those analysed
    /// are filled).
    pub band: Vec<f32>,
}

/// Bins up to `BANDS[nbands]` of `frames` Hann-windowed frames of `x`
/// (interleaved, `channels`), window `win`, hop `step`, scaled by `scale`.
fn analyse(
    x: &[f32],
    channels: usize,
    frames: usize,
    win: usize,
    step: usize,
    scale: f32,
    nbands: usize,
) -> Spectrum {
    let freqs = NFREQS * win / WIN;
    let window: Vec<f32> = (0..win)
        .map(|k| {
            0.5 - 0.5 * (2.0 * std::f64::consts::PI * k as f64 / (win - 1) as f64).cos() as f32
        })
        .collect();
    let cos: Vec<f32> = (0..win)
        .map(|k| (2.0 * std::f64::consts::PI * k as f64 / win as f64).cos() as f32)
        .collect();
    let sin: Vec<f32> = (0..win)
        .map(|k| (2.0 * std::f64::consts::PI * k as f64 / win as f64).sin() as f32)
        .collect();
    let top = BANDS[nbands];
    let mut power = vec![0.0f32; frames * freqs * channels];
    let mut band = vec![0.0f32; frames * NBANDS * channels];
    let mut seg = vec![0.0f32; win];
    for f in 0..frames {
        for c in 0..channels {
            for (k, s) in seg.iter_mut().enumerate() {
                *s = window[k] * x[(f * step + k) * channels + c];
            }
            for j in 0..top {
                let (mut re, mut im) = (0.0f32, 0.0f32);
                let mut t = 0usize;
                for &s in &seg {
                    re += cos[t] * s;
                    im -= sin[t] * s;
                    t += j;
                    if t >= win {
                        t -= win;
                    }
                }
                re *= scale;
                im *= scale;
                power[(f * freqs + j) * channels + c] = re * re + im * im + 100_000.0;
            }
            for b in 0..nbands {
                let sum: f32 = (BANDS[b]..BANDS[b + 1])
                    .map(|j| power[(f * freqs + j) * channels + c])
                    .sum();
                band[(f * NBANDS + b) * channels + c] = sum / (BANDS[b + 1] - BANDS[b]) as f32;
            }
        }
    }
    Spectrum {
        frames,
        channels,
        freqs,
        power,
        band,
    }
}

/// The reference's analysis: `reference` is the 48 kHz stereo `.dec`
/// signal in 16-bit units; `channels` 1 downmixes it.
pub fn reference_spectrum(reference: &[f32], channels: usize) -> Spectrum {
    let x: Vec<f32> = if channels == 1 {
        reference
            .chunks_exact(2)
            .map(|c| 0.5 * (c[0] + c[1]))
            .collect()
    } else {
        reference.to_vec()
    };
    let len = x.len() / channels;
    assert!(len >= WIN, "reference too short");
    let frames = (len - WIN + STEP) / STEP;
    analyse(&x, channels, frames, WIN, STEP, 1.0, NBANDS)
}

/// The quality `Q` of `test` (interleaved, `channels`, at `rate`, in 16-bit
/// units) against the reference analysed by [`reference_spectrum`]
/// with the same channel count. `None` if the lengths do not match.
pub fn quality(reference: &Spectrum, test: &[f32], rate: u32, channels: usize) -> Option<f64> {
    assert_eq!(reference.channels, channels);
    let ds = (48000 / rate) as usize;
    let ybands = match rate {
        8000 => 13,
        12000 => 15,
        16000 => 17,
        24000 => 19,
        _ => NBANDS,
    };
    let frames = reference.frames;
    let len48 = (frames - 1) * STEP + WIN;
    let ylen = test.len() / channels;
    // The reference length implied by the frame count may round down; the
    // lengths themselves are checked by the caller against the file.
    if ylen * ds < len48 {
        return None;
    }
    let y = analyse(
        test,
        channels,
        frames,
        WIN / ds,
        STEP / ds,
        ds as f32,
        ybands,
    );
    let cc = channels;
    let xf = NFREQS;
    let yf = y.freqs;
    let mut xp = reference.power.clone();
    let mut yp = y.power;
    let mut xb = reference.band.clone();

    for f in 0..frames {
        let at = |b: usize, c: usize| (f * NBANDS + b) * cc + c;
        for b in 1..NBANDS {
            for c in 0..cc {
                xb[at(b, c)] += 0.1 * xb[at(b - 1, c)];
            }
        }
        for b in (0..NBANDS - 1).rev() {
            for c in 0..cc {
                xb[at(b, c)] += 0.03 * xb[at(b + 1, c)];
            }
        }
        if f > 0 {
            for b in 0..NBANDS {
                for c in 0..cc {
                    xb[at(b, c)] += 0.5 * xb[((f - 1) * NBANDS + b) * cc + c];
                }
            }
        }
        if cc == 2 {
            for b in 0..NBANDS {
                let (l, r) = (xb[at(b, 0)], xb[at(b, 1)]);
                xb[at(b, 0)] += 0.01 * r;
                xb[at(b, 1)] += 0.01 * l;
            }
        }
        for b in 0..ybands {
            for j in BANDS[b]..BANDS[b + 1] {
                for c in 0..cc {
                    let m = 0.1 * xb[at(b, c)];
                    xp[(f * xf + j) * cc + c] += m;
                    yp[(f * yf + j) * cc + c] += m;
                }
            }
        }
    }

    // Each frame's bin plus the previous frame's (pre-sum) value.
    for j in 0..BANDS[ybands] {
        for c in 0..cc {
            let (mut px, mut py) = (xp[j * cc + c], yp[j * cc + c]);
            for f in 1..frames {
                let (ix, iy) = ((f * xf + j) * cc + c, (f * yf + j) * cc + c);
                let (ox, oy) = (xp[ix], yp[iy]);
                xp[ix] += px;
                yp[iy] += py;
                px = ox;
                py = oy;
            }
        }
    }

    let max_compare = match rate {
        48000 => BANDS[NBANDS],
        12000 => BANDS[ybands],
        _ => BANDS[ybands] - 3,
    };
    let mut err = 0.0f64;
    for f in 0..frames {
        let mut ef = 0.0f64;
        for b in 0..ybands {
            let mut eb = 0.0f64;
            for j in BANDS[b]..BANDS[b + 1].min(max_compare) {
                for c in 0..cc {
                    let r = yp[(f * yf + j) * cc + c] / xp[(f * xf + j) * cc + c];
                    let mut d = r - r.ln() - 1.0;
                    if (79..=81).contains(&j) {
                        d *= 0.1;
                    }
                    if j == 80 {
                        d *= 0.1;
                    }
                    eb += f64::from(d);
                }
            }
            eb /= ((BANDS[b + 1] - BANDS[b]) * cc) as f64;
            ef += eb * eb;
        }
        ef /= NBANDS as f64;
        ef *= ef;
        err += ef * ef;
    }
    let err = (err / frames as f64).powf(1.0 / 16.0);
    Some(100.0 * (1.0 - 0.5 * (1.0 + err).ln() / 1.13f64.ln()))
}

/// Rounds and clamps decoder output (±1.0 full scale) to 16-bit units, as
/// a 16-bit PCM file of it would hold.
pub fn to_pcm16_units(x: &[f32]) -> Vec<f32> {
    x.iter()
        .map(|&v| (v * 32768.0).round().clamp(-32768.0, 32767.0))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Identical signals score exactly 100, the quality falls as noise is
    /// added, and noise as loud as the signal fails.
    #[test]
    fn identity_and_monotonicity() {
        let n = 48000;
        let mut s = 1u32;
        let mut rnd = move || {
            s = s.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            (s >> 8) as f32 / (1u32 << 23) as f32 - 1.0
        };
        // A few tones: a sparse spectrum, as music has.
        let x: Vec<f32> = (0..2 * n)
            .map(|i| {
                let t = (i / 2) as f32 / 48000.0;
                [220.0f32, 1100.0, 5300.0]
                    .iter()
                    .map(|f| 6000.0 * (2.0 * std::f32::consts::PI * f * t).sin())
                    .sum::<f32>()
                    .round()
            })
            .collect();
        let spec = reference_spectrum(&x, 2);
        assert!((quality(&spec, &x, 48000, 2).unwrap() - 100.0).abs() < 1e-9);
        let noise: Vec<f32> = (0..2 * n).map(|_| rnd()).collect();
        let mut last = 100.0;
        for amp in [10.0f32, 100.0, 1000.0, 10000.0] {
            let y: Vec<f32> = x
                .iter()
                .zip(&noise)
                .map(|(a, b)| (a + amp * b).round())
                .collect();
            let q = quality(&spec, &y, 48000, 2).unwrap();
            assert!(q < last, "noise {amp}: Q {q} not below {last}");
            last = q;
        }
        assert!(last < 0.0, "noise as loud as the signal passes: Q = {last}");
    }
}
