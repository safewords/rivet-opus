//! The CELT decoder (CELT_SPEC §1, §9–§11; RFC 6716 §4.3): frame-level
//! symbol order, anti-collapse, synthesis (energy to amplitude,
//! denormalisation, inverse MDCT with overlap-add, post-filter,
//! de-emphasis) and the decoder state, plus packet loss concealment
//! (RFC 6716 §4.4, non-normative).

use super::Synth;
use super::bands::{self, CollapseMasks, FrameBands, SPREAD_NORMAL};
use super::energy::{self, BandEnergies};
use super::mode::mode;
use super::rate::{self, EncoderChoices};
use super::tables::*;
use crate::range::RangeDecoder;

/// Samples of post-filtered output kept per channel: the longest
/// post-filter period plus its taps, plus a 20 ms frame (CELT_SPEC §11.1),
/// with room for the concealment's pitch search.
const HISTORY: usize = 2048;

/// Shortest post-filter period (CELT_SPEC §10.5).
const MIN_PERIOD: usize = 15;

/// Concealment: lag range of the pitch search, in 48 kHz samples.
const PLC_MIN_LAG: usize = 30;
const PLC_MAX_LAG: usize = 720;
/// Concealment: length of the segment the pitch search matches.
const PLC_WINDOW: usize = 480;

/// Post-filter parameters (CELT_SPEC §10.5).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct PostFilter {
    period: usize,
    gain: f32,
    tapset: usize,
}

/// A CELT decoder (CELT_SPEC §11).
pub struct CeltDecoder {
    /// Output channels (CC).
    channels: usize,
    /// 48000 / output rate.
    downsample: usize,
    /// RFC 8251 §10: ignore the decoded stereo phase inversion.
    pub disable_inv: bool,
    synth: &'static Synth,
    old_band_e: BandEnergies,
    old_log_e: BandEnergies,
    old_log_e2: BandEnergies,
    background_log_e: BandEnergies,
    rng: u32,
    pf: PostFilter,
    pf_old: PostFilter,
    deemph: [f32; 2],
    /// Post-filtered output, newest last.
    history: [Vec<f32>; 2],
    /// The inverse MDCT's tail for the next frame.
    overlap: [Vec<f32>; 2],
    loss_count: u32,
    /// Coded band range of the last good frame (for concealment).
    last_start: usize,
    last_end: usize,
}

impl CeltDecoder {
    /// A decoder producing `channels` output channels at 48000/`downsample`
    /// Hz.
    pub fn new(channels: usize, downsample: usize) -> Self {
        Self {
            channels,
            downsample,
            disable_inv: channels == 1,
            synth: Synth::new(),
            old_band_e: [[0.0; NB_EBANDS]; 2],
            old_log_e: [[-28.0; NB_EBANDS]; 2],
            old_log_e2: [[-28.0; NB_EBANDS]; 2],
            background_log_e: [[0.0; NB_EBANDS]; 2],
            rng: 0,
            pf: PostFilter::default(),
            pf_old: PostFilter::default(),
            deemph: [0.0; 2],
            history: [vec![0.0; HISTORY], vec![0.0; HISTORY]],
            overlap: [vec![0.0; OVERLAP], vec![0.0; OVERLAP]],
            loss_count: 0,
            last_start: 0,
            last_end: NB_EBANDS,
        }
    }

    /// Back to the initial state (CELT_SPEC §11.1); the phase-inversion
    /// setting is kept.
    pub fn reset(&mut self) {
        let disable_inv = self.disable_inv;
        *self = Self::new(self.channels, self.downsample);
        self.disable_inv = disable_inv;
    }

    /// Decodes one frame of `n48` samples at 48 kHz from `ec`, coded with
    /// `c` channels over bands `start..end` (CELT_SPEC §1), writing (or
    /// with `accumulate`, adding) `n48 / downsample` interleaved samples per
    /// output channel to `out`.
    pub fn decode(
        &mut self,
        ec: &mut RangeDecoder,
        n48: usize,
        c: usize,
        start: usize,
        end: usize,
        out: &mut [f32],
        accumulate: bool,
    ) {
        let len = ec.storage();
        if len <= 1 {
            self.decode_lost(n48, out, accumulate);
            return;
        }
        let lm = lm_of(n48);
        let m = 1usize << lm;
        let n = n48;
        let total_bits = (len * 8) as i32;
        if c == 1 {
            for i in 0..NB_EBANDS {
                self.old_band_e[0][i] = self.old_band_e[0][i].max(self.old_band_e[1][i]);
            }
        }
        // §1.2 steps 1–4.
        let tell = ec.tell();
        let silence = if tell >= total_bits {
            true
        } else if tell == 1 {
            ec.bit_logp(15)
        } else {
            false
        };
        if silence {
            ec.add_bits(total_bits - ec.tell());
        }
        let mut pf_new = PostFilter::default();
        if start == 0 && ec.tell() + 16 <= total_bits && ec.bit_logp(1) {
            let octave = ec.uint(6);
            let period = (16 << octave) + ec.bits(4 + octave) as usize - 1;
            let qg = ec.bits(3);
            let tapset = if ec.tell() + 2 <= total_bits {
                ec.icdf(&TAPSET_ICDF, 2)
            } else {
                0
            };
            pf_new = PostFilter {
                period,
                gain: 0.09375 * (qg + 1) as f32,
                tapset,
            };
        }
        let transient = lm > 0 && ec.tell() + 3 <= total_bits && ec.bit_logp(3);
        let intra = ec.tell() + 3 <= total_bits && ec.bit_logp(3);
        // §1.2 steps 5–10.
        let mut scratch = [[0.0f32; NB_EBANDS]; 2];
        let zeros = [[0.0f32; NB_EBANDS]; 2];
        energy::code_coarse(
            ec,
            total_bits,
            &mut self.old_band_e,
            &zeros,
            &mut scratch,
            start,
            end,
            intra,
            c,
            lm,
        );
        let tf_res = bands::code_tf(
            ec,
            start,
            end,
            transient,
            lm,
            total_bits,
            &[0; NB_EBANDS],
            false,
        );
        let spread = if ec.tell() + 4 <= total_bits {
            ec.icdf(&SPREAD_ICDF, 5) as u32
        } else {
            SPREAD_NORMAL
        };
        let cap = rate::init_caps(lm, c);
        let (offsets, boosted_total) =
            rate::code_boosts(ec, start, end, c, lm, &cap, &[0; NB_EBANDS], len);
        let trim = if ec.tell_frac() + (6 << BITRES) <= boosted_total {
            ec.icdf(&TRIM_ICDF, 7) as i32
        } else {
            5
        };
        let mut bits = ((len as i32 * 8) << BITRES) - ec.tell_frac() - 1;
        let anti_collapse_rsv = if transient && lm >= 2 && bits >= (lm as i32 + 2) << BITRES {
            1 << BITRES
        } else {
            0
        };
        bits -= anti_collapse_rsv;
        let choices = EncoderChoices {
            intensity: 0,
            dual_stereo: false,
            prev_coded: 0,
        };
        let alloc =
            rate::compute_allocation(ec, start, end, &offsets, &cap, trim, bits, c, lm, choices);
        // §1.2 steps 11–14.
        energy::code_fine(
            ec,
            &mut self.old_band_e,
            &mut scratch,
            &alloc.fine_quant,
            start,
            end,
            c,
        );
        let mut x = vec![0.0f32; n];
        let mut y = vec![0.0f32; if c == 2 { n } else { 0 }];
        let blocks = if transient { m } else { 1 };
        let params = FrameBands {
            start,
            end,
            lm,
            blocks,
            spread,
            dual_stereo: alloc.dual_stereo,
            intensity: alloc.intensity,
            tf_res: &tf_res,
            shape_total: ((len as i32 * 8) << BITRES) - anti_collapse_rsv,
            balance: alloc.balance,
            pulses: &alloc.pulses,
            coded_bands: alloc.coded_bands,
            disable_inv: self.disable_inv,
        };
        let mut seed = self.rng;
        let collapse = bands::quant_all_bands(
            ec,
            &params,
            &mut x,
            if c == 2 { Some(y.as_mut_slice()) } else { None },
            &zeros,
            &mut seed,
        );
        let anti_collapse_on = anti_collapse_rsv > 0 && ec.bits(1) == 1;
        let left = total_bits - ec.tell();
        energy::code_final(
            ec,
            &mut self.old_band_e,
            &mut scratch,
            &alloc.fine_quant,
            &alloc.fine_priority,
            left,
            start,
            end,
            c,
        );
        if anti_collapse_on {
            let mut chans: Vec<&mut [f32]> = vec![&mut x];
            if c == 2 {
                chans.push(&mut y);
            }
            self.anti_collapse(&mut chans, &collapse, lm, start, end, &alloc.pulses, seed);
        }
        // §10.1–§10.3.
        let mut amp = [[0.0f32; NB_EBANDS]; 2];
        for ch in 0..c {
            for i in start..end {
                let lg = (self.old_band_e[ch][i] + E_MEANS[i]).min(32.0);
                amp[ch][i] = (std::f64::consts::LN_2 * f64::from(lg)).exp() as f32;
            }
        }
        if silence {
            amp = [[0.0; NB_EBANDS]; 2];
            for row in self.old_band_e.iter_mut().take(c) {
                row.fill(-28.0);
            }
        }
        let mut freq: Vec<Vec<f32>> = Vec::with_capacity(2);
        for (ch, xs) in [&x, &y].into_iter().enumerate().take(c) {
            let mut f = vec![0.0f32; n];
            for i in start..end {
                for j in m * EBANDS[i]..m * EBANDS[i + 1] {
                    f[j] = xs[j] * amp[ch][i];
                }
            }
            freq.push(f);
        }
        let bound = if self.downsample == 1 {
            m * EBANDS[end]
        } else {
            (m * EBANDS[end]).min(n / self.downsample)
        };
        for f in &mut freq {
            f[bound.min(n)..].fill(0.0);
        }
        self.synthesize(freq, n, blocks, pf_new, out, accumulate);
        // §11.2 steps 7–11.
        if c == 1 {
            self.old_band_e[1] = self.old_band_e[0];
        }
        if !transient {
            self.old_log_e2 = self.old_log_e;
            self.old_log_e = self.old_band_e;
            for ch in 0..2 {
                for i in 0..NB_EBANDS {
                    self.background_log_e[ch][i] = (self.background_log_e[ch][i]
                        + m as f32 * 0.001)
                        .min(self.old_band_e[ch][i]);
                }
            }
        } else {
            for ch in 0..2 {
                for i in 0..NB_EBANDS {
                    self.old_log_e[ch][i] = self.old_log_e[ch][i].min(self.old_band_e[ch][i]);
                }
            }
        }
        for ch in 0..2 {
            for i in (0..start).chain(end..NB_EBANDS) {
                self.old_band_e[ch][i] = 0.0;
                self.old_log_e[ch][i] = -28.0;
                self.old_log_e2[ch][i] = -28.0;
            }
        }
        self.rng = ec.range();
        self.loss_count = 0;
        self.last_start = start;
        self.last_end = end;
    }

    /// Anti-collapse (CELT_SPEC §9): noise in the short blocks of
    /// transient bands that received no pulses.
    fn anti_collapse(
        &self,
        chans: &mut [&mut [f32]],
        collapse: &CollapseMasks,
        lm: usize,
        start: usize,
        end: usize,
        pulses: &[i32; NB_EBANDS],
        mut seed: u32,
    ) {
        let m = 1usize << lm;
        let c = chans.len();
        for i in start..end {
            let n0 = EBANDS[i + 1] - EBANDS[i];
            let depth = (1 + pulses[i]) / (n0 << lm) as i32;
            let thresh = 0.5 * (-0.125 * depth as f32).exp2();
            let sqrt_1 = 1.0 / ((n0 << lm) as f32).sqrt();
            for ch in 0..c {
                let (mut p1, mut p2) = (self.old_log_e[ch][i], self.old_log_e2[ch][i]);
                if c == 1 {
                    p1 = p1.max(self.old_log_e[1][i]);
                    p2 = p2.max(self.old_log_e2[1][i]);
                }
                let ediff = (self.old_band_e[ch][i] - p1.min(p2)).max(0.0);
                let mut r = 2.0 * (-ediff).exp2();
                if lm == 3 {
                    r *= std::f32::consts::SQRT_2;
                }
                let r = r.min(thresh) * sqrt_1;
                let band = &mut chans[ch][m * EBANDS[i]..m * EBANDS[i + 1]];
                let mut filled = false;
                for k in 0..m {
                    if u32::from(collapse[i][ch]) & 1 << k == 0 {
                        for j in 0..n0 {
                            seed = bands::lcg(seed);
                            band[j * m + k] = if seed & 0x8000 != 0 { r } else { -r };
                        }
                        filled = true;
                    }
                }
                if filled {
                    bands::renormalise(band, 1.0);
                }
            }
        }
    }

    /// From coded-channel spectra to output samples (CELT_SPEC §10.3–§10.6):
    /// channel mapping, inverse MDCT with overlap-add, post-filter,
    /// de-emphasis and decimation. Updates the post-filter state.
    fn synthesize(
        &mut self,
        mut freq: Vec<Vec<f32>>,
        n: usize,
        blocks: usize,
        pf_new: PostFilter,
        out: &mut [f32],
        accumulate: bool,
    ) {
        let cc = self.channels;
        if cc == 2 && freq.len() == 1 {
            freq.push(freq[0].clone());
        } else if cc == 1 && freq.len() == 2 {
            let f1 = freq.pop().unwrap_or_default();
            for (a, b) in freq[0].iter_mut().zip(&f1) {
                *a = 0.5 * (*a + b);
            }
        }
        // The textbook inverse MDCT of CELT_SPEC §10.4. The gain
        // 1 + (π/(8n))² that section adds is left out: with it the CELT
        // test vectors lose about 10 dB of SNR (106 → 96 dB on vector 01).
        let lm = lm_of(n);
        self.pf.period = self.pf.period.max(MIN_PERIOD);
        self.pf_old.period = self.pf_old.period.max(MIN_PERIOD);
        let mut buf = vec![0.0f32; n + OVERLAP];
        for (ch, f) in freq.iter().enumerate() {
            buf[..OVERLAP].copy_from_slice(&self.overlap[ch]);
            self.synth.imdct_ola(f, &mut buf, n, blocks);
            self.overlap[ch].copy_from_slice(&buf[n..n + OVERLAP]);
            let h = &mut self.history[ch];
            h.copy_within(n.., 0);
            h[HISTORY - n..].copy_from_slice(&buf[..n]);
            let base = HISTORY - n;
            comb_filter(h, base, OVERLAP, self.pf_old, self.pf, OVERLAP);
            if lm != 0 {
                comb_filter(h, base + OVERLAP, n - OVERLAP, self.pf, pf_new, OVERLAP);
            }
            let mut mem = self.deemph[ch];
            let ds = self.downsample;
            for j in 0..n {
                let t = h[base + j] + mem;
                mem = PREEMPH_COEF * t;
                if j % ds == 0 {
                    let o = &mut out[(j / ds) * cc + ch];
                    let v = t * (1.0 / 32768.0);
                    if accumulate {
                        *o += v;
                    } else {
                        *o = v;
                    }
                }
            }
            self.deemph[ch] = mem;
        }
        self.pf_old = self.pf;
        self.pf = pf_new;
        if lm != 0 {
            self.pf_old = pf_new;
        }
    }

    /// Packet loss concealment for one frame of `n48` samples (RFC 6716
    /// §4.4; CELT_SPEC §11.4 describes what later frames depend on). Our
    /// own design: after a good CELT-only frame the output is continued
    /// periodically at the pitch found in the history, fading frame by
    /// frame, and fed through the regular synthesis (forward MDCT of the
    /// continuation, then inverse MDCT with overlap-add) so the overlap
    /// with the previous and the next frame cancels its aliasing as in
    /// normal decoding. After five losses in a row, or after a hybrid
    /// frame, it is noise shaped by decaying band energies.
    pub fn decode_lost(&mut self, n48: usize, out: &mut [f32], accumulate: bool) {
        let n = n48;
        let lm = lm_of(n);
        let m = 1usize << lm;
        let (start, end) = (self.last_start, self.last_end);
        let pf = self.pf;
        let freq: Vec<Vec<f32>> = if self.loss_count >= 5 || start != 0 {
            let decay = if self.loss_count == 0 { 1.5 } else { 0.5 };
            let mut freq = Vec::with_capacity(self.channels);
            for ch in 0..self.channels {
                let mut f = vec![0.0f32; n];
                for i in start..end {
                    let e = &mut self.old_band_e[ch][i];
                    *e = (*e - decay).max(-28.0);
                    let amp = (*e + E_MEANS[i]).min(32.0).exp2();
                    let band = &mut f[m * EBANDS[i]..m * EBANDS[i + 1]];
                    for v in band.iter_mut() {
                        self.rng = bands::lcg(self.rng);
                        *v = (self.rng as i32 >> 20) as f32;
                    }
                    bands::renormalise(band, amp);
                }
                let bound = (m * EBANDS[end]).min(n / self.downsample);
                f[bound..].fill(0.0);
                freq.push(f);
            }
            freq
        } else {
            let lag = self.pitch_lag();
            // Periodic signals fade slowly, others quickly.
            let fade = if lag.1 > 0.6 { 0.85f32 } else { 0.5 };
            let mut freq = Vec::with_capacity(self.channels);
            for ch in 0..self.channels {
                let h = &self.history[ch];
                let total = n + OVERLAP;
                // The continuation of the post-filtered output…
                let mut cont = vec![0.0f32; total];
                for j in 0..total {
                    let src = HISTORY as isize + j as isize - lag.0 as isize;
                    cont[j] = if src < HISTORY as isize {
                        h[src as usize]
                    } else {
                        cont[src as usize - HISTORY]
                    };
                }
                for (j, v) in cont.iter_mut().enumerate() {
                    *v *= fade.powf(((j + 1).min(n) as f32) / n as f32);
                }
                // …taken back before the post-filter, which synthesis
                // applies again with the current parameters.
                let at = |j: isize| {
                    if j < 0 {
                        h[(HISTORY as isize + j) as usize]
                    } else {
                        cont[j as usize]
                    }
                };
                let gains = COMB_FILTER_GAINS[pf.tapset];
                let t = pf.period.max(MIN_PERIOD) as isize;
                let mut pre = cont.clone();
                if pf.gain != 0.0 {
                    for (j, p) in pre.iter_mut().enumerate() {
                        let j = j as isize;
                        *p -= pf.gain
                            * (gains[0] * at(j - t)
                                + gains[1] * (at(j - t - 1) + at(j - t + 1))
                                + gains[2] * (at(j - t - 2) + at(j - t + 2)));
                    }
                }
                freq.push(self.synth.mdct_windowed(&pre, n));
            }
            // The next frame predicts its energies from these: never let
            // them exceed what was concealed, so a note that ended during
            // the loss does not come back too loud. An error in the state
            // decays by the inter-frame prediction coefficient per frame,
            // slowly for short frames, so those keep a further margin.
            let margin = PRED_COEF[lm] - 0.5;
            for (ch, f) in freq.iter().enumerate() {
                for i in start..end {
                    let e: f32 = f[m * EBANDS[i]..m * EBANDS[i + 1]]
                        .iter()
                        .map(|v| v * v)
                        .sum();
                    let lg = 0.5 * (e + 1e-15).log2() - E_MEANS[i] - margin;
                    for row in self.old_band_e.iter_mut().skip(ch).step_by(self.channels) {
                        row[i] = row[i].min(lg.max(-28.0));
                    }
                }
            }
            freq
        };
        self.synthesize(freq, n, 1, pf, out, accumulate);
        self.loss_count += 1;
    }

    /// The pitch lag of the recent output (normalised cross-correlation of
    /// the last [`PLC_WINDOW`] samples with the past, channels summed) and
    /// its correlation.
    fn pitch_lag(&self) -> (usize, f32) {
        let mono: Vec<f32> = (0..HISTORY)
            .map(|k| self.history[..self.channels].iter().map(|h| h[k]).sum())
            .collect();
        let seg = &mono[HISTORY - PLC_WINDOW..];
        let e_seg: f32 = seg.iter().map(|v| v * v).sum();
        let mut best = (PLC_MAX_LAG, 0.0f32);
        if e_seg <= 1e-9 {
            return best;
        }
        for lag in PLC_MIN_LAG..=PLC_MAX_LAG {
            let past = &mono[HISTORY - PLC_WINDOW - lag..HISTORY - lag];
            let (mut xy, mut yy) = (0.0f32, 0.0f32);
            for (a, b) in seg.iter().zip(past) {
                xy += a * b;
                yy += b * b;
            }
            let r = xy / (e_seg * yy).sqrt().max(1e-9);
            if r > best.1 {
                best = (lag, r);
            }
        }
        best
    }
}

/// LM of a frame of `n` samples at 48 kHz (CELT_SPEC §1.1).
fn lm_of(n: usize) -> usize {
    match n {
        120 => 0,
        240 => 1,
        480 => 2,
        _ => 3,
    }
}

/// The post-filter comb of CELT_SPEC §10.5, in place on `len` samples of
/// `h` from `base`, cross-fading over the first `ov` samples from `from`
/// to `to`.
fn comb_filter(
    h: &mut [f32],
    base: usize,
    len: usize,
    from: PostFilter,
    to: PostFilter,
    ov: usize,
) {
    let window = &mode().window;
    let a = COMB_FILTER_GAINS[from.tapset].map(|g| from.gain * g);
    let b = COMB_FILTER_GAINS[to.tapset].map(|g| to.gain * g);
    let (t0, t1) = (from.period, to.period);
    for i in 0..len {
        let k = base + i;
        let mut acc = h[k];
        let f = if i < ov { window[i] * window[i] } else { 1.0 };
        if i < ov && from.gain != 0.0 {
            let w = 1.0 - f;
            acc += (w * a[0]) * h[k - t0]
                + (w * a[1]) * (h[k - t0 - 1] + h[k - t0 + 1])
                + (w * a[2]) * (h[k - t0 - 2] + h[k - t0 + 2]);
        }
        if to.gain != 0.0 {
            acc += (f * b[0]) * h[k - t1]
                + (f * b[1]) * (h[k - t1 - 1] + h[k - t1 + 1])
                + (f * b[2]) * (h[k - t1 - 2] + h[k - t1 + 2]);
        }
        h[k] = acc;
    }
}
