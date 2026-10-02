//! The CELT decoder (RFC 6716 §4.3): frame decoding, synthesis (inverse
//! MDCT with the low-overlap window), the pitch post-filter, de-emphasis,
//! and packet loss concealment.

use super::bands::{self, FrameParams, SPREAD_NORMAL};
use super::energy;
use super::mode::{BITRES, init_caps, mode};
use super::rate::{self, EncoderChoices};
use super::tables::*;
use super::Synth;
use crate::range::RangeDecoder;

/// Samples of output history kept per channel (for the post-filter and the
/// PLC).
pub const DECODE_BUFFER: usize = 2048;

/// What a decoded CELT frame tells the Opus layer.
#[derive(Debug, Clone, Copy, Default)]
pub struct FrameInfo {
    /// The post-filter pitch period, if the post-filter is on.
    pub pitch: Option<u32>,
}

/// A CELT decoder for one Opus stream (one or two channels).
pub struct CeltDecoder {
    /// Output channels.
    channels: usize,
    /// 48000 / output rate.
    downsample: usize,
    decode_mem: Vec<Vec<f32>>,
    pub(crate) old_band_e: [f32; 2 * NB_EBANDS],
    old_log_e: [f32; 2 * NB_EBANDS],
    old_log_e2: [f32; 2 * NB_EBANDS],
    postfilter_period: usize,
    postfilter_period_old: usize,
    postfilter_gain: f32,
    postfilter_gain_old: f32,
    postfilter_tapset: usize,
    postfilter_tapset_old: usize,
    preemph_mem: [f32; 2],
    /// The final range of the last frame, also the folding seed.
    pub(crate) rng: u32,
    loss_count: u32,
    last_pitch: usize,
    /// Do not invert the phase of intensity-stereo bands (RFC 8251 §10).
    pub disable_inv: bool,
    synth: Synth,
}

impl CeltDecoder {
    /// A decoder with `channels` output channels at `48000 / downsample` Hz.
    pub fn new(channels: usize, downsample: usize) -> Self {
        let mut d = Self {
            channels,
            downsample,
            decode_mem: vec![vec![0.0; DECODE_BUFFER + OVERLAP]; channels],
            old_band_e: [0.0; 2 * NB_EBANDS],
            old_log_e: [-28.0; 2 * NB_EBANDS],
            old_log_e2: [-28.0; 2 * NB_EBANDS],
            postfilter_period: 0,
            postfilter_period_old: 0,
            postfilter_gain: 0.0,
            postfilter_gain_old: 0.0,
            postfilter_tapset: 0,
            postfilter_tapset_old: 0,
            preemph_mem: [0.0; 2],
            rng: 0,
            loss_count: 0,
            last_pitch: 0,
            disable_inv: channels == 1,
            synth: Synth::new(),
        };
        d.reset();
        d
    }

    /// Returns the decoder to its initial state (§4.5.2).
    pub fn reset(&mut self) {
        for m in &mut self.decode_mem {
            m.fill(0.0);
        }
        self.old_band_e = [0.0; 2 * NB_EBANDS];
        self.old_log_e = [-28.0; 2 * NB_EBANDS];
        self.old_log_e2 = [-28.0; 2 * NB_EBANDS];
        self.postfilter_period = 0;
        self.postfilter_period_old = 0;
        self.postfilter_gain = 0.0;
        self.postfilter_gain_old = 0.0;
        self.postfilter_tapset = 0;
        self.postfilter_tapset_old = 0;
        self.preemph_mem = [0.0; 2];
        self.rng = 0;
        self.loss_count = 0;
        self.last_pitch = 0;
    }

    /// Output channels.
    pub fn channels(&self) -> usize {
        self.channels
    }

    /// Decodes one frame of `n` samples (at 48 kHz) coded with `c` channels
    /// from `ec`, whose frame is `ec.storage()` bytes, coding bands
    /// `start..end`. Writes `n / downsample` interleaved samples per channel
    /// to `out` (full scale ±1.0); `accumulate` adds instead of writing.
    #[allow(clippy::too_many_arguments)]
    pub fn decode(
        &mut self,
        ec: &mut RangeDecoder,
        n: usize,
        c: usize,
        start: usize,
        end: usize,
        out: &mut [f32],
        accumulate: bool,
    ) -> FrameInfo {
        let len = ec.storage();
        if len <= 1 {
            self.decode_lost(n, out, accumulate);
            return FrameInfo::default();
        }
        let lm = match n {
            120 => 0,
            240 => 1,
            480 => 2,
            _ => 3,
        };
        let mm = 1usize << lm;
        let cc = self.channels;
        if c == 1 {
            for i in 0..NB_EBANDS {
                self.old_band_e[i] = self.old_band_e[i].max(self.old_band_e[NB_EBANDS + i]);
            }
        }
        let total_bits = (len * 8) as i32;
        let mut tell = ec.tell();
        let silence = if tell >= total_bits {
            true
        } else if tell == 1 {
            ec.bit_logp(15)
        } else {
            false
        };
        if silence {
            // Pretend every bit was used.
            ec.add_bits(total_bits - ec.tell());
            tell = total_bits;
        }
        let mut pf_gain = 0.0f32;
        let mut pf_pitch = 0usize;
        let mut pf_tapset = 0usize;
        if start == 0 && tell + 16 <= total_bits {
            if ec.bit_logp(1) {
                let octave = ec.uint(6);
                pf_pitch = ((16usize << octave) + ec.bits(4 + octave) as usize) - 1;
                let qg = ec.bits(3);
                if ec.tell() + 2 <= total_bits {
                    pf_tapset = ec.icdf(&TAPSET_ICDF, 2);
                }
                pf_gain = 0.09375 * (qg + 1) as f32;
            }
            tell = ec.tell();
        }
        let is_transient = if lm > 0 && tell + 3 <= total_bits {
            let t = ec.bit_logp(3);
            tell = ec.tell();
            t
        } else {
            false
        };
        let short_blocks = if is_transient { mm } else { 0 };
        let intra = if tell + 3 <= total_bits { ec.bit_logp(3) } else { false };
        energy::unquant_coarse(ec, &mut self.old_band_e, start, end, intra, c, lm);
        let tf_res = tf_decode(ec, start, end, is_transient, lm, total_bits);
        tell = ec.tell();
        let mut spread = SPREAD_NORMAL;
        if tell + 4 <= total_bits {
            spread = ec.icdf(&SPREAD_ICDF, 5) as u32;
        }
        let cap = init_caps(lm, c);
        let mut offsets = [0i32; NB_EBANDS];
        let mut dynalloc_logp = 6;
        let mut total_bits8 = total_bits << BITRES;
        let mut tell8 = ec.tell_frac();
        for i in start..end {
            let width = (c * (EBANDS[i + 1] - EBANDS[i]) << lm) as i32;
            let quanta = (width << BITRES).min((6 << BITRES).max(width));
            let mut loop_logp = dynalloc_logp;
            let mut boost = 0;
            while tell8 + ((loop_logp as i32) << BITRES) < total_bits8 && boost < cap[i] {
                let flag = ec.bit_logp(loop_logp);
                tell8 = ec.tell_frac();
                if !flag {
                    break;
                }
                boost += quanta;
                total_bits8 -= quanta;
                loop_logp = 1;
            }
            offsets[i] = boost;
            if boost > 0 {
                dynalloc_logp = 2.max(dynalloc_logp - 1);
            }
        }
        let mut alloc_trim = 5;
        if tell8 + (6 << BITRES) <= total_bits8 {
            alloc_trim = ec.icdf(&TRIM_ICDF, 7) as i32;
        }
        let mut bits = ((len * 8) as i32) << BITRES;
        bits -= ec.tell_frac() + 1;
        let anti_collapse_rsv = if is_transient && lm >= 2 && bits >= (lm as i32 + 2) << BITRES { 1 << BITRES } else { 0 };
        bits -= anti_collapse_rsv;
        let alloc =
            rate::compute_allocation(start, end, &offsets, &cap, alloc_trim, bits, c, lm, ec, EncoderChoices::default());
        let mut err = [0.0f32; 2 * NB_EBANDS];
        energy::code_fine(ec, &mut self.old_band_e, &mut err, &alloc.fine_quant, start, end, c);
        let size = n;
        let mut x = vec![0.0f32; size];
        let mut y = vec![0.0f32; if c == 2 { size } else { 0 }];
        let params = FrameParams {
            start,
            end,
            lm,
            short_blocks: if is_transient { mm } else { 1 },
            spread,
            dual_stereo: alloc.dual_stereo,
            intensity: alloc.intensity,
            tf_res: &tf_res,
            total_bits: (((len * 8) as i32) << BITRES) - anti_collapse_rsv,
            balance: alloc.balance,
            pulses: &alloc.pulses,
            coded_bands: alloc.coded_bands,
            disable_inv: self.disable_inv,
            resynth: true,
        };
        let (collapse, seed) = bands::quant_all_bands(ec, &params, &mut x, &mut y, &[], self.rng);
        let anti_collapse_on = anti_collapse_rsv > 0 && ec.bits(1) != 0;
        if std::env::var_os("OPUS_DEBUG").is_some() {
            eprintln!(
                "celt lm {lm} c {c} tr {} intra {intra} spread {spread} tf {:?} dual {} int {} coded {} ac {anti_collapse_on} pf {pf_gain} trim {alloc_trim} boosts {:?}",
                is_transient, &tf_res[start..end], alloc.dual_stereo, alloc.intensity, alloc.coded_bands, &offsets[start..end]
            );
        }
        let left = (len * 8) as i32 - ec.tell();
        energy::code_finalise(
            ec,
            &mut self.old_band_e,
            &mut err,
            &alloc.fine_quant,
            &alloc.fine_priority,
            left,
            start,
            end,
            c,
        );
        let mut xy = x;
        xy.extend_from_slice(&y);
        if anti_collapse_on {
            bands::anti_collapse(
                &mut xy,
                size,
                &collapse,
                lm,
                c,
                start,
                end,
                &self.old_band_e,
                &self.old_log_e,
                &self.old_log_e2,
                &alloc.pulses,
                seed,
            );
        }
        if silence {
            self.old_band_e = [-28.0; 2 * NB_EBANDS];
        }
        self.synthesize(&xy, n, c, start, end, is_transient, silence);
        self.postfilter_and_output(n, lm, pf_pitch, pf_gain, pf_tapset, out, accumulate);
        // Energy history.
        if c == 1 {
            let (a, b) = self.old_band_e.split_at_mut(NB_EBANDS);
            b.copy_from_slice(a);
        }
        if !is_transient {
            self.old_log_e2 = self.old_log_e;
            self.old_log_e = self.old_band_e;
        } else {
            for i in 0..2 * NB_EBANDS {
                self.old_log_e[i] = self.old_log_e[i].min(self.old_band_e[i]);
            }
        }
        for ch in 0..2 {
            for i in (0..start).chain(end..NB_EBANDS) {
                self.old_band_e[ch * NB_EBANDS + i] = 0.0;
                self.old_log_e[ch * NB_EBANDS + i] = -28.0;
                self.old_log_e2[ch * NB_EBANDS + i] = -28.0;
            }
        }
        self.rng = ec.range();
        self.loss_count = 0;
        FrameInfo { pitch: if pf_gain > 0.0 { Some(pf_pitch as u32) } else { None } }
    }

    /// Denormalizes the bands and runs the inverse MDCTs, overlap-adding into
    /// the decode buffers (which are shifted by `n` first).
    #[allow(clippy::too_many_arguments)]
    fn synthesize(&mut self, xy: &[f32], n: usize, c: usize, start: usize, end: usize, transient: bool, silence: bool) {
        let cc = self.channels;
        let lm = match n {
            120 => 0,
            240 => 1,
            480 => 2,
            _ => 3,
        };
        let mm = 1 << lm;
        let mut bound = mm * EBANDS[end];
        if self.downsample != 1 {
            bound = bound.min(n / self.downsample);
        }
        let denorm = |x: &[f32], e: &[f32]| -> Vec<f32> {
            let mut f = vec![0.0f32; n];
            if silence {
                return f;
            }
            for i in start..end {
                let lg = (e[i] + E_MEANS[i]).min(32.0);
                let g = lg.exp2();
                for j in mm * EBANDS[i]..(mm * EBANDS[i + 1]).min(bound) {
                    f[j] = x[j] * g;
                }
            }
            f
        };
        let freqs: Vec<Vec<f32>> = if cc == 2 && c == 1 {
            let f = denorm(&xy[..n], &self.old_band_e[..NB_EBANDS]);
            vec![f.clone(), f]
        } else if cc == 1 && c == 2 {
            let f1 = denorm(&xy[..n], &self.old_band_e[..NB_EBANDS]);
            let f2 = denorm(&xy[n..2 * n], &self.old_band_e[NB_EBANDS..]);
            vec![f1.iter().zip(&f2).map(|(a, b)| 0.5 * a + 0.5 * b).collect()]
        } else {
            (0..cc).map(|ch| denorm(&xy[ch * n..(ch + 1) * n], &self.old_band_e[ch * NB_EBANDS..])).collect()
        };
        let blocks = if transient { mm } else { 1 };
        for (ch, f) in freqs.iter().enumerate() {
            let mem = &mut self.decode_mem[ch];
            mem.copy_within(n.., 0);
            let out = &mut mem[DECODE_BUFFER - n..];
            self.synth.imdct_ola(f, out, n, blocks);
        }
    }

    /// Applies the post-filter to this frame's `n` samples, de-emphasizes
    /// and writes them out.
    #[allow(clippy::too_many_arguments)]
    fn postfilter_and_output(
        &mut self,
        n: usize,
        lm: usize,
        pitch: usize,
        gain: f32,
        tapset: usize,
        out: &mut [f32],
        accumulate: bool,
    ) {
        let window = &mode().window;
        self.postfilter_period = self.postfilter_period.max(COMB_MIN_PERIOD);
        self.postfilter_period_old = self.postfilter_period_old.max(COMB_MIN_PERIOD);
        for ch in 0..self.channels {
            let mem = &mut self.decode_mem[ch];
            let base = DECODE_BUFFER - n;
            comb_filter(
                mem,
                base,
                self.postfilter_period_old,
                self.postfilter_period,
                SHORT_MDCT,
                self.postfilter_gain_old,
                self.postfilter_gain,
                self.postfilter_tapset_old,
                self.postfilter_tapset,
                window,
            );
            if lm != 0 {
                comb_filter(
                    mem,
                    base + SHORT_MDCT,
                    self.postfilter_period,
                    pitch,
                    n - SHORT_MDCT,
                    self.postfilter_gain,
                    gain,
                    self.postfilter_tapset,
                    tapset,
                    window,
                );
            }
        }
        self.postfilter_period_old = self.postfilter_period;
        self.postfilter_gain_old = self.postfilter_gain;
        self.postfilter_tapset_old = self.postfilter_tapset;
        self.postfilter_period = pitch;
        self.postfilter_gain = gain;
        self.postfilter_tapset = tapset;
        if lm != 0 {
            self.postfilter_period_old = self.postfilter_period;
            self.postfilter_gain_old = self.postfilter_gain;
            self.postfilter_tapset_old = self.postfilter_tapset;
        }
        self.deemphasis(n, out, accumulate);
    }

    fn deemphasis(&mut self, n: usize, out: &mut [f32], accumulate: bool) {
        let cc = self.channels;
        let ds = self.downsample;
        for ch in 0..cc {
            let mem = &self.decode_mem[ch];
            let mut m = self.preemph_mem[ch];
            let src = &mem[DECODE_BUFFER - n..DECODE_BUFFER];
            for (j, &s) in src.iter().enumerate() {
                let v = s + DEEMPHASIS * m;
                m = v;
                if j % ds == 0 {
                    let o = &mut out[(j / ds) * cc + ch];
                    let sample = v * (1.0 / 32768.0);
                    if accumulate {
                        *o += sample;
                    } else {
                        *o = sample;
                    }
                }
            }
            self.preemph_mem[ch] = m;
        }
    }

    /// Packet loss concealment (§4.4): repeats the last pitch period of the
    /// output with a decaying gain, through the same overlap-add as decoded
    /// frames so the transition keeps time-domain aliasing cancellation.
    pub fn decode_lost(&mut self, n: usize, out: &mut [f32], accumulate: bool) {
        let cc = self.channels;
        if self.loss_count == 0 {
            self.last_pitch = self.find_pitch();
        }
        self.loss_count += 1;
        let t = self.last_pitch.max(COMB_MIN_PERIOD);
        let att0 = 0.8f32.powi(self.loss_count as i32 - 1);
        let att1 = 0.8f32.powi(self.loss_count as i32);
        let frame_att = if self.loss_count > 5 { 0.0 } else { 1.0 };
        let l = OVERLAP;
        for ch in 0..cc {
            let mem = &self.decode_mem[ch];
            // Periodic extension of the final output.
            let mut e = vec![0.0f32; n + l];
            for (k, v) in e.iter_mut().enumerate() {
                let src = DECODE_BUFFER - t + (k % t);
                let g = att0 + (att1 - att0) * (k as f32 / (n + l) as f32);
                *v = mem[src.min(DECODE_BUFFER - 1)] * g * frame_att;
            }
            // Into the MDCT domain and back so the overlap stays consistent.
            let coefs = self.synth.mdct_windowed(&e, n);
            let mem = &mut self.decode_mem[ch];
            mem.copy_within(n.., 0);
            let o = &mut mem[DECODE_BUFFER - n..];
            self.synth.imdct_ola(&coefs, o, n, 1);
        }
        // The extension is already filtered: no post-filter this frame.
        self.postfilter_gain = 0.0;
        self.postfilter_gain_old = 0.0;
        for e in self.old_band_e.iter_mut() {
            *e = (*e - 0.5).max(-28.0);
        }
        self.deemphasis(n, out, accumulate);
    }

    /// A pitch period of the recent output, by normalized autocorrelation.
    fn find_pitch(&self) -> usize {
        let w = 480;
        let mono: Vec<f32> = (0..DECODE_BUFFER)
            .map(|i| self.decode_mem.iter().map(|m| m[i]).sum::<f32>() / self.channels as f32)
            .collect();
        let tail = &mono[DECODE_BUFFER - w..];
        let mut best = (0.0f32, 240usize);
        for t in 30..=720 {
            let past = &mono[DECODE_BUFFER - w - t..DECODE_BUFFER - t];
            let xy: f32 = tail.iter().zip(past).map(|(a, b)| a * b).sum();
            let yy: f32 = past.iter().map(|b| b * b).sum::<f32>() + 1e-9;
            let score = if xy > 0.0 { xy * xy / yy } else { 0.0 };
            if score > best.0 {
                best = (score, t);
            }
        }
        best.1
    }
}

const COMB_MIN_PERIOD: usize = 15;

/// §4.3.7.1: the comb post-filter, in place on `buf[base .. base + n]`,
/// cross-fading from `(t0, g0, tapset0)` to `(t1, g1, tapset1)` over the
/// overlap.
#[allow(clippy::too_many_arguments)]
pub(crate) fn comb_filter(
    buf: &mut [f32],
    base: usize,
    t0: usize,
    t1: usize,
    n: usize,
    g0: f32,
    g1: f32,
    tapset0: usize,
    tapset1: usize,
    window: &[f32],
) {
    if g0 == 0.0 && g1 == 0.0 {
        return;
    }
    let t0 = t0.max(COMB_MIN_PERIOD);
    let t1 = t1.max(COMB_MIN_PERIOD);
    let g00 = g0 * POSTFILTER_TAPS[tapset0][0];
    let g01 = g0 * POSTFILTER_TAPS[tapset0][1];
    let g02 = g0 * POSTFILTER_TAPS[tapset0][2];
    let g10 = g1 * POSTFILTER_TAPS[tapset1][0];
    let g11 = g1 * POSTFILTER_TAPS[tapset1][1];
    let g12 = g1 * POSTFILTER_TAPS[tapset1][2];
    let mut overlap = OVERLAP.min(n);
    if g0 == g1 && t0 == t1 && tapset0 == tapset1 {
        overlap = 0;
    }
    let mut x1 = buf[base - t1 + 1];
    let mut x2 = buf[base - t1];
    let mut x3 = buf[base - t1 - 1];
    let mut x4 = buf[base - t1 - 2];
    for i in 0..overlap {
        let p = base + i;
        let x0 = buf[p - t1 + 2];
        let f = window[i] * window[i];
        let v = buf[p]
            + (1.0 - f) * g00 * buf[p - t0]
            + (1.0 - f) * g01 * (buf[p - t0 + 1] + buf[p - t0 - 1])
            + (1.0 - f) * g02 * (buf[p - t0 + 2] + buf[p - t0 - 2])
            + f * g10 * x2
            + f * g11 * (x1 + x3)
            + f * g12 * (x0 + x4);
        buf[p] = v;
        x4 = x3;
        x3 = x2;
        x2 = x1;
        x1 = x0;
    }
    if g1 == 0.0 {
        return;
    }
    for i in overlap..n {
        let p = base + i;
        let v = buf[p]
            + g10 * buf[p - t1]
            + g11 * (buf[p - t1 + 1] + buf[p - t1 - 1])
            + g12 * (buf[p - t1 + 2] + buf[p - t1 - 2]);
        buf[p] = v;
    }
}

/// §4.3.1 / §4.3.4.5: the per-band TF changes.
pub(crate) fn tf_decode(
    ec: &mut RangeDecoder,
    start: usize,
    end: usize,
    transient: bool,
    lm: usize,
    total_bits: i32,
) -> [i32; NB_EBANDS] {
    let mut tf_res = [0i32; NB_EBANDS];
    let mut budget = total_bits;
    let mut tell = ec.tell();
    let mut logp = if transient { 2 } else { 4 };
    let tf_select_rsv = lm > 0 && tell + logp + 1 <= budget;
    budget -= i32::from(tf_select_rsv);
    let mut tf_changed = 0;
    let mut curr = 0;
    for r in tf_res.iter_mut().take(end).skip(start) {
        if tell + logp <= budget {
            curr ^= i32::from(ec.bit_logp(logp as u32));
            tell = ec.tell();
            tf_changed |= curr;
        }
        *r = curr;
        logp = if transient { 4 } else { 5 };
    }
    let ti = 4 * usize::from(transient);
    let mut tf_select = 0;
    if tf_select_rsv && TF_SELECT[lm][ti + tf_changed as usize] != TF_SELECT[lm][ti + 2 + tf_changed as usize] {
        tf_select = usize::from(ec.bit_logp(1));
    }
    for r in tf_res.iter_mut().take(end).skip(start) {
        *r = i32::from(TF_SELECT[lm][ti + 2 * tf_select + *r as usize]);
    }
    tf_res
}
