//! Band shapes (CELT_SPEC §8; RFC 6716 §4.3.4, §5.3.8): the per-band
//! loop with its running balance, the recursive splitting with the theta
//! parameter, stereo (mid/side, intensity, dual), the time-frequency
//! changes, PVQ coding with spreading, and folding / noise filling.
//!
//! The procedures serve both sides through [`Coder`]: every budget, split
//! and theta resolution is computed identically. The decoder reconstructs
//! the normalised spectrum (and the folding data later bands read); the
//! encoder transforms its input into the domain each step codes in and
//! never resynthesises, which the bitstream does not require (CELT_SPEC
//! §8.1, §12).

use super::cwrs;
use super::energy::tell;
use super::mode::{get_pulses, mode};
use super::rate::{bits2pulses, pulses2bits};
use super::tables::{
    BITRES, EBANDS, NB_EBANDS, QTHETA_OFFSET, QTHETA_OFFSET_TWOPHASE, SPREAD_FACTOR,
    TF_SELECT_TABLE,
};
use crate::range::Coder;

/// Spreading choices (CELT_SPEC §1.2 step 7).
pub(crate) const SPREAD_NONE: u32 = 0;
pub(crate) const SPREAD_LIGHT: u32 = 1;
pub(crate) const SPREAD_NORMAL: u32 = 2;
pub(crate) const SPREAD_AGGRESSIVE: u32 = 3;

/// The folding / anti-collapse generator (CELT_SPEC §8.2.7).
pub(crate) fn lcg(seed: u32) -> u32 {
    seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223)
}

/// The time-frequency resolution changes (CELT_SPEC §4; RFC 6716 §4.3.1,
/// §4.3.4.5): one flag per band, coded as changes, and `tf_select`.
/// `flags` and `select` are the encoder's choices; a flag that cannot be
/// coded takes the running value. Returns the per-band tf_change.
pub(crate) fn code_tf<C: Coder>(
    ec: &mut C,
    start: usize,
    end: usize,
    transient: bool,
    lm: usize,
    total_bits: i32,
    flags: &[i32; NB_EBANDS],
    select: bool,
) -> [i32; NB_EBANDS] {
    let mut budget = total_bits;
    let mut t = tell(ec);
    let mut logp: i32 = if transient { 2 } else { 4 };
    // `t + logp + 1 <= budget`.
    let tf_select_rsv = lm > 0 && t + logp < budget;
    if tf_select_rsv {
        budget -= 1;
    }
    let mut curr = 0;
    let mut changed = 0;
    let mut res = [0i32; NB_EBANDS];
    for i in start..end {
        if t + logp <= budget {
            curr ^= i32::from(ec.bit_logp(flags[i] != curr, logp as u32));
            t = tell(ec);
            changed |= curr;
        }
        res[i] = curr;
        logp = if transient { 4 } else { 5 };
    }
    let row = &TF_SELECT_TABLE[lm];
    let ti = 4 * usize::from(transient);
    let changed = changed as usize;
    let sel = tf_select_rsv && row[ti + changed] != row[ti + 2 + changed] && ec.bit_logp(select, 1);
    for r in &mut res[start..end] {
        *r = i32::from(row[ti + 2 * usize::from(sel) + *r as usize]);
    }
    res
}

/// The frame-level inputs of the band loop (CELT_SPEC §8.1).
pub(crate) struct FrameBands<'a> {
    pub start: usize,
    pub end: usize,
    pub lm: usize,
    /// B: M for transient frames, else 1.
    pub blocks: usize,
    pub spread: u32,
    pub dual_stereo: bool,
    pub intensity: usize,
    pub tf_res: &'a [i32; NB_EBANDS],
    /// `8·8·len − anti_collapse_rsv`.
    pub shape_total: i32,
    pub balance: i32,
    pub pulses: &'a [i32; NB_EBANDS],
    pub coded_bands: usize,
    /// RFC 8251 §10: ignore the decoded phase inversion.
    pub disable_inv: bool,
}

/// Collapse masks per band and channel (CELT_SPEC §8.1 step 10).
pub(crate) type CollapseMasks = [[u8; 2]; NB_EBANDS];

/// The band loop of CELT_SPEC §8.1. `x` (and `y` for stereo) hold the
/// whole frame's spectrum: on the encoder the normalised input, on the
/// decoder the normalised output. `amps` are the encoder's band
/// amplitudes per channel (its intensity downmix weights); the decoder
/// ignores them. `seed` is the folding generator state.
pub(crate) fn quant_all_bands<C: Coder>(
    ec: &mut C,
    p: &FrameBands,
    x: &mut [f32],
    mut y: Option<&mut [f32]>,
    amps: &[[f32; NB_EBANDS]; 2],
    seed: &mut u32,
) -> CollapseMasks {
    let m = 1usize << p.lm;
    let lm = p.lm as i32;
    let stereo = y.is_some();
    let base = m * EBANDS[p.start];
    let norm_len = m * EBANDS[NB_EBANDS];
    let (mut norm, mut norm2) = if C::ENCODE {
        (Vec::new(), Vec::new())
    } else {
        (
            vec![0.0f32; norm_len],
            vec![0.0f32; if stereo { norm_len } else { 0 }],
        )
    };
    let mut collapse: CollapseMasks = [[0; 2]; NB_EBANDS];
    let mut balance = p.balance;
    let mut lowband_offset = 0usize;
    let mut update_lowband = true;
    let mut dual = p.dual_stereo;
    for i in p.start..p.end {
        let o = m * EBANDS[i];
        let n = m * EBANDS[i + 1] - o;
        let tell = ec.tell_frac();
        if i != p.start {
            balance -= tell;
        }
        let remaining_bits = p.shape_total - tell - 1;
        let b = if i < p.coded_bands {
            let curr_balance = balance / (p.coded_bands - i).min(3) as i32;
            (remaining_bits + 1)
                .min(p.pulses[i] + curr_balance)
                .clamp(0, 16383)
        } else {
            0
        };
        let tf_change = p.tf_res[i];
        // Folding source (decoder only; it never affects the bitstream).
        let mut lowband_at = None;
        let (mut x_cm, mut y_cm) = ((1u32 << p.blocks) - 1, (1u32 << p.blocks) - 1);
        if !C::ENCODE {
            if (o as i32 - n as i32 >= base as i32 || i == p.start + 1)
                && (update_lowband || lowband_offset == 0)
            {
                lowband_offset = i;
            }
            if i == p.start + 1 {
                // RFC 8251 §9: extend the first band's folding data to the
                // width of the second.
                let n1 = m * (EBANDS[p.start + 1] - EBANDS[p.start]);
                let n2 = n;
                if n2 > n1 {
                    norm.copy_within(base + 2 * n1 - n2..base + n1, base + n1);
                    if stereo {
                        norm2.copy_within(base + 2 * n1 - n2..base + n1, base + n1);
                    }
                }
            }
            if lowband_offset != 0
                && (p.spread != SPREAD_AGGRESSIVE || p.blocks > 1 || tf_change < 0)
            {
                let eff = base.max((m * EBANDS[lowband_offset]).saturating_sub(n));
                let mut fold_start = lowband_offset;
                loop {
                    fold_start -= 1;
                    if m * EBANDS[fold_start] <= eff {
                        break;
                    }
                }
                let mut fold_end = lowband_offset;
                while fold_end < i && m * EBANDS[fold_end] < eff + n {
                    fold_end += 1;
                }
                let (mut xc, mut yc) = (0u32, 0u32);
                for f in fold_start..fold_end.max(fold_start + 1) {
                    xc |= u32::from(collapse[f][0]);
                    yc |= u32::from(collapse[f][usize::from(stereo)]);
                }
                x_cm = xc;
                y_cm = yc;
                lowband_at = Some(eff);
            }
        }
        if dual && i == p.intensity {
            dual = false;
            if !C::ENCODE {
                for k in base..o {
                    norm[k] = 0.5 * (norm[k] + norm2[k]);
                }
            }
        }
        let mut ctx = Ctx {
            ec: &mut *ec,
            band: i,
            spread: p.spread,
            intensity: p.intensity,
            tf_change,
            remaining_bits,
            seed: *seed,
            disable_inv: p.disable_inv,
            amps: [amps[0][i], amps[1][i]],
        };
        let lowband_of = |buf: &Vec<f32>| lowband_at.map(|e: usize| buf[e..e + n].to_vec());
        let lowband_x = if C::ENCODE { None } else { lowband_of(&norm) };
        let xb = &mut x[o..o + n];
        let out_x = if C::ENCODE {
            None
        } else {
            Some(&mut norm[o..o + n])
        };
        match y.as_deref_mut() {
            Some(yfull) if dual => {
                let lowband_y = if C::ENCODE { None } else { lowband_of(&norm2) };
                x_cm = ctx.mono_band(
                    xb,
                    b >> 1,
                    p.blocks,
                    lowband_x.as_deref(),
                    lm,
                    out_x,
                    1.0,
                    x_cm,
                );
                let out_y = if C::ENCODE {
                    None
                } else {
                    Some(&mut norm2[o..o + n])
                };
                y_cm = ctx.mono_band(
                    &mut yfull[o..o + n],
                    b >> 1,
                    p.blocks,
                    lowband_y.as_deref(),
                    lm,
                    out_y,
                    1.0,
                    y_cm,
                );
            }
            Some(yfull) => {
                x_cm = ctx.stereo_band(
                    xb,
                    &mut yfull[o..o + n],
                    b,
                    p.blocks,
                    lowband_x.as_deref(),
                    lm,
                    out_x,
                    x_cm | y_cm,
                );
                y_cm = x_cm;
            }
            None => {
                x_cm = ctx.mono_band(
                    xb,
                    b,
                    p.blocks,
                    lowband_x.as_deref(),
                    lm,
                    out_x,
                    1.0,
                    x_cm | y_cm,
                );
                y_cm = x_cm;
            }
        }
        *seed = ctx.seed;
        collapse[i] = [x_cm as u8, y_cm as u8];
        balance += p.pulses[i] + tell;
        update_lowband = b > (n as i32) << BITRES;
    }
    collapse
}

/// The state shared by the nested calls coding one band (CELT_SPEC §8.2).
struct Ctx<'a, C: Coder> {
    ec: &'a mut C,
    band: usize,
    spread: u32,
    intensity: usize,
    tf_change: i32,
    /// The single counter of CELT_SPEC §8.1, decremented as bits are used.
    remaining_bits: i32,
    seed: u32,
    disable_inv: bool,
    /// Encoder: the band's amplitude in each channel.
    amps: [f32; 2],
}

/// The outcome of the theta step (CELT_SPEC §8.2.4).
struct Theta {
    itheta: i32,
    qalloc: i32,
    imid: i32,
    iside: i32,
    delta: i32,
    inv: bool,
    fill: u32,
}

impl<C: Coder> Ctx<'_, C> {
    /// A mono band at level 0 (CELT_SPEC §8.2.1, §8.2.2, §8.2.8): the TF
    /// changes around the recursive coding, and the folding output.
    fn mono_band(
        &mut self,
        x: &mut [f32],
        b: i32,
        blocks: usize,
        lowband: Option<&[f32]>,
        lm: i32,
        lowband_out: Option<&mut [f32]>,
        gain: f32,
        fill: u32,
    ) -> u32 {
        let n = x.len();
        if n == 1 {
            return self.single(x, None, lowband_out);
        }
        let m = mode();
        let long_blocks = blocks == 1;
        let mut tf_change = self.tf_change;
        let recombine = tf_change.max(0) as usize;
        let mut fill = fill;
        let mut blocks = blocks;
        let mut n_b = n / blocks;
        let mut lb: Option<Vec<f32>> = lowband.map(<[f32]>::to_vec);
        for k in 0..recombine {
            if C::ENCODE {
                haar(x, n >> k, 1 << k);
            }
            if let Some(l) = lb.as_deref_mut() {
                haar(l, n >> k, 1 << k);
            }
            fill = m.bit_interleave[(fill & 15) as usize]
                | m.bit_interleave[(fill >> 4) as usize] << 2;
        }
        blocks >>= recombine;
        n_b <<= recombine;
        let mut time_divide = 0;
        while n_b & 1 == 0 && tf_change < 0 {
            if C::ENCODE {
                haar(x, n_b, blocks);
            }
            if let Some(l) = lb.as_deref_mut() {
                haar(l, n_b, blocks);
            }
            fill |= fill << blocks;
            blocks <<= 1;
            n_b >>= 1;
            time_divide += 1;
            tf_change += 1;
        }
        let b0 = blocks;
        let n_b0 = n_b;
        if b0 > 1 {
            if C::ENCODE {
                deinterleave(x, n_b >> recombine, b0 << recombine, long_blocks);
            }
            if let Some(l) = lb.as_deref_mut() {
                deinterleave(l, n_b >> recombine, b0 << recombine, long_blocks);
            }
        }
        let mut cm = self.partition(x, b, b0, lb.as_deref(), lm, gain, fill);
        if C::ENCODE {
            return cm;
        }
        if b0 > 1 {
            interleave(x, n_b0 >> recombine, b0 << recombine, long_blocks);
        }
        let mut nb = n_b0;
        let mut bb = b0;
        for _ in 0..time_divide {
            bb >>= 1;
            nb <<= 1;
            cm |= cm >> bb;
            haar(x, nb, bb);
        }
        for k in 0..recombine {
            cm = m.bit_deinterleave[cm as usize];
            haar(x, n >> k, 1 << k);
        }
        bb <<= recombine;
        if let Some(out) = lowband_out {
            let s = (n as f32).sqrt();
            for (o, v) in out.iter_mut().zip(x.iter()) {
                *o = s * v;
            }
        }
        cm & ((1 << bb) - 1)
    }

    /// One-coefficient bands: just the signs (CELT_SPEC §8.2.1).
    fn single(
        &mut self,
        x: &mut [f32],
        y: Option<&mut [f32]>,
        lowband_out: Option<&mut [f32]>,
    ) -> u32 {
        for v in std::iter::once(&mut *x).chain(y) {
            let mut sign = 0;
            if self.remaining_bits >= 1 << BITRES {
                sign = self.ec.bits(u32::from(v[0] < 0.0), 1);
                self.remaining_bits -= 1 << BITRES;
            }
            if !C::ENCODE {
                v[0] = if sign != 0 { -1.0 } else { 1.0 };
            }
        }
        if let Some(out) = lowband_out {
            out[0] = x[0];
        }
        1
    }

    /// A mono partition at any level: split in halves while the budget
    /// exceeds the largest codebook, else PVQ (CELT_SPEC §8.2.3, §8.2.6,
    /// §8.2.7). `blocks` on entry is this call's B0.
    fn partition(
        &mut self,
        x: &mut [f32],
        b: i32,
        blocks: usize,
        lowband: Option<&[f32]>,
        lm: i32,
        gain: f32,
        fill: u32,
    ) -> u32 {
        let n = x.len();
        let b0 = blocks;
        let cache = mode().cache(self.band, lm);
        if lm != -1 && b > i32::from(cache[usize::from(cache[0])]) + 12 && n > 2 {
            let half = n / 2;
            let lm = lm - 1;
            let mut fill = fill;
            if blocks == 1 {
                fill = (fill & 1) | (fill << 1);
            }
            let blocks = (blocks + 1) >> 1;
            let (xm, xs) = x.split_at_mut(half);
            let th = self.theta(xm, xs, b, b0, blocks, lm, false, fill);
            let b = b - th.qalloc;
            let mut delta = th.delta;
            if b0 > 1 && th.itheta & 0x3fff != 0 {
                if th.itheta > 8192 {
                    delta -= delta >> (4 - lm);
                } else {
                    delta = (delta + (((half as i32) << BITRES) >> (5 - lm))).min(0);
                }
            }
            let mut mbits = ((b - delta) / 2).min(b).max(0);
            let mut sbits = b - mbits;
            self.remaining_bits -= th.qalloc;
            let mid = th.imid as f32 / 32768.0;
            let side = th.iside as f32 / 32768.0;
            let (lb_mid, lb_side) = match lowband {
                Some(l) => (Some(&l[..half]), Some(&l[half..])),
                None => (None, None),
            };
            let r0 = self.remaining_bits;
            let (cm_mid, cm_side);
            if mbits >= sbits {
                cm_mid = self.partition(xm, mbits, blocks, lb_mid, lm, gain * mid, th.fill);
                let rebalance = mbits - (r0 - self.remaining_bits);
                if rebalance > 3 << BITRES && th.itheta != 0 {
                    sbits += rebalance - (3 << BITRES);
                }
                cm_side = self.partition(
                    xs,
                    sbits,
                    blocks,
                    lb_side,
                    lm,
                    gain * side,
                    th.fill >> blocks,
                );
            } else {
                cm_side = self.partition(
                    xs,
                    sbits,
                    blocks,
                    lb_side,
                    lm,
                    gain * side,
                    th.fill >> blocks,
                );
                let rebalance = sbits - (r0 - self.remaining_bits);
                if rebalance > 3 << BITRES && th.itheta != 16384 {
                    mbits += rebalance - (3 << BITRES);
                }
                cm_mid = self.partition(xm, mbits, blocks, lb_mid, lm, gain * mid, th.fill);
            }
            return cm_mid | cm_side << (b0 >> 1);
        }
        let mut q = bits2pulses(self.band, lm, b);
        let mut cost = pulses2bits(self.band, lm, q);
        self.remaining_bits -= cost;
        while self.remaining_bits < 0 && q > 0 {
            self.remaining_bits += cost;
            q -= 1;
            cost = pulses2bits(self.band, lm, q);
            self.remaining_bits -= cost;
        }
        if q > 0 {
            return self.pvq(x, get_pulses(q), blocks, gain);
        }
        if C::ENCODE {
            return 0;
        }
        let all = (1u32 << blocks) - 1;
        let fill = fill & all;
        if fill == 0 {
            x.fill(0.0);
            return 0;
        }
        match lowband {
            None => {
                for v in x.iter_mut() {
                    self.seed = lcg(self.seed);
                    *v = (self.seed as i32 >> 20) as f32;
                }
                renormalise(x, gain);
                all
            }
            Some(l) => {
                for (v, &s) in x.iter_mut().zip(l) {
                    self.seed = lcg(self.seed);
                    *v = s + if self.seed & 0x8000 != 0 {
                        1.0 / 256.0
                    } else {
                        -1.0 / 256.0
                    };
                }
                renormalise(x, gain);
                fill
            }
        }
    }

    /// A stereo band (CELT_SPEC §8.2.4, §8.2.5, §8.2.6, §8.2.8).
    fn stereo_band(
        &mut self,
        x: &mut [f32],
        y: &mut [f32],
        b: i32,
        blocks: usize,
        lowband: Option<&[f32]>,
        lm: i32,
        lowband_out: Option<&mut [f32]>,
        fill: u32,
    ) -> u32 {
        let n = x.len();
        if n == 1 {
            return self.single(x, Some(y), lowband_out);
        }
        let orig_fill = fill;
        let th = self.theta(x, y, b, blocks, blocks, lm, true, fill);
        let b = b - th.qalloc;
        let mid = th.imid as f32 / 32768.0;
        let side = th.iside as f32 / 32768.0;
        let cm;
        if n == 2 {
            let sbits = if th.itheta != 0 && th.itheta != 16384 {
                1 << BITRES
            } else {
                0
            };
            let mbits = b - sbits;
            self.remaining_bits -= th.qalloc + sbits;
            let swap = th.itheta > 8192;
            let (x2, y2) = if swap {
                (&mut *y, &mut *x)
            } else {
                (&mut *x, &mut *y)
            };
            let mut sign = 0;
            if sbits > 0 {
                sign = self
                    .ec
                    .bits(u32::from(x2[0] * y2[1] - x2[1] * y2[0] < 0.0), 1);
            }
            cm = self.mono_band(x2, mbits, blocks, lowband, lm, lowband_out, 1.0, orig_fill);
            if C::ENCODE {
                return cm;
            }
            let s = 1.0 - 2.0 * sign as f32;
            y2[0] = -s * x2[1];
            y2[1] = s * x2[0];
            for k in 0..2 {
                let (a, c) = (mid * x[k], side * y[k]);
                x[k] = a - c;
                y[k] = a + c;
            }
        } else {
            let mut mbits = ((b - th.delta) / 2).min(b).max(0);
            let mut sbits = b - mbits;
            self.remaining_bits -= th.qalloc;
            let r0 = self.remaining_bits;
            if mbits >= sbits {
                let cm_mid =
                    self.mono_band(x, mbits, blocks, lowband, lm, lowband_out, 1.0, th.fill);
                let rebalance = mbits - (r0 - self.remaining_bits);
                if rebalance > 3 << BITRES && th.itheta != 0 {
                    sbits += rebalance - (3 << BITRES);
                }
                cm = cm_mid
                    | self.mono_band(y, sbits, blocks, None, lm, None, side, th.fill >> blocks);
            } else {
                let cm_side =
                    self.mono_band(y, sbits, blocks, None, lm, None, side, th.fill >> blocks);
                let rebalance = sbits - (r0 - self.remaining_bits);
                if rebalance > 3 << BITRES && th.itheta != 16384 {
                    mbits += rebalance - (3 << BITRES);
                }
                cm = cm_side
                    | self.mono_band(x, mbits, blocks, lowband, lm, lowband_out, 1.0, th.fill);
            }
            if C::ENCODE {
                return cm;
            }
            stereo_merge(x, y, mid);
        }
        if th.inv {
            for v in y.iter_mut() {
                *v = -*v;
            }
        }
        cm
    }

    /// The theta step (CELT_SPEC §8.2.4): the split of the budget between
    /// the two halves (mono split) or mid and side (stereo). `b0` selects
    /// the PDF, `blocks` (the current B) the fill masks. On the encoder,
    /// a stereo `x`, `y` (left, right) are replaced by mid and side, or by
    /// the intensity downmix when no angle is coded.
    fn theta(
        &mut self,
        x: &mut [f32],
        y: &mut [f32],
        b: i32,
        b0: usize,
        blocks: usize,
        lm: i32,
        stereo: bool,
        fill: u32,
    ) -> Theta {
        let n = x.len() as i32;
        let pulse_cap = mode().log_n[self.band] + (lm << BITRES);
        let offset = (pulse_cap >> 1)
            - if stereo && n == 2 {
                QTHETA_OFFSET_TWOPHASE
            } else {
                QTHETA_OFFSET
            };
        let mut qn = compute_qn(n, b, offset, pulse_cap, stereo);
        if stereo && self.band >= self.intensity {
            qn = 1;
        }
        let mut target = 0i32;
        if C::ENCODE && qn != 1 {
            if stereo {
                for (l, r) in x.iter_mut().zip(y.iter_mut()) {
                    let (m, s) = (*l + *r, *r - *l);
                    *l = m;
                    *r = s;
                }
            }
            let ex: f32 = x.iter().map(|v| v * v).sum();
            let ey: f32 = y.iter().map(|v| v * v).sum();
            let angle = ey.sqrt().atan2(ex.sqrt());
            target =
                ((angle / std::f32::consts::FRAC_PI_2 * 16384.0).round() as i32).clamp(0, 16384);
        }
        let t0 = self.ec.tell_frac();
        let mut itheta = 0;
        let mut inv = false;
        if qn != 1 {
            let q = (target * qn + 8192) >> 14;
            itheta = if stereo && n > 2 {
                // Step PDF: weight 3 up to qn/2, 1 above.
                let x0 = qn / 2;
                let ft = 3 * (x0 + 1) + x0;
                let v = if C::ENCODE {
                    q
                } else {
                    let fs = self.ec.decode_fs(ft as u32) as i32;
                    if fs < (x0 + 1) * 3 {
                        fs / 3
                    } else {
                        x0 + 1 + (fs - (x0 + 1) * 3)
                    }
                };
                let (fl, fh) = if v <= x0 {
                    (3 * v, 3 * (v + 1))
                } else {
                    (v - 1 - x0 + (x0 + 1) * 3, v - x0 + (x0 + 1) * 3)
                };
                self.ec.code(fl as u32, fh as u32, ft as u32);
                v
            } else if b0 > 1 || stereo {
                self.ec.uint(q as u32, (qn + 1) as u32) as i32
            } else {
                // Triangular PDF.
                let h = qn >> 1;
                let ft = (h + 1) * (h + 1);
                let v = if C::ENCODE {
                    q
                } else {
                    let fm = self.ec.decode_fs(ft as u32) as i32;
                    if fm < (h * (h + 1)) >> 1 {
                        (isqrt(8 * fm as u32 + 1) as i32 - 1) >> 1
                    } else {
                        (2 * (qn + 1) - isqrt(8 * (ft - fm - 1) as u32 + 1) as i32) >> 1
                    }
                };
                let (fl, fs) = if v <= h {
                    ((v * (v + 1)) >> 1, v + 1)
                } else {
                    (ft - (((qn + 1 - v) * (qn + 2 - v)) >> 1), qn + 1 - v)
                };
                self.ec.code(fl as u32, (fl + fs) as u32, ft as u32);
                v
            };
            itheta = itheta * 16384 / qn;
        } else if stereo {
            if b > 2 << BITRES && self.remaining_bits > 2 << BITRES {
                let want =
                    C::ENCODE && x.iter().zip(y.iter()).map(|(l, r)| l * r).sum::<f32>() < 0.0;
                inv = self.ec.bit_logp(want, 2);
            }
            if C::ENCODE {
                // The intensity downmix, weighted by the channel amplitudes.
                let gr = if inv { -self.amps[1] } else { self.amps[1] };
                for (l, r) in x.iter_mut().zip(y.iter()) {
                    *l = self.amps[0] * *l + gr * r;
                }
            }
            inv &= !self.disable_inv;
        }
        let qalloc = self.ec.tell_frac() - t0;
        let all = (1u32 << blocks) - 1;
        let (imid, iside, delta, fill) = match itheta {
            0 => (32767, 0, -16384, fill & all),
            16384 => (0, 32767, 16384, fill & (all << blocks)),
            _ => {
                let imid = cosx(itheta);
                let iside = cosx(16384 - itheta);
                (
                    imid,
                    iside,
                    frac_mul16((n - 1) << 7, log2tan(iside, imid)),
                    fill,
                )
            }
        };
        Theta {
            itheta,
            qalloc,
            imid,
            iside,
            delta,
            inv,
            fill,
        }
    }

    /// PVQ with `k` pulses (CELT_SPEC §8.3): the encoder searches and codes
    /// the vector; the decoder decodes, normalises to `gain` and spreads.
    fn pvq(&mut self, x: &mut [f32], k: usize, blocks: usize, gain: f32) -> u32 {
        let n = x.len();
        let ft = cwrs::v(n, k) as u32;
        let mut yv_buf = [0i32; 176];
        let mut yv_heap = Vec::new();
        let yv: &mut [i32] = if n <= yv_buf.len() {
            &mut yv_buf[..n]
        } else {
            yv_heap.resize(n, 0);
            &mut yv_heap
        };
        if C::ENCODE {
            exp_rotation(x, blocks, k, self.spread, true);
            pvq_search(x, k, yv);
            self.ec.uint(cwrs::encode(yv, k), ft);
            return 0;
        }
        let idx = self.ec.uint(0, ft);
        cwrs::decode(idx, n, k, yv);
        let ryy: i64 = yv.iter().map(|&v| i64::from(v) * i64::from(v)).sum();
        let g = gain / (ryy as f32).sqrt();
        for (v, &q) in x.iter_mut().zip(yv.iter()) {
            *v = g * q as f32;
        }
        exp_rotation(x, blocks, k, self.spread, false);
        if blocks <= 1 {
            return 1;
        }
        let n0 = n / blocks;
        (0..blocks)
            .filter(|&bk| yv[bk * n0..(bk + 1) * n0].iter().any(|&v| v != 0))
            .fold(0, |cm, bk| cm | 1 << bk)
    }
}

/// The theta resolution (CELT_SPEC §8.2.4 step 2).
fn compute_qn(n: i32, b: i32, offset: i32, pulse_cap: i32, stereo: bool) -> i32 {
    let mut n2 = 2 * n - 1;
    if stereo && n == 2 {
        n2 -= 1;
    }
    let qb = (b - pulse_cap - (4 << BITRES))
        .min((b + n2 * offset) / n2)
        .min(8 << BITRES);
    if qb < (1 << BITRES >> 1) {
        return 1;
    }
    let qn = mode().exp2_table8[(qb & 7) as usize] >> (14 - (qb >> 3));
    ((qn + 1) >> 1) << 1
}

/// `(16384 + a·c) >> 15` (CELT_SPEC §8.2.4).
fn frac_mul16(a: i32, c: i32) -> i32 {
    (16384 + a * c) >> 15
}

/// The bit-exact cosine of CELT_SPEC §8.2.4, `x` in 1/16384 of π/2.
fn cosx(x: i32) -> i32 {
    let t = (4096 + x * x) >> 13;
    let r = (32767 - t) + frac_mul16(t, -7651 + frac_mul16(t, 8277 + frac_mul16(-626, t)));
    1 + r
}

/// The bit-exact log2(s/c) of CELT_SPEC §8.2.4, in 1/2048.
fn log2tan(s: i32, c: i32) -> i32 {
    let lc = 32 - (c as u32).leading_zeros() as i32;
    let ls = 32 - (s as u32).leading_zeros() as i32;
    let c = c << (15 - lc);
    let s = s << (15 - ls);
    (ls - lc) * 2048 + frac_mul16(s, frac_mul16(s, -2597) + 7932)
        - frac_mul16(c, frac_mul16(c, -2597) + 7932)
}

/// ⌊√x⌋, exact.
fn isqrt(x: u32) -> u32 {
    let x = u64::from(x);
    let mut r = (x as f64).sqrt() as u64;
    while r * r > x {
        r -= 1;
    }
    while (r + 1) * (r + 1) <= x {
        r += 1;
    }
    r as u32
}

/// The Haar step of CELT_SPEC §8.2.2 on `n` rows of `stride` interleaved
/// columns.
fn haar(v: &mut [f32], n: usize, stride: usize) {
    // 0.70710678 (CELT_SPEC §8.2.2).
    const S: f32 = std::f32::consts::FRAC_1_SQRT_2;
    for i in 0..stride {
        for j in 0..n / 2 {
            let (p, q) = (stride * 2 * j + i, stride * (2 * j + 1) + i);
            let a = S * v[p];
            let c = S * v[q];
            v[p] = a + c;
            v[q] = a - c;
        }
    }
}

/// The block order of `deinterleave` (CELT_SPEC §8.2.2).
fn block_order(stride: usize, hadamard: bool) -> impl Fn(usize) -> usize {
    let m = mode();
    move |i| {
        if hadamard {
            m.ordery[stride - 2 + i]
        } else {
            i
        }
    }
}

/// Coefficients interleaved by block become contiguous per block
/// (CELT_SPEC §8.2.2).
fn deinterleave(v: &mut [f32], n0: usize, stride: usize, hadamard: bool) {
    let ord = block_order(stride, hadamard);
    let mut t = vec![0.0f32; n0 * stride];
    for i in 0..stride {
        for j in 0..n0 {
            t[ord(i) * n0 + j] = v[j * stride + i];
        }
    }
    v[..n0 * stride].copy_from_slice(&t);
}

/// The inverse of [`deinterleave`].
fn interleave(v: &mut [f32], n0: usize, stride: usize, hadamard: bool) {
    let ord = block_order(stride, hadamard);
    let mut t = vec![0.0f32; n0 * stride];
    for i in 0..stride {
        for j in 0..n0 {
            t[j * stride + i] = v[ord(i) * n0 + j];
        }
    }
    v[..n0 * stride].copy_from_slice(&t);
}

/// Scales `x` to norm `gain` (CELT_SPEC §8.3.3).
pub(crate) fn renormalise(x: &mut [f32], gain: f32) {
    let e = 1e-15 + x.iter().map(|v| v * v).sum::<f32>();
    let g = gain / e.sqrt();
    for v in x {
        *v *= g;
    }
}

/// The stereo merge of CELT_SPEC §8.5: unit-norm left and right from the
/// unit-norm mid `x` and the scaled side `y`.
fn stereo_merge(x: &mut [f32], y: &mut [f32], mid: f32) {
    let xp = mid * x.iter().zip(y.iter()).map(|(a, b)| a * b).sum::<f32>();
    let sy: f32 = y.iter().map(|v| v * v).sum();
    let el = mid * mid + sy - 2.0 * xp;
    let er = mid * mid + sy + 2.0 * xp;
    if er < 6e-4 || el < 6e-4 {
        y.copy_from_slice(x);
        return;
    }
    let lg = 1.0 / el.sqrt();
    let rg = 1.0 / er.sqrt();
    for (a, b) in x.iter_mut().zip(y.iter_mut()) {
        let l = mid * *a;
        let r = *b;
        *a = lg * (l - r);
        *b = rg * (l + r);
    }
}

/// The spreading rotation of CELT_SPEC §8.3.4: the decoder's direction,
/// or its exact inverse for the encoder.
fn exp_rotation(x: &mut [f32], blocks: usize, k: usize, spread: u32, inverse: bool) {
    let n = x.len();
    if 2 * k >= n || spread == SPREAD_NONE {
        return;
    }
    let f = SPREAD_FACTOR[spread as usize - 1] as f32;
    let g = n as f32 / (n as f32 + f * k as f32);
    let theta = 0.5 * g * g;
    let c = (std::f32::consts::FRAC_PI_2 * theta).cos();
    let s = (std::f32::consts::FRAC_PI_2 * (1.0 - theta)).cos();
    let mut stride2 = 0;
    if n >= 8 * blocks {
        stride2 = 1;
        while (stride2 * stride2 + stride2) * blocks + (blocks >> 2) < n {
            stride2 += 1;
        }
    }
    let len = n / blocks;
    for blk in x.chunks_exact_mut(len).take(blocks) {
        if inverse {
            rot_inverse(blk, 1, c, s);
            if stride2 > 0 {
                rot_inverse(blk, stride2, s, c);
            }
        } else {
            if stride2 > 0 {
                rot(blk, stride2, s, c);
            }
            rot(blk, 1, c, s);
        }
    }
}

/// One pass pair of Givens rotations at distance `d` (CELT_SPEC §8.3.4).
fn rot(v: &mut [f32], d: usize, c: f32, s: f32) {
    let l = v.len();
    let step = |v: &mut [f32], k: usize| {
        let (a, e) = (v[k], v[k + d]);
        v[k + d] = c * e + s * a;
        v[k] = c * a - s * e;
    };
    for k in 0..l.saturating_sub(d) {
        step(v, k);
    }
    for k in (0..l.saturating_sub(2 * d)).rev() {
        step(v, k);
    }
}

/// The exact inverse of [`rot`]: the same rotations undone in reverse
/// order.
fn rot_inverse(v: &mut [f32], d: usize, c: f32, s: f32) {
    let l = v.len();
    let undo = |v: &mut [f32], k: usize| {
        let (a, e) = (v[k], v[k + d]);
        v[k] = c * a + s * e;
        v[k + d] = c * e - s * a;
    };
    for k in 0..l.saturating_sub(2 * d) {
        undo(v, k);
    }
    for k in (0..l.saturating_sub(d)).rev() {
        undo(v, k);
    }
}

/// The encoder's pulse search: the vector of `k` pulses whose direction is
/// closest to `x` (largest normalised correlation), by projection and then
/// greedy placement of the remaining pulses.
fn pvq_search(x: &[f32], k: usize, y: &mut [i32]) {
    let n = x.len();
    let ax: Vec<f32> = x.iter().map(|v| v.abs()).collect();
    let sum: f32 = ax.iter().sum();
    y.fill(0);
    if sum.is_nan() || sum <= 1e-15 {
        y[0] = k as i32;
        return;
    }
    let mut placed = 0usize;
    if k > n / 2 {
        let r = (k as f32 - 1.0) / sum;
        for (q, &a) in y.iter_mut().zip(&ax) {
            *q = (a * r).floor() as i32;
            placed += *q as usize;
        }
    }
    let mut xy: f32 = ax.iter().zip(y.iter()).map(|(a, &q)| a * q as f32).sum();
    let mut yy: f32 = y.iter().map(|&q| (q * q) as f32).sum();
    while placed < k {
        let mut best = 0;
        let mut best_num = -1.0f32;
        let mut best_den = 1.0f32;
        for j in 0..n {
            let num = xy + ax[j];
            let num = num * num;
            let den = yy + 2.0 * y[j] as f32 + 1.0;
            if num * best_den > best_num * den {
                best = j;
                best_num = num;
                best_den = den;
            }
        }
        xy += ax[best];
        yy += 2.0 * y[best] as f32 + 1.0;
        y[best] += 1;
        placed += 1;
    }
    for (q, &v) in y.iter_mut().zip(x) {
        if v < 0.0 {
            *q = -*q;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interleave_inverts_deinterleave() {
        for (n0, stride) in [(3, 2), (5, 4), (2, 8), (1, 16)] {
            for hadamard in [false, true] {
                let orig: Vec<f32> = (0..n0 * stride).map(|i| i as f32).collect();
                let mut v = orig.clone();
                deinterleave(&mut v, n0, stride, hadamard);
                interleave(&mut v, n0, stride, hadamard);
                assert_eq!(v, orig);
            }
        }
    }

    #[test]
    fn rotation_inverse_is_exact() {
        for (n, k, blocks, spread) in [(16, 2, 1, 1), (48, 5, 2, 2), (96, 10, 1, 3), (8, 1, 1, 2)] {
            let orig: Vec<f32> = (0..n).map(|i| ((i * 37 % 11) as f32 - 5.0) / 7.0).collect();
            let mut v = orig.clone();
            exp_rotation(&mut v, blocks, k, spread, true);
            exp_rotation(&mut v, blocks, k, spread, false);
            for (a, b) in v.iter().zip(&orig) {
                assert!((a - b).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn isqrt_is_exact() {
        for x in (0..100_000u32).chain([u32::MAX, u32::MAX - 1, 1 << 30]) {
            let r = isqrt(x);
            assert!(
                u64::from(r) * u64::from(r) <= u64::from(x)
                    && u64::from(r + 1) * u64::from(r + 1) > u64::from(x)
            );
        }
    }

    #[test]
    fn cosx_and_log2tan_ranges() {
        for x in 1..16384 {
            let c = cosx(x);
            assert!(
                (200..=32767).contains(&c) || !(64..=16320).contains(&x),
                "cosx({x}) = {c}"
            );
        }
        assert!(log2tan(32767, 200).abs() <= 15059);
    }
}
