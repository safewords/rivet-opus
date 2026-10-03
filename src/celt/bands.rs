//! Band shape coding (RFC 6716 §4.3.4): PVQ with recursive splitting,
//! spreading rotations, time-frequency changes, stereo (mid/side, intensity,
//! dual), folding of lower bands into uncoded ones, and anti-collapse
//! (§4.3.5). Encoder and decoder share this code through [`Coder`]; the
//! encoder additionally searches the codebooks.

use super::cwrs;
use super::mode::{BITRES, Mode, get_pulses, mode};
use super::tables::{
    BIT_DEINTERLEAVE, BIT_INTERLEAVE, EBANDS, NB_EBANDS, ORDERY, SPREAD_FACTOR,
};
use crate::range::Coder;

/// Spread values (§4.3.4.3).
pub const SPREAD_NONE: u32 = 0;
pub const SPREAD_LIGHT: u32 = 1;
pub const SPREAD_NORMAL: u32 = 2;
pub const SPREAD_AGGRESSIVE: u32 = 3;

const EPSILON: f32 = 1e-15;

/// The pseudo-random generator used for folding noise and anti-collapse.
#[inline]
pub fn lcg_rand(seed: u32) -> u32 {
    seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223)
}

/// Scales `x` to the L2 norm `gain`.
pub fn renormalise(x: &mut [f32], gain: f32) {
    let e: f32 = EPSILON + x.iter().map(|v| v * v).sum::<f32>();
    let g = gain / e.sqrt();
    for v in x {
        *v *= g;
    }
}

fn exp_rotation1(x: &mut [f32], len: usize, stride: usize, c: f32, s: f32) {
    let ms = -s;
    for i in 0..len - stride {
        let x1 = x[i];
        let x2 = x[i + stride];
        x[i + stride] = c * x2 + s * x1;
        x[i] = c * x1 + ms * x2;
    }
    if len > 2 * stride {
        for i in (0..=len - 2 * stride - 1).rev() {
            let x1 = x[i];
            let x2 = x[i + stride];
            x[i + stride] = c * x2 + s * x1;
            x[i] = c * x1 + ms * x2;
        }
    }
}

/// §4.3.4.3: the spreading rotation; `dir` +1 before quantization
/// (encoder), -1 after decoding.
pub fn exp_rotation(x: &mut [f32], len: usize, dir: i32, stride: usize, k: usize, spread: u32) {
    if 2 * k >= len || spread == SPREAD_NONE {
        return;
    }
    let factor = SPREAD_FACTOR[spread as usize - 1] as f32;
    let gain = len as f32 / (len as f32 + factor * k as f32);
    let theta = 0.5 * gain * gain;
    let c = (0.5 * std::f32::consts::PI * theta).cos();
    let s = (0.5 * std::f32::consts::PI * (1.0 - theta)).cos();
    let mut stride2 = 0;
    if len >= 8 * stride {
        stride2 = 1;
        while (stride2 * stride2 + stride2) * stride + (stride >> 2) < len {
            stride2 += 1;
        }
    }
    let l = len / stride;
    for i in 0..stride {
        let xs = &mut x[i * l..(i + 1) * l];
        if dir < 0 {
            if stride2 > 0 {
                exp_rotation1(xs, l, stride2, s, c);
            }
            exp_rotation1(xs, l, 1, c, s);
        } else {
            exp_rotation1(xs, l, 1, c, -s);
            if stride2 > 0 {
                exp_rotation1(xs, l, stride2, s, -c);
            }
        }
    }
}

fn extract_collapse_mask(iy: &[i32], n: usize, b: usize) -> u32 {
    if b <= 1 {
        return 1;
    }
    let n0 = n / b;
    let mut mask = 0;
    for i in 0..b {
        if iy[i * n0..(i + 1) * n0].iter().any(|&v| v != 0) {
            mask |= 1 << i;
        }
    }
    mask
}

/// A greedy search for the PVQ codeword of `k` pulses nearest in angle to
/// `x`; returns the codeword.
pub fn pvq_search(x: &[f32], k: usize) -> Vec<i32> {
    let n = x.len();
    let mut iy = vec![0i32; n];
    let ax: Vec<f32> = x.iter().map(|v| v.abs()).collect();
    let mut left = k as i32;
    let mut xy = 0.0f32;
    let mut yy = 0.0f32;
    if k > n / 2 {
        let sum: f32 = ax.iter().sum();
        if sum > 1e-15 {
            let r = (k as f32 - 1.0) / sum;
            for j in 0..n {
                let p = (ax[j] * r).floor() as i32;
                iy[j] = p;
                left -= p;
                xy += ax[j] * p as f32;
                yy += (p * p) as f32;
            }
        } else {
            iy[0] = k as i32;
            left = 0;
        }
    }
    while left > 0 {
        let mut best = 0;
        let mut best_num = -1.0f32;
        let mut best_den = 1.0f32;
        for j in 0..n {
            let num = xy + ax[j];
            let num = num * num;
            let den = yy + 2.0 * iy[j] as f32 + 1.0;
            if num * best_den > best_num * den {
                best_num = num;
                best_den = den;
                best = j;
            }
        }
        xy += ax[best];
        yy += 2.0 * iy[best] as f32 + 1.0;
        iy[best] += 1;
        left -= 1;
    }
    for j in 0..n {
        if x[j] < 0.0 {
            iy[j] = -iy[j];
        }
    }
    iy
}

fn haar1(x: &mut [f32], n0: usize, stride: usize) {
    let n0 = n0 >> 1;
    let s = std::f32::consts::FRAC_1_SQRT_2;
    for i in 0..stride {
        for j in 0..n0 {
            let a = s * x[stride * 2 * j + i];
            let b = s * x[stride * (2 * j + 1) + i];
            x[stride * 2 * j + i] = a + b;
            x[stride * (2 * j + 1) + i] = a - b;
        }
    }
}

fn deinterleave_hadamard(x: &mut [f32], n0: usize, stride: usize, hadamard: bool) {
    let n = n0 * stride;
    let mut tmp = vec![0.0f32; n];
    if hadamard {
        let ord = &ORDERY[stride - 2..];
        for i in 0..stride {
            for j in 0..n0 {
                tmp[ord[i] * n0 + j] = x[j * stride + i];
            }
        }
    } else {
        for i in 0..stride {
            for j in 0..n0 {
                tmp[i * n0 + j] = x[j * stride + i];
            }
        }
    }
    x[..n].copy_from_slice(&tmp);
}

fn interleave_hadamard(x: &mut [f32], n0: usize, stride: usize, hadamard: bool) {
    let n = n0 * stride;
    let mut tmp = vec![0.0f32; n];
    if hadamard {
        let ord = &ORDERY[stride - 2..];
        for i in 0..stride {
            for j in 0..n0 {
                tmp[j * stride + i] = x[ord[i] * n0 + j];
            }
        }
    } else {
        for i in 0..stride {
            for j in 0..n0 {
                tmp[j * stride + i] = x[i * n0 + j];
            }
        }
    }
    x[..n].copy_from_slice(&tmp);
}

/// `cos(pi/2 * x / 16384)` in Q15, bit-exactly.
fn bitexact_cos(x: i32) -> i32 {
    let frac_mul16 = |a: i32, b: i32| (16384 + a * b) >> 15;
    let tmp = (4096 + x * x) >> 13;
    let x2 = tmp;
    let x2 = (32767 - x2) + frac_mul16(x2, -7651 + frac_mul16(x2, 8277 + frac_mul16(-626, x2)));
    1 + x2
}

/// `log2(isin / icos)` in Q11, bit-exactly.
fn bitexact_log2tan(isin: i32, icos: i32) -> i32 {
    let frac_mul16 = |a: i32, b: i32| (16384 + a * b) >> 15;
    let lc = crate::range::ilog(icos as u32);
    let ls = crate::range::ilog(isin as u32);
    let icos = icos << (15 - lc);
    let isin = isin << (15 - ls);
    (ls - lc) * (1 << 11) + frac_mul16(isin, frac_mul16(isin, -2597) + 7932)
        - frac_mul16(icos, frac_mul16(icos, -2597) + 7932)
}

fn compute_qn(n: usize, b: i32, offset: i32, pulse_cap: i32, stereo: bool) -> i32 {
    const EXP2_TABLE8: [i32; 8] = [16384, 17866, 19483, 21247, 23170, 25267, 27554, 30048];
    let mut n2 = 2 * n as i32 - 1;
    if stereo && n == 2 {
        n2 -= 1;
    }
    // Signed division truncating toward zero.
    let qb = ((b + n2 * offset) / n2).min(b - pulse_cap - (4 << BITRES)).min(8 << BITRES);
    if qb < (1 << BITRES >> 1) {
        1
    } else {
        let qn = EXP2_TABLE8[(qb & 7) as usize] >> (14 - (qb >> BITRES));
        ((qn + 1) >> 1) << 1
    }
}

/// Shared state of one frame's band coding.
pub struct BandCtx<'a, E: Coder> {
    pub ec: &'a mut E,
    m: &'static Mode,
    i: usize,
    intensity: usize,
    spread: u32,
    tf_change: i32,
    remaining_bits: i32,
    pub seed: u32,
    resynth: bool,
    disable_inv: bool,
    /// Encoder: linear band amplitudes, for intensity stereo.
    band_e: &'a [f32],
}

struct SplitCtx {
    inv: bool,
    imid: i32,
    iside: i32,
    delta: i32,
    itheta: i32,
    qalloc: i32,
}

impl<E: Coder> BandCtx<'_, E> {
    /// §4.3.4.4 / §4.3.4.5: codes the split angle between two halves (or
    /// between mid and side).
    #[allow(clippy::too_many_arguments)]
    fn compute_theta(
        &mut self,
        x: &mut [f32],
        y: &mut [f32],
        n: usize,
        b: &mut i32,
        bb: usize,
        b0: usize,
        lm: i32,
        stereo: bool,
        fill: &mut u32,
    ) -> SplitCtx {
        let i = self.i;
        let pulse_cap = self.m.log_n[i] + lm * (1 << BITRES);
        let offset = (pulse_cap >> 1) - if stereo && n == 2 { 16 } else { 4 };
        let mut qn = compute_qn(n, *b, offset, pulse_cap, stereo);
        if stereo && i >= self.intensity {
            qn = 1;
        }
        let mut itheta = 0i32;
        if E::ENCODE {
            itheta = stereo_itheta(x, y, stereo, n);
        }
        let tell = self.ec.tell_frac();
        let mut inv = false;
        if qn != 1 {
            if E::ENCODE {
                itheta = (itheta * qn + 8192) >> 14;
            }
            let qn_u = qn as u32;
            if stereo && n > 2 {
                let p0 = 3u32;
                let x0 = qn_u / 2;
                let ft = p0 * (x0 + 1) + x0;
                let xv = if E::ENCODE {
                    itheta as u32
                } else {
                    let fs = self.ec.decode_fs(ft);
                    if fs < (x0 + 1) * p0 { fs / p0 } else { x0 + 1 + (fs - (x0 + 1) * p0) }
                };
                let (fl, fh) = if xv <= x0 {
                    (p0 * xv, p0 * (xv + 1))
                } else {
                    ((xv - 1 - x0) + (x0 + 1) * p0, (xv - x0) + (x0 + 1) * p0)
                };
                self.ec.code(fl, fh, ft);
                itheta = xv as i32;
            } else if b0 > 1 || stereo {
                itheta = self.ec.uint(itheta as u32, qn_u + 1) as i32;
            } else {
                let half = qn_u >> 1;
                let ft = (half + 1) * (half + 1);
                let (fl, fs, it) = if E::ENCODE {
                    let it = itheta as u32;
                    let fs = if it <= half { it + 1 } else { qn_u + 1 - it };
                    let fl = if it <= half { it * (it + 1) >> 1 } else { ft - ((qn_u + 1 - it) * (qn_u + 2 - it) >> 1) };
                    (fl, fs, it)
                } else {
                    let fm = self.ec.decode_fs(ft);
                    if fm < (half * (half + 1) >> 1) {
                        let it = (isqrt32(8 * fm + 1) - 1) >> 1;
                        (it * (it + 1) >> 1, it + 1, it)
                    } else {
                        let it = (2 * (qn_u + 1) - isqrt32(8 * (ft - fm - 1) + 1)) >> 1;
                        (ft - ((qn_u + 1 - it) * (qn_u + 2 - it) >> 1), qn_u + 1 - it, it)
                    }
                };
                self.ec.code(fl, fl + fs, ft);
                itheta = it as i32;
            }
            itheta = itheta * 16384 / qn;
            if E::ENCODE && stereo {
                if itheta == 0 {
                    intensity_stereo(self.band_e, x, y, i, n);
                } else if itheta != 16384 {
                    stereo_split(x, y, n);
                }
            }
        } else if stereo {
            if E::ENCODE {
                inv = itheta > 8192 && !self.disable_inv;
                if inv {
                    for v in y.iter_mut().take(n) {
                        *v = -*v;
                    }
                }
                intensity_stereo(self.band_e, x, y, i, n);
            }
            if *b > 2 << BITRES && self.remaining_bits > 2 << BITRES {
                inv = self.ec.bit_logp(inv, 2);
            } else {
                inv = false;
            }
            if self.disable_inv {
                inv = false;
            }
            itheta = 0;
        }
        let qalloc = self.ec.tell_frac() - tell;
        *b -= qalloc;
        let (imid, iside, delta);
        if itheta == 0 {
            imid = 32767;
            iside = 0;
            *fill &= (1 << bb) - 1;
            delta = -16384;
        } else if itheta == 16384 {
            imid = 0;
            iside = 32767;
            *fill &= ((1 << bb) - 1) << bb;
            delta = 16384;
        } else {
            imid = bitexact_cos(itheta);
            iside = bitexact_cos(16384 - itheta);
            delta = (16384 + (((n as i32 - 1) << 7) * bitexact_log2tan(iside, imid))) >> 15;
        }
        SplitCtx { inv, imid, iside, delta, itheta, qalloc }
    }

    fn quant_band_n1(&mut self, x: &mut [f32], y: Option<&mut [f32]>, lowband_out: Option<&mut [f32]>) -> u32 {
        let code = |v: &mut f32, ctx: &mut Self| {
            let mut sign = false;
            if ctx.remaining_bits >= 1 << BITRES {
                sign = ctx.ec.bits(u32::from(*v < 0.0), 1) != 0;
                ctx.remaining_bits -= 1 << BITRES;
            }
            if ctx.resynth {
                *v = if sign { -1.0 } else { 1.0 };
            }
        };
        code(&mut x[0], self);
        if let Some(y) = y {
            code(&mut y[0], self);
        }
        if let Some(lo) = lowband_out {
            lo[0] = x[0];
        }
        1
    }

    /// Codes a PVQ vector or splits it in two (§4.3.4.4).
    #[allow(clippy::too_many_arguments)]
    fn quant_partition(
        &mut self,
        x: &mut [f32],
        n: usize,
        mut b: i32,
        mut bb: usize,
        lowband: Option<&[f32]>,
        mut lm: i32,
        gain: f32,
        mut fill: u32,
    ) -> u32 {
        let i = self.i;
        let b0 = bb;
        let cache = self.m.cache(i, lm);
        let max_bits = i32::from(cache[usize::from(cache[0])]);
        if lm != -1 && b > max_bits + 12 && n > 2 {
            let n = n >> 1;
            let (xa, ya) = x.split_at_mut(n);
            lm -= 1;
            if bb == 1 {
                fill = (fill & 1) | (fill << 1);
            }
            bb = (bb + 1) >> 1;
            let sctx = self.compute_theta(xa, ya, n, &mut b, bb, b0, lm, false, &mut fill);
            let mid = sctx.imid as f32 / 32768.0;
            let side = sctx.iside as f32 / 32768.0;
            let mut delta = sctx.delta;
            let itheta = sctx.itheta;
            if b0 > 1 && (itheta & 0x3fff) != 0 {
                if itheta > 8192 {
                    delta -= delta >> (4 - lm);
                } else {
                    delta = (delta + ((n as i32) << BITRES >> (5 - lm))).min(0);
                }
            }
            let mut mbits = 0.max(b.min((b - delta) / 2));
            let mut sbits = b - mbits;
            self.remaining_bits -= sctx.qalloc;
            let (low1, low2) = match lowband {
                Some(l) => (Some(&l[..n]), Some(&l[n..2 * n])),
                None => (None, None),
            };
            let rebalance = self.remaining_bits;
            let mut cm;
            if mbits >= sbits {
                cm = self.quant_partition(xa, n, mbits, bb, low1, lm, gain * mid, fill);
                let rebalance = mbits - (rebalance - self.remaining_bits);
                if rebalance > 3 << BITRES && itheta != 0 {
                    sbits += rebalance - (3 << BITRES);
                }
                cm |= self.quant_partition(ya, n, sbits, bb, low2, lm, gain * side, fill >> bb) << (b0 >> 1);
            } else {
                cm = self.quant_partition(ya, n, sbits, bb, low2, lm, gain * side, fill >> bb) << (b0 >> 1);
                let rebalance = sbits - (rebalance - self.remaining_bits);
                if rebalance > 3 << BITRES && itheta != 16384 {
                    mbits += rebalance - (3 << BITRES);
                }
                cm |= self.quant_partition(xa, n, mbits, bb, low1, lm, gain * mid, fill);
            }
            return cm;
        }
        let mut q = self.m.bits2pulses(i, lm, b);
        let mut curr_bits = self.m.pulses2bits(i, lm, q);
        self.remaining_bits -= curr_bits;
        while self.remaining_bits < 0 && q > 0 {
            self.remaining_bits += curr_bits;
            q -= 1;
            curr_bits = self.m.pulses2bits(i, lm, q);
            self.remaining_bits -= curr_bits;
        }
        if q != 0 {
            let k = get_pulses(q);
            let x = &mut x[..n];
            if E::ENCODE {
                exp_rotation(x, n, 1, bb, k, self.spread);
                let iy = pvq_search(x, k);
                let idx = cwrs::encode(&iy, k);
                self.ec.uint(idx, cwrs::v(n, k) as u32);
                if self.resynth {
                    normalise_residual(&iy, x, gain);
                    exp_rotation(x, n, -1, bb, k, self.spread);
                }
                extract_collapse_mask(&iy, n, bb)
            } else {
                let idx = self.ec.uint(0, cwrs::v(n, k) as u32);
                let mut iy = vec![0i32; n];
                cwrs::decode(idx, n, k, &mut iy);
                normalise_residual(&iy, x, gain);
                exp_rotation(x, n, -1, bb, k, self.spread);
                extract_collapse_mask(&iy, n, bb)
            }
        } else {
            let mut cm = 0;
            if self.resynth {
                let cm_mask = (1u32 << bb) - 1;
                fill &= cm_mask;
                let x = &mut x[..n];
                if fill == 0 {
                    x.fill(0.0);
                } else {
                    match lowband {
                        None => {
                            for v in x.iter_mut() {
                                self.seed = lcg_rand(self.seed);
                                *v = ((self.seed as i32) >> 20) as f32;
                            }
                            cm = cm_mask;
                        }
                        Some(lb) => {
                            for (v, &l) in x.iter_mut().zip(lb) {
                                self.seed = lcg_rand(self.seed);
                                let tmp = if self.seed & 0x8000 != 0 { 1.0 / 256.0 } else { -1.0 / 256.0 };
                                *v = l + tmp;
                            }
                            cm = fill;
                        }
                    }
                    renormalise(x, gain);
                }
            }
            cm
        }
    }

    /// Codes one band of one channel (§4.3.4), with its TF change.
    #[allow(clippy::too_many_arguments)]
    fn quant_band(
        &mut self,
        x: &mut [f32],
        n: usize,
        b: i32,
        mut bb: usize,
        lowband: Option<&[f32]>,
        lm: i32,
        lowband_out: Option<&mut [f32]>,
        gain: f32,
        mut fill: u32,
    ) -> u32 {
        let n0 = n;
        let mut n_b = n / bb;
        let long_blocks = bb == 1;
        let mut tf_change = self.tf_change;
        if n == 1 {
            return self.quant_band_n1(x, None, lowband_out);
        }
        let mut recombine = 0;
        if tf_change > 0 {
            recombine = tf_change;
        }
        let mut lowband_buf: Option<Vec<f32>> = lowband.map(|l| l[..n].to_vec());
        for k in 0..recombine {
            if E::ENCODE {
                haar1(x, n >> k, 1 << k);
            }
            if let Some(l) = lowband_buf.as_mut() {
                haar1(l, n >> k, 1 << k);
            }
            fill = BIT_INTERLEAVE[(fill & 0xF) as usize] | (BIT_INTERLEAVE[(fill >> 4) as usize] << 2);
        }
        bb >>= recombine;
        n_b <<= recombine;
        let mut time_divide = 0;
        while n_b & 1 == 0 && tf_change < 0 {
            if E::ENCODE {
                haar1(x, n_b, bb);
            }
            if let Some(l) = lowband_buf.as_mut() {
                haar1(l, n_b, bb);
            }
            fill |= fill << bb;
            bb <<= 1;
            n_b >>= 1;
            time_divide += 1;
            tf_change += 1;
        }
        let b0 = bb;
        let n_b0 = n_b;
        if b0 > 1 {
            if E::ENCODE {
                deinterleave_hadamard(x, n_b >> recombine, b0 << recombine, long_blocks);
            }
            if let Some(l) = lowband_buf.as_mut() {
                deinterleave_hadamard(l, n_b >> recombine, b0 << recombine, long_blocks);
            }
        }
        let mut cm = self.quant_partition(x, n, b, bb, lowband_buf.as_deref(), lm, gain, fill);
        if self.resynth {
            if b0 > 1 {
                interleave_hadamard(x, n_b >> recombine, b0 << recombine, long_blocks);
            }
            let mut n_b = n_b0;
            bb = b0;
            for _ in 0..time_divide {
                bb >>= 1;
                n_b <<= 1;
                cm |= cm >> bb;
                haar1(x, n_b, bb);
            }
            for k in 0..recombine {
                cm = BIT_DEINTERLEAVE[cm as usize];
                haar1(x, n0 >> k, 1 << k);
            }
            bb <<= recombine;
            if let Some(out) = lowband_out {
                let s = (n0 as f32).sqrt();
                for j in 0..n0 {
                    out[j] = s * x[j];
                }
            }
            cm &= (1 << bb) - 1;
        }
        cm
    }

    /// Codes one band of a stereo pair jointly (§4.3.4.4).
    #[allow(clippy::too_many_arguments)]
    fn quant_band_stereo(
        &mut self,
        x: &mut [f32],
        y: &mut [f32],
        n: usize,
        mut b: i32,
        bb: usize,
        lowband: Option<&[f32]>,
        lm: i32,
        lowband_out: Option<&mut [f32]>,
        mut fill: u32,
    ) -> u32 {
        let orig_fill = fill;
        if n == 1 {
            return self.quant_band_n1(x, Some(y), lowband_out);
        }
        let sctx = self.compute_theta(x, y, n, &mut b, bb, bb, lm, true, &mut fill);
        let mid = sctx.imid as f32 / 32768.0;
        let side = sctx.iside as f32 / 32768.0;
        let itheta = sctx.itheta;
        let mut cm;
        if n == 2 {
            let mut mbits = b;
            let mut sbits = 0;
            if itheta != 0 && itheta != 16384 {
                sbits = 1 << BITRES;
            }
            mbits -= sbits;
            let c = itheta > 8192;
            self.remaining_bits -= sctx.qalloc + sbits;
            let (x2, y2): (&mut [f32], &mut [f32]) = if c { (&mut *y, &mut *x) } else { (&mut *x, &mut *y) };
            let mut sign = 0;
            if sbits != 0 {
                let want = u32::from(x2[0] * y2[1] - x2[1] * y2[0] < 0.0);
                sign = self.ec.bits(want, 1);
            }
            let sign = 1.0 - 2.0 * sign as f32;
            cm = self.quant_band(x2, n, mbits, bb, lowband, lm, lowband_out, 1.0, orig_fill);
            y2[0] = -sign * x2[1];
            y2[1] = sign * x2[0];
            if self.resynth {
                x[0] *= mid;
                x[1] *= mid;
                y[0] *= side;
                y[1] *= side;
                let t = x[0];
                x[0] = t - y[0];
                y[0] += t;
                let t = x[1];
                x[1] = t - y[1];
                y[1] += t;
            }
        } else {
            let delta = sctx.delta;
            let mut mbits = 0.max(b.min((b - delta) / 2));
            let mut sbits = b - mbits;
            self.remaining_bits -= sctx.qalloc;
            let rebalance = self.remaining_bits;
            if mbits >= sbits {
                cm = self.quant_band(x, n, mbits, bb, lowband, lm, lowband_out, 1.0, fill);
                let rebalance = mbits - (rebalance - self.remaining_bits);
                if rebalance > 3 << BITRES && itheta != 0 {
                    sbits += rebalance - (3 << BITRES);
                }
                cm |= self.quant_band(y, n, sbits, bb, None, lm, None, side, fill >> bb);
            } else {
                cm = self.quant_band(y, n, sbits, bb, None, lm, None, side, fill >> bb);
                let rebalance = sbits - (rebalance - self.remaining_bits);
                if rebalance > 3 << BITRES && itheta != 16384 {
                    mbits += rebalance - (3 << BITRES);
                }
                cm |= self.quant_band(x, n, mbits, bb, lowband, lm, lowband_out, 1.0, fill);
            }
        }
        if self.resynth {
            if n != 2 {
                stereo_merge(x, y, mid, n);
            }
            if sctx.inv {
                for v in y.iter_mut().take(n) {
                    *v = -*v;
                }
            }
        }
        cm
    }
}

fn isqrt32(v: u32) -> u32 {
    let mut r = (v as f64).sqrt() as u32;
    while r * r > v {
        r -= 1;
    }
    while (r + 1) * (r + 1) <= v {
        r += 1;
    }
    r
}

fn normalise_residual(iy: &[i32], x: &mut [f32], gain: f32) {
    let ryy: f32 = iy.iter().map(|&v| (v * v) as f32).sum();
    let g = gain / ryy.sqrt();
    for (xv, &yv) in x.iter_mut().zip(iy) {
        *xv = g * yv as f32;
    }
}

fn stereo_merge(x: &mut [f32], y: &mut [f32], mid: f32, n: usize) {
    let mut xp = 0.0f32;
    let mut side = 0.0f32;
    for j in 0..n {
        xp += y[j] * x[j];
        side += y[j] * y[j];
    }
    xp *= mid;
    let mid2 = mid * mid;
    let el = mid2 + side - 2.0 * xp;
    let er = mid2 + side + 2.0 * xp;
    if er < 6e-4 || el < 6e-4 {
        y[..n].copy_from_slice(&x[..n]);
        return;
    }
    let lgain = 1.0 / el.sqrt();
    let rgain = 1.0 / er.sqrt();
    for j in 0..n {
        let l = mid * x[j];
        let r = y[j];
        x[j] = lgain * (l - r);
        y[j] = rgain * (l + r);
    }
}

fn stereo_itheta(x: &[f32], y: &[f32], stereo: bool, n: usize) -> i32 {
    let (mut emid, mut eside) = (EPSILON, EPSILON);
    if stereo {
        for i in 0..n {
            let m = x[i] + y[i];
            let s = x[i] - y[i];
            emid += m * m;
            eside += s * s;
        }
    } else {
        for i in 0..n {
            emid += x[i] * x[i];
            eside += y[i] * y[i];
        }
    }
    let a = eside.sqrt().atan2(emid.sqrt());
    (0.5 + 16384.0 * std::f32::consts::FRAC_2_PI * a).floor() as i32
}

fn intensity_stereo(band_e: &[f32], x: &mut [f32], y: &[f32], i: usize, n: usize) {
    let left = band_e[i];
    let right = band_e[i + NB_EBANDS];
    let norm = EPSILON + (EPSILON + left * left + right * right).sqrt();
    let a1 = left / norm;
    let a2 = right / norm;
    for j in 0..n {
        x[j] = a1 * x[j] + a2 * y[j];
    }
}

fn stereo_split(x: &mut [f32], y: &mut [f32], n: usize) {
    let s = std::f32::consts::FRAC_1_SQRT_2;
    for j in 0..n {
        let l = s * x[j];
        let r = s * y[j];
        x[j] = l + r;
        y[j] = r - l;
    }
}

/// Everything [`quant_all_bands`] needs about the frame.
pub struct FrameParams<'a> {
    pub start: usize,
    pub end: usize,
    pub lm: usize,
    /// Number of short MDCTs (1 for a long block).
    pub short_blocks: usize,
    pub spread: u32,
    pub dual_stereo: bool,
    pub intensity: usize,
    pub tf_res: &'a [i32; NB_EBANDS],
    /// Total bits for the frame in 1/8 bits (less the anti-collapse bit).
    pub total_bits: i32,
    pub balance: i32,
    pub pulses: &'a [i32; NB_EBANDS],
    pub coded_bands: usize,
    pub disable_inv: bool,
    /// Whether to reconstruct the normalized spectrum (always on decode).
    pub resynth: bool,
}

/// Codes the normalized shapes of all bands (§4.3.4). `x` and `y` hold the
/// channels' normalized spectra (`y` empty for mono); on return they hold
/// the decoded shapes. Returns the collapse masks `[band * C + c]` and the
/// updated seed.
pub fn quant_all_bands<E: Coder>(
    ec: &mut E,
    p: &FrameParams,
    x: &mut [f32],
    y: &mut [f32],
    band_e: &[f32],
    seed: u32,
) -> (Vec<u32>, u32) {
    let m = mode();
    let lm = p.lm;
    let mm = 1usize << lm;
    let bb = if p.short_blocks > 1 { mm } else { 1 };
    let c = if y.is_empty() { 1 } else { 2 };
    let norm_offset = mm * EBANDS[p.start];
    let norm_len = mm * EBANDS[NB_EBANDS - 1] - norm_offset;
    let mut norm = vec![0.0f32; norm_len.max(1)];
    let mut norm2 = vec![0.0f32; if c == 2 { norm_len.max(1) } else { 0 }];
    let mut collapse = vec![0u32; NB_EBANDS * c];
    let mut lowband_offset = 0usize;
    let mut update_lowband = true;
    let mut balance = p.balance;
    let mut dual_stereo = p.dual_stereo;
    let mut ctx = BandCtx {
        ec,
        m,
        i: 0,
        intensity: p.intensity,
        spread: p.spread,
        tf_change: 0,
        remaining_bits: 0,
        seed,
        resynth: p.resynth,
        disable_inv: p.disable_inv,
        band_e,
    };
    for i in p.start..p.end {
        ctx.i = i;
        let last = i == p.end - 1;
        let lo = mm * EBANDS[i];
        let n = mm * EBANDS[i + 1] - lo;
        let tell = ctx.ec.tell_frac();
        if i != p.start {
            balance -= tell;
        }
        let remaining_bits = p.total_bits - tell - 1;
        ctx.remaining_bits = remaining_bits;
        let b = if i < p.coded_bands {
            let curr_balance = balance / (p.coded_bands - i).min(3) as i32;
            0.max(16383.min((remaining_bits + 1).min(p.pulses[i] + curr_balance)))
        } else {
            0
        };
        if p.resynth
            && (lo as isize - n as isize >= (mm * EBANDS[p.start]) as isize || i == p.start + 1)
            && (update_lowband || lowband_offset == 0)
        {
            lowband_offset = i;
        }
        if i == p.start + 1 {
            // RFC 8251 §9: repeat part of the first band so the second can
            // always fold.
            let n1 = mm * (EBANDS[p.start + 1] - EBANDS[p.start]);
            let n2 = mm * (EBANDS[p.start + 2] - EBANDS[p.start + 1]);
            if n2 > n1 {
                let src = 2 * n1 - n2;
                norm.copy_within(src..src + n2 - n1, n1);
                if c == 2 {
                    norm2.copy_within(src..src + n2 - n1, n1);
                }
            }
        }
        ctx.tf_change = p.tf_res[i];
        let mut effective_lowband: Option<usize> = None;
        let (mut x_cm, mut y_cm);
        if lowband_offset != 0 && (p.spread != SPREAD_AGGRESSIVE || bb > 1 || ctx.tf_change < 0) {
            let el = (mm * EBANDS[lowband_offset]).saturating_sub(norm_offset + n);
            effective_lowband = Some(el);
            let mut fold_start = lowband_offset;
            loop {
                fold_start -= 1;
                if mm * EBANDS[fold_start] <= el + norm_offset {
                    break;
                }
            }
            let mut fold_end = lowband_offset - 1;
            loop {
                fold_end += 1;
                if !(fold_end < i && mm * EBANDS[fold_end] < el + norm_offset + n) {
                    break;
                }
            }
            x_cm = 0;
            y_cm = 0;
            let mut fold_i = fold_start;
            loop {
                x_cm |= collapse[fold_i * c];
                y_cm |= collapse[fold_i * c + c - 1];
                fold_i += 1;
                if fold_i >= fold_end {
                    break;
                }
            }
        } else {
            x_cm = (1 << bb) - 1;
            y_cm = x_cm;
        }
        if dual_stereo && i == p.intensity {
            dual_stereo = false;
            if p.resynth {
                for j in 0..lo - norm_offset {
                    norm[j] = 0.5 * (norm[j] + norm2[j]);
                }
            }
        }
        let out_start = lo - norm_offset;
        let xs = &mut x[lo..lo + n];
        if dual_stereo {
            let ys = &mut y[lo..lo + n];
            let lowband: Option<Vec<f32>> = effective_lowband.map(|e| norm[e..e + n].to_vec());
            let mut out = vec![0.0f32; n];
            x_cm = ctx.quant_band(xs, n, b / 2, bb, lowband.as_deref(), lm as i32, if last { None } else { Some(&mut out) }, 1.0, x_cm);
            if !last && out_start + n <= norm.len() {
                norm[out_start..out_start + n].copy_from_slice(&out);
            }
            let lowband2: Option<Vec<f32>> = effective_lowband.map(|e| norm2[e..e + n].to_vec());
            y_cm = ctx.quant_band(ys, n, b / 2, bb, lowband2.as_deref(), lm as i32, if last { None } else { Some(&mut out) }, 1.0, y_cm);
            if !last && out_start + n <= norm2.len() {
                norm2[out_start..out_start + n].copy_from_slice(&out);
            }
        } else {
            let lowband: Option<Vec<f32>> = effective_lowband.map(|e| norm[e..e + n].to_vec());
            let mut out = vec![0.0f32; n];
            if c == 2 {
                let ys = &mut y[lo..lo + n];
                x_cm = ctx.quant_band_stereo(xs, ys, n, b, bb, lowband.as_deref(), lm as i32, if last { None } else { Some(&mut out) }, x_cm | y_cm);
            } else {
                x_cm = ctx.quant_band(xs, n, b, bb, lowband.as_deref(), lm as i32, if last { None } else { Some(&mut out) }, 1.0, x_cm | y_cm);
            }
            if !last && out_start + n <= norm.len() {
                norm[out_start..out_start + n].copy_from_slice(&out);
            }
            y_cm = x_cm;
        }
        collapse[i * c] = x_cm;
        collapse[i * c + c - 1] = y_cm;
        balance += p.pulses[i] + tell;
        update_lowband = b > (n as i32) << BITRES;
    }
    (collapse, ctx.seed)
}

/// §4.3.5: fills collapsed short blocks with noise at an energy derived
/// from the two previous frames.
#[allow(clippy::too_many_arguments)]
pub fn anti_collapse(
    x: &mut [f32],
    size: usize,
    collapse: &[u32],
    lm: usize,
    c: usize,
    start: usize,
    end: usize,
    log_e: &[f32],
    prev1: &[f32],
    prev2: &[f32],
    pulses: &[i32; NB_EBANDS],
    mut seed: u32,
) {
    for i in start..end {
        let n0 = EBANDS[i + 1] - EBANDS[i];
        let depth = ((1 + pulses[i]) as usize / n0) >> lm;
        let thresh = 0.5 * (-0.125 * depth as f32).exp2();
        let sqrt_1 = 1.0 / ((n0 << lm) as f32).sqrt();
        for ch in 0..c {
            let mut p1 = prev1[ch * NB_EBANDS + i];
            let mut p2 = prev2[ch * NB_EBANDS + i];
            if c == 1 {
                p1 = p1.max(prev1[NB_EBANDS + i]);
                p2 = p2.max(prev2[NB_EBANDS + i]);
            }
            let ediff = (log_e[ch * NB_EBANDS + i] - p1.min(p2)).max(0.0);
            let mut r = 2.0 * (-ediff).exp2();
            if lm == 3 {
                r *= std::f32::consts::SQRT_2;
            }
            r = r.min(thresh) * sqrt_1;
            let band = &mut x[ch * size + (EBANDS[i] << lm)..ch * size + (EBANDS[i + 1] << lm)];
            let mut renorm = false;
            for k in 0..1usize << lm {
                if collapse[i * c + ch] & (1 << k) == 0 {
                    for j in 0..n0 {
                        seed = lcg_rand(seed);
                        band[(j << lm) + k] = if seed & 0x8000 != 0 { r } else { -r };
                    }
                    renorm = true;
                }
            }
            if renorm {
                renormalise(band, 1.0);
            }
        }
    }
}
