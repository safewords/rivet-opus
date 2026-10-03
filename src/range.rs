//! The range coder of RFC 6716 Â§4.1 (decoding) and Â§5.1 (encoding).
//!
//! Both sides keep the RFC's state exactly â€” `rng`, `val`, the whole-bit
//! counter `nbits_total` â€” because the CELT bit allocation is driven by
//! [`RangeDecoder::tell_frac`] and must agree bit for bit between encoder and
//! decoder. Raw bits (Â§4.1.4) are packed backwards from the end of the frame.

/// `ilog(n)` of RFC 6716 Â§1.1.10: the number of bits needed to write `n`
/// (0 for 0).
#[inline]
pub(crate) fn ilog(n: u32) -> i32 {
    32 - n.leading_zeros() as i32
}

/// The range decoder of RFC 6716 Â§4.1.
#[derive(Clone)]
pub struct RangeDecoder<'a> {
    buf: &'a [u8],
    /// Bytes of `buf` that belong to this frame (Â§4.5.1.3 may shrink it).
    storage: usize,
    /// Bytes read from the front.
    offs: usize,
    /// Bytes read from the back by the raw-bit reader.
    end_offs: usize,
    end_window: u32,
    nend_bits: i32,
    nbits_total: i32,
    rng: u32,
    val: u32,
    /// The last byte read; its low bit is still to be used.
    rem: u32,
    /// The divisor saved by [`Self::decode`] for [`Self::update`].
    ext: u32,
    error: bool,
}

impl<'a> RangeDecoder<'a> {
    /// Starts decoding `buf` (Â§4.1.1).
    pub fn new(buf: &'a [u8]) -> Self {
        let mut d = Self {
            buf,
            storage: buf.len(),
            offs: 0,
            end_offs: 0,
            end_window: 0,
            nend_bits: 0,
            nbits_total: 9,
            rng: 128,
            val: 0,
            rem: 0,
            ext: 0,
            error: false,
        };
        d.rem = d.read_byte();
        d.val = 127 - (d.rem >> 1);
        d.normalize();
        d
    }

    #[inline]
    fn read_byte(&mut self) -> u32 {
        if self.offs < self.storage {
            let b = self.buf[self.offs];
            self.offs += 1;
            u32::from(b)
        } else {
            0
        }
    }

    #[inline]
    fn read_byte_from_end(&mut self) -> u32 {
        if self.end_offs < self.storage {
            self.end_offs += 1;
            u32::from(self.buf[self.storage - self.end_offs])
        } else {
            0
        }
    }

    /// Â§4.1.2.1.
    #[inline]
    fn normalize(&mut self) {
        while self.rng <= 1 << 23 {
            self.nbits_total += 8;
            self.rng <<= 8;
            let prev = self.rem;
            self.rem = self.read_byte();
            let sym = ((prev << 8) | self.rem) >> 1;
            self.val = ((self.val << 8).wrapping_add(255 & !sym)) & 0x7FFF_FFFF;
        }
    }

    /// The first step of Â§4.1.2: the value `fs` in `[0, ft)` that locates the
    /// next symbol.
    #[inline]
    pub fn decode(&mut self, ft: u32) -> u32 {
        self.ext = self.rng / ft;
        let s = self.val / self.ext;
        ft - (s + 1).min(ft)
    }

    /// [`Self::decode`] with `ft = 1 << bits` (Â§4.1.3.1).
    #[inline]
    pub fn decode_bin(&mut self, bits: u32) -> u32 {
        self.ext = self.rng >> bits;
        let s = self.val / self.ext;
        (1u32 << bits) - (s + 1).min(1u32 << bits)
    }

    /// The second step of Â§4.1.2, with the symbol's `(fl, fh, ft)`.
    #[inline]
    pub fn update(&mut self, fl: u32, fh: u32, ft: u32) {
        let s = self.ext * (ft - fh);
        self.val -= s;
        self.rng = if fl > 0 { self.ext * (fh - fl) } else { self.rng - s };
        self.normalize();
    }

    /// One binary symbol whose "1" has probability `2^-logp` (Â§4.1.3.2).
    #[inline]
    pub fn bit_logp(&mut self, logp: u32) -> bool {
        let s = self.rng >> logp;
        let one = self.val < s;
        if one {
            self.rng = s;
        } else {
            self.val -= s;
            self.rng -= s;
        }
        self.normalize();
        one
    }

    /// A symbol from an inverse-CDF table with `ft = 1 << ftb` (Â§4.1.3.3).
    #[inline]
    pub fn icdf(&mut self, icdf: &[u8], ftb: u32) -> usize {
        let r = self.rng >> ftb;
        let mut s = self.rng;
        let mut t;
        let mut k = 0usize;
        loop {
            t = s;
            s = r * u32::from(icdf[k]);
            if self.val >= s {
                break;
            }
            k += 1;
        }
        self.val -= s;
        self.rng = t - s;
        self.normalize();
        k
    }

    /// A uniformly distributed integer in `[0, ft)` (Â§4.1.5). An index past
    /// the end marks the frame corrupt and saturates to `ft - 1`.
    pub fn uint(&mut self, ft: u32) -> u32 {
        debug_assert!(ft > 1);
        let ftm1 = ft - 1;
        let ftb = ilog(ftm1);
        if ftb > 8 {
            let sh = (ftb - 8) as u32;
            let ft1 = (ftm1 >> sh) + 1;
            let s = self.decode(ft1);
            self.update(s, s + 1, ft1);
            let t = (s << sh) | self.bits(sh);
            if t <= ftm1 {
                t
            } else {
                self.error = true;
                ftm1
            }
        } else {
            let s = self.decode(ft);
            self.update(s, s + 1, ft);
            s
        }
    }

    /// `n` raw bits from the end of the frame (Â§4.1.4), `n <= 25`.
    pub fn bits(&mut self, n: u32) -> u32 {
        if n == 0 {
            return 0;
        }
        let mut window = self.end_window;
        let mut available = self.nend_bits;
        if available < n as i32 {
            loop {
                window |= self.read_byte_from_end() << available;
                available += 8;
                if available > 24 {
                    break;
                }
            }
        }
        let ret = window & ((1u32 << n) - 1);
        window >>= n;
        available -= n as i32;
        self.end_window = window;
        self.nend_bits = available;
        self.nbits_total += n as i32;
        ret
    }

    /// Bits used so far, rounded up (Â§4.1.6.1).
    #[inline]
    pub fn tell(&self) -> i32 {
        self.nbits_total - ilog(self.rng)
    }

    /// Bits used so far in 1/8 bits (Â§4.1.6.2).
    pub fn tell_frac(&self) -> i32 {
        tell_frac(self.nbits_total, self.rng)
    }

    /// The current range: after the last symbol of a frame, the "final range"
    /// a conforming decoder must reproduce (RFC 6716 Â§6).
    pub fn range(&self) -> u32 {
        self.rng
    }

    /// The frame's size in bytes (as shrunk by [`Self::shrink`]).
    pub fn storage(&self) -> usize {
        self.storage
    }

    /// Drops `n` bytes from the end of the frame (Â§4.5.1.3: the redundant
    /// CELT frame is cut off before the CELT layer reads its raw bits).
    pub fn shrink(&mut self, n: usize) {
        self.storage = self.storage.saturating_sub(n);
    }

    /// Whether an out-of-range value was decoded.
    pub fn error(&self) -> bool {
        self.error
    }

    /// Counts `n` more bits as used (a CELT silence frame consumes the rest
    /// of the frame, §4.3).
    pub fn add_bits(&mut self, n: i32) {
        self.nbits_total += n.max(0);
    }
}

fn tell_frac(nbits_total: i32, rng: u32) -> i32 {
    let nbits = nbits_total << 3;
    let mut l = ilog(rng);
    let mut r = rng >> (l - 16);
    for _ in 0..3 {
        r = (r * r) >> 15;
        let b = (r >> 16) as i32;
        l = 2 * l + b;
        r >>= b;
    }
    nbits - l
}

/// The range encoder of RFC 6716 Â§5.1, writing into a buffer of fixed size.
#[derive(Clone)]
pub struct RangeEncoder {
    buf: Vec<u8>,
    offs: usize,
    end_offs: usize,
    end_window: u32,
    nend_bits: i32,
    nbits_total: i32,
    rng: u32,
    val: u32,
    rem: i32,
    ext: u32,
    error: bool,
}

impl RangeEncoder {
    /// An encoder for a frame of exactly `size` bytes.
    pub fn new(size: usize) -> Self {
        Self {
            buf: vec![0; size],
            offs: 0,
            end_offs: 0,
            end_window: 0,
            nend_bits: 0,
            nbits_total: 33,
            rng: 1 << 31,
            val: 0,
            rem: -1,
            ext: 0,
            error: false,
        }
    }

    fn write_byte(&mut self, b: u32) {
        if self.offs + self.end_offs >= self.buf.len() {
            self.error = true;
            return;
        }
        self.buf[self.offs] = b as u8;
        self.offs += 1;
    }

    fn write_byte_at_end(&mut self, b: u32) {
        if self.offs + self.end_offs >= self.buf.len() {
            self.error = true;
            return;
        }
        self.end_offs += 1;
        let n = self.buf.len();
        self.buf[n - self.end_offs] = b as u8;
    }

    /// Â§5.1.1.2.
    fn carry_out(&mut self, c: u32) {
        if c == 255 {
            self.ext += 1;
            return;
        }
        let carry = c >> 8;
        if self.rem >= 0 {
            self.write_byte(self.rem as u32 + carry);
        }
        if self.ext > 0 {
            let sym = (255 + carry) & 255;
            for _ in 0..self.ext {
                self.write_byte(sym);
            }
            self.ext = 0;
        }
        self.rem = (c & 255) as i32;
    }

    /// Â§5.1.1.1.
    #[inline]
    fn normalize(&mut self) {
        while self.rng <= 1 << 23 {
            self.carry_out(self.val >> 23);
            self.val = (self.val << 8) & 0x7FFF_FFFF;
            self.rng <<= 8;
            self.nbits_total += 8;
        }
    }

    /// Encodes the symbol `(fl, fh, ft)` (Â§5.1.1).
    #[inline]
    pub fn encode(&mut self, fl: u32, fh: u32, ft: u32) {
        let r = self.rng / ft;
        if fl > 0 {
            self.val += self.rng - r * (ft - fl);
            self.rng = r * (fh - fl);
        } else {
            self.rng -= r * (ft - fh);
        }
        self.normalize();
    }

    /// [`Self::encode`] with `ft = 1 << bits`.
    #[inline]
    pub fn encode_bin(&mut self, fl: u32, fh: u32, bits: u32) {
        let r = self.rng >> bits;
        if fl > 0 {
            self.val += self.rng - r * ((1 << bits) - fl);
            self.rng = r * (fh - fl);
        } else {
            self.rng -= r * ((1 << bits) - fh);
        }
        self.normalize();
    }

    /// One binary symbol whose "1" has probability `2^-logp` (Â§5.1.2.2).
    #[inline]
    pub fn bit_logp(&mut self, bit: bool, logp: u32) {
        let r = self.rng >> logp;
        if bit {
            self.val += self.rng - r;
            self.rng = r;
        } else {
            self.rng -= r;
        }
        self.normalize();
    }

    /// Symbol `s` of an inverse-CDF table with `ft = 1 << ftb` (Â§5.1.2.3).
    #[inline]
    pub fn icdf(&mut self, s: usize, icdf: &[u8], ftb: u32) {
        let r = self.rng >> ftb;
        if s > 0 {
            self.val += self.rng - r * u32::from(icdf[s - 1]);
            self.rng = r * (u32::from(icdf[s - 1]) - u32::from(icdf[s]));
        } else {
            self.rng -= r * u32::from(icdf[s]);
        }
        self.normalize();
    }

    /// A uniformly distributed integer `t` in `[0, ft)` (Â§5.1.4).
    pub fn uint(&mut self, t: u32, ft: u32) {
        debug_assert!(ft > 1 && t < ft);
        let ftm1 = ft - 1;
        let ftb = ilog(ftm1);
        if ftb > 8 {
            let sh = (ftb - 8) as u32;
            let ft1 = (ftm1 >> sh) + 1;
            let hi = t >> sh;
            self.encode(hi, hi + 1, ft1);
            self.bits(t & ((1 << sh) - 1), sh);
        } else {
            self.encode(t, t + 1, ft);
        }
    }

    /// `n` raw bits at the end of the frame (Â§5.1.3), `n <= 25`.
    pub fn bits(&mut self, value: u32, n: u32) {
        if n == 0 {
            return;
        }
        let mut window = self.end_window;
        let mut used = self.nend_bits;
        if used + n as i32 > 32 {
            loop {
                self.write_byte_at_end(window & 255);
                window >>= 8;
                used -= 8;
                if used < 8 {
                    break;
                }
            }
        }
        window |= value << used;
        used += n as i32;
        self.end_window = window;
        self.nend_bits = used;
        self.nbits_total += n as i32;
    }

    /// Bits used so far, rounded up.
    #[inline]
    pub fn tell(&self) -> i32 {
        self.nbits_total - ilog(self.rng)
    }

    /// Bits used so far in 1/8 bits.
    pub fn tell_frac(&self) -> i32 {
        tell_frac(self.nbits_total, self.rng)
    }

    /// The current range (the "final range" once the frame is done).
    pub fn range(&self) -> u32 {
        self.rng
    }

    /// The frame size in bytes.
    pub fn storage(&self) -> usize {
        self.buf.len()
    }

    /// Bytes written from the front so far (range-coded data only).
    pub fn range_bytes(&self) -> usize {
        self.offs
    }

    /// Whether the frame overflowed.
    pub fn error(&self) -> bool {
        self.error
    }

    /// Shrinks the frame to `size` bytes, moving any raw bits already
    /// written to the new end.
    pub fn shrink(&mut self, size: usize) {
        debug_assert!(self.offs + self.end_offs <= size);
        let old = self.buf.len();
        let tail: Vec<u8> = self.buf[old - self.end_offs..].to_vec();
        self.buf.truncate(size);
        let n = self.buf.len();
        self.buf[n - self.end_offs..].copy_from_slice(&tail);
    }

    /// Finalizes the frame (Â§5.1.5) and returns its bytes.
    pub fn finish(mut self) -> Vec<u8> {
        self.done();
        self.buf
    }

    /// The Â§5.1.5 termination.
    fn done(&mut self) {
        // The value in [val, val+rng) with the most trailing zero bits.
        let mut l = 32 - ilog(self.rng);
        let mut msk: u32 = 0x7FFF_FFFF >> l;
        let mut end = (self.val.wrapping_add(msk)) & !msk;
        if (end | msk) >= self.val.wrapping_add(self.rng) {
            l += 1;
            msk >>= 1;
            end = (self.val.wrapping_add(msk)) & !msk;
        }
        while l > 0 {
            self.carry_out(end >> 23);
            end = (end << 8) & 0x7FFF_FFFF;
            l -= 8;
        }
        if self.rem >= 0 || self.ext > 0 {
            self.carry_out(0);
        }
        // Flush the raw bits.
        let mut window = self.end_window;
        let mut used = self.nend_bits;
        while used >= 8 {
            self.write_byte_at_end(window & 255);
            window >>= 8;
            used -= 8;
        }
        if !self.error {
            let n = self.buf.len();
            for b in &mut self.buf[self.offs..n - self.end_offs] {
                *b = 0;
            }
            if used > 0 {
                if self.end_offs >= n {
                    self.error = true;
                } else {
                    let l = -l;
                    if self.offs + self.end_offs >= n && l < used {
                        window &= (1 << l) - 1;
                        self.error = true;
                    }
                    self.buf[n - self.end_offs - 1] |= window as u8;
                }
            }
        }
    }
}

/// What the CELT band and allocation code needs from either side of the
/// range coder: on the encoder each call writes the value it is given and
/// returns it; on the decoder each call ignores it and returns what it read.
pub(crate) trait Coder {
    /// True for the encoder.
    const ENCODE: bool;
    fn tell_frac(&self) -> i32;
    fn bit_logp(&mut self, value: bool, logp: u32) -> bool;
    fn uint(&mut self, value: u32, ft: u32) -> u32;
    fn bits(&mut self, value: u32, n: u32) -> u32;
    /// Decoder only: `decode(ft)`.
    fn decode_fs(&mut self, ft: u32) -> u32;
    /// Encoder: encodes `(fl, fh, ft)`; decoder: updates with it.
    fn code(&mut self, fl: u32, fh: u32, ft: u32);
}

impl Coder for RangeDecoder<'_> {
    const ENCODE: bool = false;
    fn tell_frac(&self) -> i32 {
        RangeDecoder::tell_frac(self)
    }
    fn bit_logp(&mut self, _: bool, logp: u32) -> bool {
        RangeDecoder::bit_logp(self, logp)
    }
    fn uint(&mut self, _: u32, ft: u32) -> u32 {
        RangeDecoder::uint(self, ft)
    }
    fn bits(&mut self, _: u32, n: u32) -> u32 {
        RangeDecoder::bits(self, n)
    }
    fn decode_fs(&mut self, ft: u32) -> u32 {
        self.decode(ft)
    }
    fn code(&mut self, fl: u32, fh: u32, ft: u32) {
        self.update(fl, fh, ft)
    }
}

impl Coder for RangeEncoder {
    const ENCODE: bool = true;
    fn tell_frac(&self) -> i32 {
        RangeEncoder::tell_frac(self)
    }
    fn bit_logp(&mut self, value: bool, logp: u32) -> bool {
        RangeEncoder::bit_logp(self, value, logp);
        value
    }
    fn uint(&mut self, value: u32, ft: u32) -> u32 {
        RangeEncoder::uint(self, value, ft);
        value
    }
    fn bits(&mut self, value: u32, n: u32) -> u32 {
        RangeEncoder::bits(self, value, n);
        value
    }
    fn decode_fs(&mut self, _: u32) -> u32 {
        unreachable!("decode_fs on an encoder")
    }
    fn code(&mut self, fl: u32, fh: u32, ft: u32) {
        self.encode(fl, fh, ft)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tiny deterministic generator for the tests.
    struct Lcg(u32);
    impl Lcg {
        fn next(&mut self) -> u32 {
            self.0 = self.0.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            self.0
        }
    }

    #[derive(Debug, Clone, Copy)]
    enum Op {
        Sym(u32, u32),
        Bin(u32, u32),
        Logp(bool, u32),
        Icdf(usize),
        Uint(u32, u32),
        Bits(u32, u32),
    }

    const ICDF: [u8; 4] = [200, 120, 30, 0];

    fn ops(seed: u32, n: usize) -> Vec<Op> {
        let mut g = Lcg(seed);
        (0..n)
            .map(|_| match g.next() >> 29 {
                0 => {
                    let ft = 2 + g.next() % 60000;
                    Op::Sym(g.next() % ft, ft)
                }
                1 => {
                    let bits = 1 + g.next() % 15;
                    Op::Bin(g.next() % (1 << bits), bits)
                }
                2 | 3 => Op::Logp(g.next().is_multiple_of(4), 1 + g.next() % 15),
                4 => Op::Icdf((g.next() % 4) as usize),
                5 => {
                    let ft = 2 + (g.next() >> (g.next() % 31));
                    let ft = ft.max(2);
                    Op::Uint(g.next() % ft, ft)
                }
                _ => {
                    let n = 1 + g.next() % 25;
                    Op::Bits(g.next() & ((1 << n) - 1), n)
                }
            })
            .collect()
    }

    /// Everything encoded decodes to itself, and the encoder's range and
    /// `tell_frac` agree with the decoder's after every symbol (Â§5.1: "the
    /// value of rng in the encoder should exactly match").
    #[test]
    fn round_trip_and_matching_state() {
        for seed in 0..200 {
            let ops = ops(seed, 1 + (seed as usize * 7) % 300);
            let mut enc = RangeEncoder::new(4000);
            let mut trace = Vec::new();
            for &op in &ops {
                match op {
                    Op::Sym(s, ft) => enc.encode(s, s + 1, ft),
                    Op::Bin(s, b) => enc.encode_bin(s, s + 1, b),
                    Op::Logp(v, l) => enc.bit_logp(v, l),
                    Op::Icdf(s) => enc.icdf(s, &ICDF, 8),
                    Op::Uint(t, ft) => enc.uint(t, ft),
                    Op::Bits(v, n) => enc.bits(v, n),
                }
                trace.push((enc.range(), enc.tell_frac()));
            }
            let final_range = enc.range();
            let bytes = enc.finish();
            let mut dec = RangeDecoder::new(&bytes);
            for (i, &op) in ops.iter().enumerate() {
                match op {
                    Op::Sym(s, ft) => {
                        let f = dec.decode(ft);
                        assert_eq!(f, s, "seed {seed} op {i}");
                        dec.update(s, s + 1, ft);
                    }
                    Op::Bin(s, b) => {
                        let f = dec.decode_bin(b);
                        assert_eq!(f, s);
                        dec.update(s, s + 1, 1 << b);
                    }
                    Op::Logp(v, l) => assert_eq!(dec.bit_logp(l), v),
                    Op::Icdf(s) => assert_eq!(dec.icdf(&ICDF, 8), s),
                    Op::Uint(t, ft) => assert_eq!(dec.uint(ft), t),
                    Op::Bits(v, n) => assert_eq!(dec.bits(n), v),
                }
                assert_eq!((dec.range(), dec.tell_frac()), trace[i], "seed {seed} op {i}");
            }
            assert_eq!(dec.range(), final_range);
            assert!(!dec.error());
        }
    }

    /// Â§4.1.6.1: a fresh decoder reports one bit used, and `tell` is the
    /// ceiling of `tell_frac / 8`.
    #[test]
    fn tell_starts_at_one_bit() {
        let dec = RangeDecoder::new(&[0x55, 0xAA, 0x12]);
        assert_eq!(dec.tell(), 1);
        let enc = RangeEncoder::new(10);
        assert_eq!(enc.tell(), 1);
        let mut d = RangeDecoder::new(&[0x37, 0x99, 0x01, 0xFE, 0x44]);
        for _ in 0..6 {
            d.bit_logp(3);
            assert_eq!(d.tell(), (d.tell_frac() + 7) >> 3);
        }
    }

    /// The smallest frame that holds the symbols: every encoded frame of
    /// exactly-sized storage still decodes.
    #[test]
    fn tight_frames_decode() {
        for seed in 0..100 {
            let ops: Vec<bool> = (0..40).map(|i| (seed * 31 + i * 7) % 5 == 0).collect();
            let mut enc = RangeEncoder::new(64);
            for &b in &ops {
                enc.bit_logp(b, 2);
            }
            let used = (enc.tell() + 7) as usize / 8;
            enc.shrink(used);
            let bytes = enc.finish();
            assert_eq!(bytes.len(), used);
            let mut dec = RangeDecoder::new(&bytes);
            for &b in &ops {
                assert_eq!(dec.bit_logp(2), b);
            }
        }
    }
}
