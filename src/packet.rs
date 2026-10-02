//! Opus packet framing (RFC 6716 §3): the TOC byte, the four frame-count
//! codes, frame lengths and Opus padding — parsing with every rule [R1]–[R7]
//! enforced, and building packets that obey them.

use crate::error::{Result, invalid};

/// The operating mode of a frame (RFC 6716 §3.1, Table 2).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Mode {
    /// The LP (SILK) layer only.
    Silk,
    /// SILK below 8 kHz and CELT above.
    Hybrid,
    /// The MDCT (CELT) layer only.
    Celt,
}

/// Audio bandwidth (RFC 6716 §2.1.3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Bandwidth {
    /// Narrowband, 4 kHz.
    Narrow,
    /// Medium-band, 6 kHz.
    Medium,
    /// Wideband, 8 kHz.
    Wide,
    /// Super-wideband, 12 kHz.
    SuperWide,
    /// Fullband, 20 kHz.
    Full,
}

impl Bandwidth {
    /// The audio bandwidth in Hz.
    pub fn hz(self) -> u32 {
        match self {
            Bandwidth::Narrow => 4000,
            Bandwidth::Medium => 6000,
            Bandwidth::Wide => 8000,
            Bandwidth::SuperWide => 12000,
            Bandwidth::Full => 20000,
        }
    }

    /// The number of CELT bands coded at this bandwidth (RFC 6716 Table 55).
    pub(crate) fn celt_end_band(self) -> usize {
        match self {
            Bandwidth::Narrow => 13,
            Bandwidth::Medium | Bandwidth::Wide => 17,
            Bandwidth::SuperWide => 19,
            Bandwidth::Full => 21,
        }
    }
}

/// What a TOC byte says (RFC 6716 §3.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Toc {
    /// The configuration number, 0–31.
    pub config: u8,
    /// The `s` bit.
    pub stereo: bool,
    /// The frame count code `c`, 0–3.
    pub code: u8,
}

impl Toc {
    /// Splits a TOC byte.
    pub fn from_byte(b: u8) -> Self {
        Self { config: b >> 3, stereo: b & 4 != 0, code: b & 3 }
    }

    /// The TOC byte.
    pub fn to_byte(self) -> u8 {
        (self.config << 3) | (u8::from(self.stereo) << 2) | (self.code & 3)
    }

    /// The configuration number for a mode, bandwidth and frame size (in
    /// 48 kHz samples), if Table 2 has one.
    pub fn config_for(mode: Mode, bw: Bandwidth, frame_size: usize) -> Option<u8> {
        let size_index = |sizes: [usize; 4]| sizes.iter().position(|&s| s == frame_size);
        match mode {
            Mode::Silk => {
                let base = match bw {
                    Bandwidth::Narrow => 0,
                    Bandwidth::Medium => 4,
                    Bandwidth::Wide => 8,
                    _ => return None,
                };
                Some(base + size_index([480, 960, 1920, 2880])? as u8)
            }
            Mode::Hybrid => {
                let base = match bw {
                    Bandwidth::SuperWide => 12,
                    Bandwidth::Full => 14,
                    _ => return None,
                };
                Some(base + [480, 960].iter().position(|&s| s == frame_size)? as u8)
            }
            Mode::Celt => {
                let base = match bw {
                    Bandwidth::Narrow => 16,
                    Bandwidth::Medium => return None,
                    Bandwidth::Wide => 20,
                    Bandwidth::SuperWide => 24,
                    Bandwidth::Full => 28,
                };
                Some(base + size_index([120, 240, 480, 960])? as u8)
            }
        }
    }

    /// The operating mode.
    pub fn mode(self) -> Mode {
        match self.config {
            0..=11 => Mode::Silk,
            12..=15 => Mode::Hybrid,
            _ => Mode::Celt,
        }
    }

    /// The audio bandwidth.
    pub fn bandwidth(self) -> Bandwidth {
        match self.config {
            0..=3 => Bandwidth::Narrow,
            4..=7 => Bandwidth::Medium,
            8..=11 => Bandwidth::Wide,
            12 | 13 => Bandwidth::SuperWide,
            14 | 15 => Bandwidth::Full,
            16..=19 => Bandwidth::Narrow,
            20..=23 => Bandwidth::Wide,
            24..=27 => Bandwidth::SuperWide,
            _ => Bandwidth::Full,
        }
    }

    /// The duration of one frame in samples at 48 kHz.
    pub fn frame_size(self) -> usize {
        match self.config {
            0..=11 => [480, 960, 1920, 2880][usize::from(self.config & 3)],
            12..=15 => [480, 960][usize::from(self.config & 1)],
            _ => [120, 240, 480, 960][usize::from(self.config & 3)],
        }
    }

    /// The number of channels the packet codes.
    pub fn channels(self) -> usize {
        if self.stereo { 2 } else { 1 }
    }
}

/// A parsed packet: its TOC and its frames, borrowed from the input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Packet<'a> {
    /// The TOC byte.
    pub toc: Toc,
    /// The compressed frames, in order; a frame may be empty (DTX or lost,
    /// §3.2.1).
    pub frames: Vec<&'a [u8]>,
    /// Opus padding bytes (code 3 only).
    pub padding: usize,
}

impl Packet<'_> {
    /// Samples per channel at 48 kHz the packet decodes to.
    pub fn samples(&self) -> usize {
        self.frames.len() * self.toc.frame_size()
    }
}

/// Reads a §3.2.1 frame length; returns `(length, bytes used)`.
fn frame_length(data: &[u8]) -> Result<(usize, usize)> {
    match data {
        [] => Err(invalid("frame length missing")),
        [b, ..] if *b < 252 => Ok((usize::from(*b), 1)),
        [_] => Err(invalid("two-byte frame length cut short")),
        [b0, b1, ..] => Ok((usize::from(*b1) * 4 + usize::from(*b0), 2)),
    }
}

/// The longest frame (§3.2.1, [R2]).
pub const MAX_FRAME_BYTES: usize = 1275;

/// Parses a packet (RFC 6716 §3.2), enforcing [R1]–[R7].
pub fn parse(data: &[u8]) -> Result<Packet<'_>> {
    let Some((&first, mut rest)) = data.split_first() else {
        return Err(invalid("empty packet [R1]"));
    };
    let toc = Toc::from_byte(first);
    let mut frames = Vec::new();
    let mut padding = 0;
    match toc.code {
        0 => {
            if rest.len() > MAX_FRAME_BYTES {
                return Err(invalid("frame longer than 1275 bytes [R2]"));
            }
            frames.push(rest);
        }
        1 => {
            if rest.len() % 2 != 0 {
                return Err(invalid("code 1 packet with an odd payload [R3]"));
            }
            let half = rest.len() / 2;
            if half > MAX_FRAME_BYTES {
                return Err(invalid("frame longer than 1275 bytes [R2]"));
            }
            frames.push(&rest[..half]);
            frames.push(&rest[half..]);
        }
        2 => {
            let (n1, used) = frame_length(rest).map_err(|_| invalid("code 2 length [R4]"))?;
            rest = &rest[used..];
            if n1 > rest.len() {
                return Err(invalid("code 2 first frame longer than the packet [R4]"));
            }
            if rest.len() - n1 > MAX_FRAME_BYTES {
                return Err(invalid("frame longer than 1275 bytes [R2]"));
            }
            frames.push(&rest[..n1]);
            frames.push(&rest[n1..]);
        }
        _ => {
            let Some((&fc, r)) = rest.split_first() else {
                return Err(invalid("code 3 packet without a frame count [R6]"));
            };
            rest = r;
            let vbr = fc & 0x80 != 0;
            let has_padding = fc & 0x40 != 0;
            let m = usize::from(fc & 0x3F);
            if m == 0 {
                return Err(invalid("code 3 packet with zero frames [R5]"));
            }
            if m * toc.frame_size() > 5760 {
                return Err(invalid("packet longer than 120 ms [R5]"));
            }
            if has_padding {
                // RFC 8251 §4: count the padding down from the remaining
                // length, never up into an overflow.
                let mut remaining = rest.len() as isize;
                let mut i = 0;
                loop {
                    if remaining <= 0 {
                        return Err(invalid("padding length cut short [R6]"));
                    }
                    let p = rest[i];
                    i += 1;
                    remaining -= 1;
                    let add = if p == 255 { 254 } else { usize::from(p) };
                    remaining -= add as isize;
                    padding += add;
                    if p != 255 {
                        break;
                    }
                }
                if remaining < 0 {
                    return Err(invalid("padding longer than the packet [R6]"));
                }
                rest = &rest[i..i + remaining as usize];
            }
            if vbr {
                let mut lens = Vec::with_capacity(m);
                for _ in 0..m - 1 {
                    let (n, used) = frame_length(rest).map_err(|_| invalid("code 3 frame lengths [R7]"))?;
                    rest = &rest[used..];
                    lens.push(n);
                }
                let total: usize = lens.iter().sum();
                if total > rest.len() {
                    return Err(invalid("code 3 frame lengths exceed the packet [R7]"));
                }
                let last = rest.len() - total;
                if last > MAX_FRAME_BYTES {
                    return Err(invalid("frame longer than 1275 bytes [R2]"));
                }
                lens.push(last);
                for n in lens {
                    frames.push(&rest[..n]);
                    rest = &rest[n..];
                }
            } else {
                if rest.len() % m != 0 {
                    return Err(invalid("CBR code 3 payload not a multiple of the frame count [R6]"));
                }
                let n = rest.len() / m;
                if n > MAX_FRAME_BYTES {
                    return Err(invalid("frame longer than 1275 bytes [R2]"));
                }
                for k in 0..m {
                    frames.push(&rest[k * n..(k + 1) * n]);
                }
            }
        }
    }
    Ok(Packet { toc, frames, padding })
}

/// The number of samples per channel at 48 kHz in a packet, without
/// decoding it.
pub fn packet_samples(data: &[u8]) -> Result<usize> {
    Ok(parse(data)?.samples())
}

fn push_length(out: &mut Vec<u8>, n: usize) {
    if n < 252 {
        out.push(n as u8);
    } else {
        let b0 = 252 + (n & 3);
        out.push(b0 as u8);
        out.push(((n - b0) / 4) as u8);
    }
}

/// Builds a packet from frames that share `toc` (whose `code` is ignored):
/// code 0 for one frame, code 1 or 2 for two, code 3 otherwise — or code 3
/// whenever `pad_to` asks for a total size the plain codes cannot reach.
/// Every frame must be at most 1275 bytes and the packet at most 120 ms.
pub fn build(toc: Toc, frames: &[&[u8]], pad_to: Option<usize>) -> Result<Vec<u8>> {
    if frames.is_empty() {
        return Err(invalid("a packet needs at least one frame"));
    }
    if frames.iter().any(|f| f.len() > MAX_FRAME_BYTES) {
        return Err(invalid("frame longer than 1275 bytes"));
    }
    if frames.len() * toc.frame_size() > 5760 || frames.len() > 48 {
        return Err(invalid("packet longer than 120 ms"));
    }
    let all_equal = frames.iter().all(|f| f.len() == frames[0].len());
    let mut out = Vec::new();
    let plain = match frames.len() {
        1 => {
            out.push(Toc { code: 0, ..toc }.to_byte());
            out.extend_from_slice(frames[0]);
            true
        }
        2 if all_equal => {
            out.push(Toc { code: 1, ..toc }.to_byte());
            out.extend_from_slice(frames[0]);
            out.extend_from_slice(frames[1]);
            true
        }
        2 => {
            out.push(Toc { code: 2, ..toc }.to_byte());
            push_length(&mut out, frames[0].len());
            out.extend_from_slice(frames[0]);
            out.extend_from_slice(frames[1]);
            true
        }
        _ => false,
    };
    let target = pad_to.unwrap_or(0);
    if plain && out.len() >= target {
        return Ok(out);
    }
    // Code 3, VBR unless the frames are all the same size.
    let body: usize = frames.iter().map(|f| f.len()).sum();
    let mut lengths = Vec::new();
    if !all_equal {
        for f in &frames[..frames.len() - 1] {
            push_length(&mut lengths, f.len());
        }
    }
    let unpadded = 2 + lengths.len() + body;
    let mut pad_header = Vec::new();
    let mut pad_bytes = 0;
    if target > unpadded {
        // P = (header bytes) + (padding bytes) = target - unpadded.
        let p = target - unpadded;
        // Each 255 header byte adds 254 bytes of padding plus itself.
        let mut left = p;
        while left > 255 {
            pad_header.push(255u8);
            left -= 255;
        }
        // The last header byte v adds v padding bytes plus itself.
        pad_header.push((left - 1) as u8);
        pad_bytes = p - pad_header.len();
    }
    out.clear();
    out.push(Toc { code: 3, ..toc }.to_byte());
    let fc = (u8::from(!all_equal) << 7) | (u8::from(!pad_header.is_empty()) << 6) | frames.len() as u8;
    out.push(fc);
    out.extend_from_slice(&pad_header);
    out.extend_from_slice(&lengths);
    for f in frames {
        out.extend_from_slice(f);
    }
    out.resize(out.len() + pad_bytes, 0);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn toc_table_2() {
        let t = Toc::from_byte(1 << 3);
        assert_eq!((t.mode(), t.bandwidth(), t.frame_size()), (Mode::Silk, Bandwidth::Narrow, 960));
        let t = Toc::from_byte(29 << 3 | 1);
        assert_eq!((t.mode(), t.bandwidth(), t.frame_size(), t.code), (Mode::Celt, Bandwidth::Full, 240, 1));
        let t = Toc::from_byte(15 << 3);
        assert_eq!((t.mode(), t.bandwidth(), t.frame_size()), (Mode::Hybrid, Bandwidth::Full, 960));
        for config in 0..32u8 {
            let t = Toc { config, stereo: false, code: 0 };
            assert_eq!(Toc::config_for(t.mode(), t.bandwidth(), t.frame_size()), Some(config));
        }
    }

    #[test]
    fn rules_r1_to_r7() {
        assert!(parse(&[]).is_err(), "R1");
        assert!(parse(&[1]).is_ok());
        assert!(parse(&[1 | 1, 0, 0, 0]).is_err(), "R3: odd payload");
        assert!(parse(&[1 | 1, 0, 0]).unwrap().frames.len() == 2);
        assert!(parse(&[2]).is_err(), "R4: no length");
        assert!(parse(&[2, 253]).is_err(), "R4: two-byte length cut");
        assert!(parse(&[2, 5, 0]).is_err(), "R4: length past end");
        assert_eq!(parse(&[2, 0]).unwrap().frames, vec![&[][..], &[][..]]);
        assert!(parse(&[3]).is_err(), "R6: no count byte");
        assert!(parse(&[3, 0]).is_err(), "R5: zero frames");
        assert!(parse(&[(3 << 3) | 3, 3]).is_err(), "R5: 180 ms");
        assert!(parse(&[3, 0x41]).is_err(), "R6: padding byte missing");
        assert!(parse(&[3, 0x41, 2, 0]).is_err(), "R6: padding past end");
        assert_eq!(parse(&[3, 0x41, 1, 0]).unwrap().frames, vec![&[][..]]);
        assert!(parse(&[3, 2, 1]).is_err(), "R6: not a multiple");
        assert!(parse(&[3, 0x82, 5, 1]).is_err(), "R7");
        let p = parse(&[3, 0x83, 1, 2, 9, 8, 8, 7]).unwrap();
        assert_eq!(p.frames, vec![&[9][..], &[8, 8][..], &[7][..]]);
        let big = vec![0u8; 1277];
        assert!(parse(&big).is_err(), "R2");
    }

    #[test]
    fn built_packets_parse_back() {
        let toc = Toc { config: 31, stereo: true, code: 0 };
        let a = vec![1u8; 300];
        let b = vec![2u8; 17];
        let c = vec![3u8; 600];
        let cases: Vec<Vec<&[u8]>> = vec![
            vec![&a],
            vec![&a, &a],
            vec![&a, &b],
            vec![&b, &c, &a],
            vec![&b, &b, &b, &b],
            vec![&[], &b],
        ];
        for frames in &cases {
            for pad in [None, Some(700), Some(1000), Some(2000)] {
                let p = build(toc, frames, pad).unwrap();
                if let Some(t) = pad {
                    assert!(p.len() >= t.min(p.len()));
                    if t >= 2 + frames.iter().map(|f| f.len() + 2).sum::<usize>() {
                        assert_eq!(p.len(), t, "{frames:?}");
                    }
                }
                let parsed = parse(&p).unwrap();
                assert_eq!(&parsed.frames, frames);
            }
        }
    }
}
