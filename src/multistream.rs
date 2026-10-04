//! Multistream Opus (RFC 7845 §5.1.1, RFC 6716 Appendix B): the `OpusHead`
//! identification header, channel mapping families 0, 1 and 255, and a
//! decoder that splits a multistream packet into its streams.

use crate::decoder::Decoder;
use crate::encoder::{Encoder, EncoderConfig};
use crate::error::{Error, Result};
use crate::packet::{self, Bandwidth};

/// The contents of an `OpusHead` identification header (RFC 7845 §5.1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OpusHead {
    /// The version (1 for RFC 7845; any value up to 15 is compatible).
    pub version: u8,
    /// Output channels.
    pub channels: u8,
    /// Samples at 48 kHz to discard from the start of the decoded output.
    pub pre_skip: u16,
    /// The sample rate of the original input, in Hz (0: unspecified).
    pub input_sample_rate: u32,
    /// Output gain, Q7.8 dB.
    pub output_gain: i16,
    /// Channel mapping family.
    pub family: u8,
    /// Number of Opus streams.
    pub streams: u8,
    /// Number of those that are coupled (stereo).
    pub coupled: u8,
    /// For each output channel, the decoded channel it takes (255:
    /// silence).
    pub mapping: Vec<u8>,
}

/// The Vorbis-order layouts of mapping family 1 (RFC 7845 §5.1.1.2) as
/// coded by [`MultistreamEncoder`](crate::MultistreamEncoder): `(streams,
/// coupled, mapping)` for 1 to 8 channels. Pairs that belong together
/// (front, rear, side) share a coupled stream; centre and LFE are mono.
pub fn family1_layout(channels: u8) -> Option<(u8, u8, &'static [u8])> {
    Some(match channels {
        1 => (1, 0, &[0]),
        2 => (1, 1, &[0, 1]),
        // L C R
        3 => (2, 1, &[0, 2, 1]),
        // FL FR RL RR
        4 => (2, 2, &[0, 1, 2, 3]),
        // FL C FR RL RR
        5 => (3, 2, &[0, 4, 1, 2, 3]),
        // FL C FR RL RR LFE
        6 => (4, 2, &[0, 4, 1, 2, 3, 5]),
        // FL C FR SL SR RC LFE
        7 => (4, 3, &[0, 4, 1, 2, 3, 5, 6]),
        // FL C FR SL SR RL RR LFE
        8 => (5, 3, &[0, 6, 1, 2, 3, 4, 5, 7]),
        _ => return None,
    })
}

impl OpusHead {
    /// A head for `channels` with family 0 (1–2 channels) or 1 (3–8).
    pub fn new(channels: u8, pre_skip: u16, input_sample_rate: u32) -> Result<Self> {
        let (streams, coupled, mapping) = family1_layout(channels).ok_or_else(|| {
            Error::Config(format!(
                "{channels} channels: families 0 and 1 cover 1 to 8"
            ))
        })?;
        Ok(Self {
            version: 1,
            channels,
            pre_skip,
            input_sample_rate,
            output_gain: 0,
            family: if channels <= 2 { 0 } else { 1 },
            streams,
            coupled,
            mapping: mapping.to_vec(),
        })
    }

    /// Parses an identification header, with or without the `OpusHead`
    /// magic (MP4 `dOps` and Matroska `CodecPrivate` forms differ in it).
    pub fn parse(data: &[u8]) -> Result<Self> {
        let bad = |why: String| Error::InvalidHead(why);
        let body = data.strip_prefix(b"OpusHead").unwrap_or(data);
        if body.len() < 11 {
            return Err(bad(format!("{} bytes, need 11", body.len())));
        }
        let version = body[0];
        if version >= 16 {
            return Err(bad(format!("version {version} is incompatible")));
        }
        let channels = body[1];
        if channels == 0 {
            return Err(bad("zero channels".into()));
        }
        let pre_skip = u16::from_le_bytes([body[2], body[3]]);
        let input_sample_rate = u32::from_le_bytes([body[4], body[5], body[6], body[7]]);
        let output_gain = i16::from_le_bytes([body[8], body[9]]);
        let family = body[10];
        if family == 0 {
            if channels > 2 {
                return Err(bad(format!("family 0 with {channels} channels")));
            }
            return Ok(Self {
                version,
                channels,
                pre_skip,
                input_sample_rate,
                output_gain,
                family,
                streams: 1,
                coupled: channels - 1,
                mapping: (0..channels).collect(),
            });
        }
        let need = 13 + usize::from(channels);
        if body.len() < need {
            return Err(bad(format!(
                "family {family} needs {need} bytes, has {}",
                body.len()
            )));
        }
        let (streams, coupled) = (body[11], body[12]);
        if streams == 0 || coupled > streams || usize::from(streams) + usize::from(coupled) > 255 {
            return Err(bad(format!("{streams} streams, {coupled} coupled")));
        }
        let mapping = body[13..need].to_vec();
        let decoded = streams + coupled;
        if let Some(&m) = mapping.iter().find(|&&m| m != 255 && m >= decoded) {
            return Err(bad(format!(
                "mapping index {m} with {decoded} decoded channels"
            )));
        }
        if family == 1 && channels > 8 {
            return Err(bad(format!("family 1 with {channels} channels")));
        }
        Ok(Self {
            version,
            channels,
            pre_skip,
            input_sample_rate,
            output_gain,
            family,
            streams,
            coupled,
            mapping,
        })
    }

    /// The header with the `OpusHead` magic (the Ogg packet).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut v = b"OpusHead".to_vec();
        v.extend_from_slice(&self.body());
        v
    }

    /// The header without the magic (the `dOps` / `CodecPrivate` body
    /// layout, little-endian as in Ogg).
    pub fn body(&self) -> Vec<u8> {
        let mut v = vec![self.version, self.channels];
        v.extend_from_slice(&self.pre_skip.to_le_bytes());
        v.extend_from_slice(&self.input_sample_rate.to_le_bytes());
        v.extend_from_slice(&self.output_gain.to_le_bytes());
        v.push(self.family);
        if self.family != 0 {
            v.push(self.streams);
            v.push(self.coupled);
            v.extend_from_slice(&self.mapping);
        }
        v
    }

    /// The linear factor of the output gain.
    pub fn gain_factor(&self) -> f32 {
        10f32.powf(f32::from(self.output_gain) / (20.0 * 256.0))
    }
}

/// A decoder for multistream packets: one [`Decoder`] per stream, the
/// first `coupled` stereo, their channels routed to the outputs by the
/// mapping.
pub struct MultistreamDecoder {
    decoders: Vec<Decoder>,
    coupled: usize,
    mapping: Vec<u8>,
    gain: f32,
}

impl MultistreamDecoder {
    /// A decoder for `streams` streams of which `coupled` are stereo, with
    /// one `mapping` entry per output channel.
    pub fn new(sample_rate: u32, streams: u8, coupled: u8, mapping: &[u8]) -> Result<Self> {
        if streams == 0 || coupled > streams {
            return Err(Error::Config(format!(
                "{streams} streams, {coupled} coupled"
            )));
        }
        let decoded = streams + coupled;
        if mapping.is_empty() || mapping.iter().any(|&m| m != 255 && m >= decoded) {
            return Err(Error::Config(
                "mapping refers to a channel that is not decoded".into(),
            ));
        }
        let decoders = (0..streams)
            .map(|s| Decoder::new(sample_rate, if s < coupled { 2 } else { 1 }))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            decoders,
            coupled: usize::from(coupled),
            mapping: mapping.to_vec(),
            gain: 1.0,
        })
    }

    /// A decoder for the stream an `OpusHead` describes, applying its
    /// output gain.
    pub fn from_head(head: &OpusHead, sample_rate: u32) -> Result<Self> {
        let mut d = Self::new(sample_rate, head.streams, head.coupled, &head.mapping)?;
        d.gain = head.gain_factor();
        Ok(d)
    }

    /// Output channels.
    pub fn channels(&self) -> usize {
        self.mapping.len()
    }

    /// See [`Decoder::set_phase_inversion_disabled`].
    pub fn set_phase_inversion_disabled(&mut self, disabled: bool) {
        for d in &mut self.decoders {
            d.set_phase_inversion_disabled(disabled);
        }
    }

    /// Resets every stream's decoder.
    pub fn reset(&mut self) {
        for d in &mut self.decoders {
            d.reset();
        }
    }

    /// The XOR of the streams' final ranges.
    pub fn final_range(&self) -> u32 {
        self.decoders.iter().fold(0, |a, d| a ^ d.final_range())
    }

    /// Decodes one multistream packet (`None`: lost) into interleaved
    /// samples, one per output channel.
    pub fn decode(&mut self, packet: Option<&[u8]>) -> Result<Vec<f32>> {
        let n = self.decoders.len();
        let mut outs: Vec<Vec<f32>> = Vec::with_capacity(n);
        match packet {
            None | Some([]) => {
                for d in &mut self.decoders {
                    outs.push(d.decode(None)?);
                }
            }
            Some(mut data) => {
                for s in 0..n {
                    let piece: &[u8] = if s + 1 < n {
                        let (_, used) = packet::parse_self_delimited(data)?;
                        let (head, rest) = data.split_at(used);
                        data = rest;
                        head
                    } else {
                        data
                    };
                    let pcm = if s + 1 < n {
                        // Rebuild the regular framing for the stream decoder.
                        let (p, _) = packet::parse_self_delimited(piece)?;
                        let regular = packet::build(p.toc, &p.frames, None)?;
                        self.decoders[s].decode(Some(&regular))?
                    } else {
                        self.decoders[s].decode(Some(piece))?
                    };
                    outs.push(pcm);
                }
            }
        }
        let frames = outs[0].len() / self.decoders[0].channels();
        for (s, o) in outs.iter().enumerate() {
            if o.len() / self.decoders[s].channels() != frames {
                return Err(Error::InvalidPacket(
                    "streams of different durations".into(),
                ));
            }
        }
        let c = self.mapping.len();
        let mut out = vec![0.0f32; frames * c];
        for (oc, &m) in self.mapping.iter().enumerate() {
            if m == 255 {
                continue;
            }
            let m = usize::from(m);
            let (s, ch, nch) = if m < 2 * self.coupled {
                (m / 2, m % 2, 2)
            } else {
                (m - self.coupled, 0, 1)
            };
            for i in 0..frames {
                out[i * c + oc] = outs[s][i * nch + ch] * self.gain;
            }
        }
        Ok(out)
    }
}

/// An encoder for 1 to 8 channels: one [`Encoder`] per stream of the
/// family 0 (mono, stereo) or family 1 (surround, Vorbis order) layout of
/// [`family1_layout`], producing multistream packets (RFC 7845 §5.1.1:
/// every stream but the last self-delimited, RFC 6716 Appendix B).
pub struct MultistreamEncoder {
    encoders: Vec<Encoder>,
    coupled: usize,
    mapping: Vec<u8>,
    head: OpusHead,
}

impl MultistreamEncoder {
    /// An encoder for `cfg.channels` (1–8) input channels in Vorbis order;
    /// `cfg.bitrate` is the total for all streams.
    pub fn new(cfg: EncoderConfig) -> Result<Self> {
        let ch =
            u8::try_from(cfg.channels).map_err(|_| Error::Config("too many channels".into()))?;
        let (streams, coupled, mapping) = family1_layout(ch).ok_or_else(|| {
            Error::Config(format!("{ch} channels: families 0 and 1 cover 1 to 8"))
        })?;
        // Share the rate: a coupled stream counts 1.5, a mono one 1, the
        // LFE (the last channel of 5.1 and 7.1) 0.25.
        let lfe = if ch == 6 || ch == 8 {
            Some(usize::from(mapping[usize::from(ch) - 1]) - usize::from(coupled))
        } else {
            None
        };
        let weight = |s: usize| -> f64 {
            if s < usize::from(coupled) {
                1.5
            } else if Some(s) == lfe {
                0.25
            } else {
                1.0
            }
        };
        let total: f64 = (0..usize::from(streams)).map(weight).sum();
        let mut encoders = Vec::with_capacity(usize::from(streams));
        for s in 0..usize::from(streams) {
            let channels = if s < usize::from(coupled) { 2 } else { 1 };
            let bitrate =
                ((f64::from(cfg.bitrate) * weight(s) / total) as u32).clamp(6000, 510_000);
            let mut c = EncoderConfig {
                channels,
                bitrate,
                ..cfg
            };
            if Some(s) == lfe {
                c.max_bandwidth = Some(Bandwidth::Narrow);
                c.mode = Some(crate::packet::Mode::Celt);
            }
            encoders.push(Encoder::new(c)?);
        }
        let pre_skip = encoders[0].lookahead() as u16;
        let mut head = OpusHead::new(ch, pre_skip, cfg.sample_rate)?;
        head.pre_skip = pre_skip;
        Ok(Self {
            encoders,
            coupled: usize::from(coupled),
            mapping: mapping.to_vec(),
            head,
        })
    }

    /// The identification header describing the stream.
    pub fn head(&self) -> &OpusHead {
        &self.head
    }

    /// Samples at 48 kHz of encoder delay (the head's pre-skip).
    pub fn lookahead(&self) -> usize {
        self.encoders[0].lookahead()
    }

    /// Input samples per channel that make one packet.
    pub fn frame_samples(&self) -> usize {
        self.encoders[0].frame_samples()
    }

    /// The XOR of the streams' final ranges.
    pub fn final_range(&self) -> u32 {
        self.encoders.iter().fold(0, |a, e| a ^ e.final_range())
    }

    /// Encodes one packet from interleaved input in Vorbis channel order.
    pub fn encode(&mut self, pcm: &[f32]) -> Result<Vec<u8>> {
        let c = self.mapping.len();
        let n = self.frame_samples();
        if pcm.len() != n * c {
            return Err(Error::BadArgument(format!(
                "{} samples given, a packet takes {}",
                pcm.len(),
                n * c
            )));
        }
        let mut inputs: Vec<Vec<f32>> = (0..self.encoders.len())
            .map(|s| vec![0.0; n * if s < self.coupled { 2 } else { 1 }])
            .collect();
        for (oc, &m) in self.mapping.iter().enumerate() {
            let m = usize::from(m);
            let (s, ch, nch) = if m < 2 * self.coupled {
                (m / 2, m % 2, 2)
            } else {
                (m - self.coupled, 0, 1)
            };
            for i in 0..n {
                inputs[s][i * nch + ch] = pcm[i * c + oc];
            }
        }
        let last = self.encoders.len() - 1;
        let mut out = Vec::new();
        for (s, e) in self.encoders.iter_mut().enumerate() {
            let p = e.encode(&inputs[s])?;
            if s < last {
                out.extend_from_slice(&packet::to_self_delimited(&p)?);
            } else {
                out.extend_from_slice(&p);
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn heads_round_trip() {
        let stereo = [1, 2, 0x38, 0x01, 0x80, 0xBB, 0, 0, 0, 0, 0];
        let h = OpusHead::parse(&stereo).unwrap();
        assert_eq!(
            (h.channels, h.pre_skip, h.family, h.streams, h.coupled),
            (2, 312, 0, 1, 1)
        );
        assert_eq!(h.body(), stereo);
        let mut surround = vec![1, 6, 0x38, 0x01, 0x80, 0xBB, 0, 0, 0, 0, 1, 4, 2];
        surround.extend_from_slice(&[0, 4, 1, 2, 3, 5]);
        let h = OpusHead::parse(&surround).unwrap();
        assert_eq!((h.channels, h.family, h.streams, h.coupled), (6, 1, 4, 2));
        assert_eq!(h.mapping, vec![0, 4, 1, 2, 3, 5]);
        assert_eq!(OpusHead::parse(&h.to_bytes()).unwrap(), h);
        assert!(OpusHead::parse(&surround[..15]).is_err());
        for ch in 1..=8 {
            let h = OpusHead::new(ch, 312, 44100).unwrap();
            assert_eq!(OpusHead::parse(&h.to_bytes()).unwrap(), h);
            // Every decoded channel is used exactly once.
            let mut seen = h.mapping.clone();
            seen.sort_unstable();
            assert_eq!(seen, (0..h.streams + h.coupled).collect::<Vec<_>>());
        }
        let mut bad = surround.clone();
        bad[13] = 9;
        assert!(
            OpusHead::parse(&bad).is_err(),
            "index past the decoded channels"
        );
    }
}

#[cfg(test)]
mod round_trip {
    use super::*;

    fn power(pcm: &[f32], ch: usize, c: usize, freq: f32) -> f32 {
        let w = 2.0 * std::f32::consts::PI * freq / 48_000.0;
        let (mut s1, mut s2) = (0.0f32, 0.0f32);
        for x in pcm.iter().skip(c).step_by(ch) {
            let s = x + 2.0 * w.cos() * s1 - s2;
            s2 = s1;
            s1 = s;
        }
        s1 * s1 + s2 * s2 - 2.0 * w.cos() * s1 * s2
    }

    /// Every layout from mono to 7.1: each channel's own tone comes back in
    /// its own channel, well above the others'.
    #[test]
    fn surround_channels_stay_apart() {
        for ch in 1..=8usize {
            let freqs: Vec<f32> = (0..ch).map(|k| 300.0 + 170.0 * k as f32).collect();
            let cfg = EncoderConfig {
                channels: ch,
                bitrate: 64_000 * ch as u32,
                ..EncoderConfig::default()
            };
            let mut enc = MultistreamEncoder::new(cfg).unwrap();
            let head = OpusHead::parse(&enc.head().to_bytes()).unwrap();
            assert_eq!(head.family, if ch <= 2 { 0 } else { 1 });
            let mut dec = MultistreamDecoder::from_head(&head, 48000).unwrap();
            let n = enc.frame_samples();
            let mut out = Vec::new();
            for k in 0..60 {
                let pcm: Vec<f32> = (0..n * ch)
                    .map(|i| {
                        let (t, c) = ((k * n + i / ch) as f32 / 48000.0, i % ch);
                        0.25 * (2.0 * std::f32::consts::PI * freqs[c] * t).sin()
                    })
                    .collect();
                let p = enc.encode(&pcm).unwrap();
                out.extend(dec.decode(Some(&p)).unwrap());
                assert_eq!(dec.final_range(), enc.final_range());
            }
            let steady = &out[ch * 4800..];
            for c in 0..ch {
                let own = power(steady, ch, c, freqs[c]);
                for (o, &g) in freqs.iter().enumerate() {
                    // The LFE stream is band-limited: skip tones above it.
                    if o != c && g < 3000.0 {
                        assert!(
                            own > 30.0 * power(steady, ch, c, g),
                            "{ch} channels: channel {c} carries channel {o}'s tone"
                        );
                    }
                }
            }
        }
    }
}
