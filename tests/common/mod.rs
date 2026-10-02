//! Shared helpers for the round-trip tests: test signals, encode/decode,
//! alignment and measurement.
#![allow(dead_code)]

use opus::{Decoder, Encoder, EncoderConfig};

/// A deterministic generator.
pub struct Lcg(pub u32);
impl Lcg {
    pub fn next_f32(&mut self) -> f32 {
        self.0 = self.0.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        (self.0 >> 8) as f32 / (1u32 << 23) as f32 - 1.0
    }
}

/// Music-like: a few harmonic tones with vibrato and a slow envelope,
/// different in each channel, plus a little noise.
pub fn music(rate: u32, channels: usize, seconds: f32) -> Vec<f32> {
    let n = (rate as f32 * seconds) as usize;
    let mut g = Lcg(7);
    let mut out = vec![0.0f32; n * channels];
    for ch in 0..channels {
        let base = [220.0f32, 277.2, 329.6][ch % 3];
        for i in 0..n {
            let t = i as f32 / rate as f32;
            let vib = 1.0 + 0.003 * (2.0 * std::f32::consts::PI * 5.0 * t).sin();
            let mut s = 0.0;
            for h in 1..8 {
                let f = base * h as f32 * vib;
                if f < rate as f32 / 2.0 {
                    s += (0.5 / h as f32) * (2.0 * std::f32::consts::PI * f * t + h as f32).sin();
                }
            }
            let env = 0.6 + 0.4 * (2.0 * std::f32::consts::PI * 0.7 * t).sin();
            out[i * channels + ch] = 0.3 * env * s + 0.01 * g.next_f32();
        }
    }
    out
}

/// Speech-like: a glottal pulse train with a moving pitch through two
/// resonances, in syllables separated by pauses, with unvoiced bursts.
pub fn speech(rate: u32, channels: usize, seconds: f32) -> Vec<f32> {
    let n = (rate as f32 * seconds) as usize;
    let mut g = Lcg(11);
    let mut mono = vec![0.0f32; n];
    let (mut y1, mut y2, mut z1, mut z2) = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
    let mut phase = 0.0f32;
    for (i, m) in mono.iter_mut().enumerate() {
        let t = i as f32 / rate as f32;
        let syl = (t * 4.0) % 1.0;
        let voiced = syl < 0.55;
        let unvoiced = (0.6..0.7).contains(&syl);
        let f0 = 120.0 + 30.0 * (2.0 * std::f32::consts::PI * 0.5 * t).sin();
        phase += f0 / rate as f32;
        let mut e = 0.0;
        if voiced && phase >= 1.0 {
            e = 1.0;
        }
        if phase >= 1.0 {
            phase -= 1.0;
        }
        if unvoiced {
            e += 0.15 * g.next_f32();
        }
        // Two resonators (formants near 700 Hz and 1800 Hz, below the
        // Nyquist rate of every Opus rate).
        let res = |f: f32, bwd: f32| {
            let r = (-std::f32::consts::PI * bwd / rate as f32).exp();
            let th = 2.0 * std::f32::consts::PI * f / rate as f32;
            (2.0 * r * th.cos(), -r * r)
        };
        let (a1, a2) = res(700.0, 120.0);
        let y = e + a1 * y1 + a2 * y2;
        y2 = y1;
        y1 = y;
        let (b1, b2) = res(1800.0, 200.0);
        let z = y + b1 * z1 + b2 * z2;
        z2 = z1;
        z1 = z;
        *m = 0.02 * z;
    }
    let peak = mono.iter().fold(0.0f32, |a, v| a.max(v.abs())).max(1e-9);
    let mut out = vec![0.0f32; n * channels];
    for i in 0..n {
        for ch in 0..channels {
            out[i * channels + ch] = 0.5 * mono[i] / peak * if ch == 0 { 1.0 } else { 0.8 };
        }
    }
    out
}

/// Transients: noise bursts and clicks over a quiet tone.
pub fn transients(rate: u32, channels: usize, seconds: f32) -> Vec<f32> {
    let n = (rate as f32 * seconds) as usize;
    let mut g = Lcg(3);
    let mut out = vec![0.0f32; n * channels];
    for i in 0..n {
        let t = i as f32 / rate as f32;
        let burst = (t * 3.0) % 1.0 < 0.05;
        for ch in 0..channels {
            let tone = 0.05 * (2.0 * std::f32::consts::PI * 440.0 * t).sin();
            out[i * channels + ch] = tone + if burst { 0.4 * g.next_f32() } else { 0.0 };
        }
    }
    out
}

pub struct RoundTrip {
    pub packets: Vec<Vec<u8>>,
    pub decoded: Vec<f32>,
    pub bitrate: f64,
    pub snr: f64,
    pub lookahead: usize,
}

/// Encodes `pcm` with `cfg`, checks every packet against RFC 6716 §3,
/// decodes at `cfg.sample_rate`, and measures the SNR against the input
/// (after the encoder's lookahead), skipping the first 100 ms.
pub fn round_trip(cfg: EncoderConfig, pcm: &[f32]) -> RoundTrip {
    let c = cfg.channels;
    let mut enc = Encoder::new(cfg).unwrap();
    let mut dec = Decoder::new(cfg.sample_rate, c).unwrap();
    let fs = enc.frame_samples() * c;
    let mut packets = Vec::new();
    let mut decoded = Vec::new();
    for chunk in pcm.chunks(fs) {
        if chunk.len() < fs {
            break;
        }
        let p = enc.encode(chunk).unwrap();
        let parsed = opus::packet::parse(&p).expect("packet breaks RFC 6716 section 3");
        assert_eq!(parsed.samples(), cfg.frame_size, "packet duration");
        assert!(parsed.frames.iter().all(|f| f.len() <= 1275));
        let out = dec.decode(Some(&p)).unwrap();
        assert_eq!(dec.final_range(), enc.final_range(), "encoder and decoder final ranges differ");
        decoded.extend_from_slice(&out);
        packets.push(p);
    }
    let delay = enc.lookahead() * cfg.sample_rate as usize / 48000;
    let skip = cfg.sample_rate as usize / 10;
    let (mut s, mut e) = (0.0f64, 0.0f64);
    let frames = decoded.len() / c;
    for i in skip..frames.saturating_sub(delay) {
        for ch in 0..c {
            let x = f64::from(pcm[i * c + ch]);
            let y = f64::from(decoded[(i + delay) * c + ch]);
            s += x * x;
            e += (x - y) * (x - y);
        }
    }
    let bytes: usize = packets.iter().map(|p| p.len()).sum();
    let seconds = packets.len() as f64 * cfg.frame_size as f64 / 48000.0;
    RoundTrip { bitrate: bytes as f64 * 8.0 / seconds, snr: 10.0 * (s / e.max(1e-30)).log10(), packets, decoded, lookahead: delay }
}
