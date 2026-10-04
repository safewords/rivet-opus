//! Stereo separation: hard-panned tones through the encoder and decoder,
//! measuring how much of each side's tone the other side carries.
//!
//! `cargo test --release --test stereo_separation -- --ignored --nocapture`
//! prints the full tables; the regular tests hold the thresholds.
//!
//! The tones at multiples of 200 Hz sit on CELT band edges (RFC 6716
//! Table 55; the first eight bands are 200 Hz wide), so each is split
//! between two bands and shares one with its neighbour on the other side.
//! A signal at a multiple of 50 Hz repeats every 20 ms frame, so the coding
//! error repeats too and lands exactly on such frequencies: what is
//! measured at the other side's tone is part crosstalk, part the channel's
//! own coding noise. The "alone" figures (the other side silent) are the
//! latter.

mod common;

use common::Lcg;
use opus::{
    Application, Decoder, Encoder, EncoderConfig, Mode, MultistreamDecoder, MultistreamEncoder,
    OpusHead,
};
use std::f32::consts::PI;

const LEVEL: f32 = 0.25;
const SECONDS: f32 = 2.0;

/// Amplitude of `freq` (a multiple of 10 Hz) in channel `c` of interleaved
/// `pcm`, over a whole number of 100 ms from the middle half (Goertzel).
fn amplitude(pcm: &[f32], ch: usize, c: usize, freq: f32) -> f32 {
    let frames = pcm.len() / ch;
    let start = frames / 4;
    let len = (frames / 2) / 4800 * 4800;
    let k = 2.0 * (2.0 * std::f64::consts::PI * f64::from(freq) / 48_000.0).cos();
    let (mut s1, mut s2) = (0.0f64, 0.0f64);
    for i in start..start + len {
        let s = f64::from(pcm[i * ch + c]) + k * s1 - s2;
        s2 = s1;
        s1 = s;
    }
    let power = s1 * s1 + s2 * s2 - k * s1 * s2;
    (2.0 * power.max(0.0).sqrt() / len as f64) as f32
}

/// Relative to the tone level.
fn db(x: f32) -> f32 {
    20.0 * (x.max(1e-9) / LEVEL).log10()
}

/// Tone `fl` left, `fr` right (0 for silence), plus (if `noise`)
/// independent white noise in each channel some 26 dB down.
fn signal(fl: f32, fr: f32, noise: bool) -> Vec<f32> {
    let n = (48_000.0 * SECONDS) as usize;
    let (mut gl, mut gr) = (Lcg(1), Lcg(2));
    let tone = |f: f32, t: f32| {
        if f > 0.0 {
            LEVEL * (2.0 * PI * f * t).sin()
        } else {
            0.0
        }
    };
    let mut out = Vec::with_capacity(2 * n);
    for i in 0..n {
        let t = i as f32 / 48_000.0;
        let (nl, nr) = if noise {
            (0.025 * gl.next_f32(), 0.025 * gr.next_f32())
        } else {
            (0.0, 0.0)
        };
        out.push(tone(fl, t) + nl);
        out.push(tone(fr, t) + nr);
    }
    out
}

/// Encodes and decodes; `None` if the encoder and decoder disagree on a
/// packet's final range.
fn round_trip(cfg: EncoderConfig, pcm: &[f32]) -> Option<Vec<f32>> {
    let mut enc = Encoder::new(cfg).unwrap();
    let mut dec = Decoder::new(48_000, 2).unwrap();
    let fs = enc.frame_samples() * 2;
    let mut out = Vec::new();
    for chunk in pcm.chunks_exact(fs) {
        let p = enc.encode(chunk).unwrap();
        out.extend(dec.decode(Some(&p)).unwrap());
        if dec.final_range() != enc.final_range() {
            return None;
        }
    }
    Some(out)
}

struct Leak {
    /// The worse of left's tone in the right channel and right's in the left.
    both: f32,
    /// The same with the other side silent: each channel's own coding noise
    /// at the other's frequency.
    alone: f32,
    /// What the input itself carries there (the broadband case).
    input: f32,
}

fn leak(cfg: EncoderConfig, fl: f32, fr: f32, noise: bool) -> Option<Leak> {
    let x = signal(fl, fr, noise);
    let y = round_trip(cfg, &x)?;
    let yl = round_trip(cfg, &signal(fl, 0.0, noise))?;
    let yr = round_trip(cfg, &signal(0.0, fr, noise))?;
    Some(Leak {
        both: db(amplitude(&y, 2, 1, fl).max(amplitude(&y, 2, 0, fr))),
        alone: db(amplitude(&yr, 2, 1, fl).max(amplitude(&yl, 2, 0, fr))),
        input: db(amplitude(&x, 2, 1, fl).max(amplitude(&x, 2, 0, fr))),
    })
}

fn cfg(bitrate: u32, frame_size: usize, mode: Option<Mode>) -> EncoderConfig {
    EncoderConfig {
        channels: 2,
        bitrate,
        frame_size,
        mode,
        application: if mode == Some(Mode::Hybrid) {
            Application::Voip
        } else {
            Application::Audio
        },
        ..EncoderConfig::default()
    }
}

/// 5.1 at rivet's default 320 kb/s: a coupled stream's share (weights 1.5
/// per coupled stream, 1 for the centre, 0.25 for the LFE).
const COUPLED_51: u32 = 320_000 * 6 / 17;

#[test]
#[ignore = "a report; run with --ignored --nocapture"]
fn separation_table() {
    let pairs: [(f32, f32, bool, &str); 7] = [
        (440.0, 660.0, false, "440/660"),
        (400.0, 600.0, false, "400/600"),
        (410.0, 610.0, false, "410/610"),
        (1000.0, 1200.0, false, "1000/1200"),
        (3000.0, 5000.0, false, "3k/5k"),
        (2900.0, 3100.0, false, "2.9k/3.1k"),
        (440.0, 660.0, true, "440/660+noise"),
    ];
    let rates = [
        24_000, 32_000, 48_000, 64_000, 96_000, COUPLED_51, 128_000, 192_000, 256_000,
    ];
    let variants: [(&str, usize, Option<Mode>); 5] = [
        ("auto 20 ms", 960, None),
        ("CELT 10 ms", 480, Some(Mode::Celt)),
        ("CELT 5 ms", 240, Some(Mode::Celt)),
        ("CELT 40 ms", 1920, Some(Mode::Celt)),
        ("hybrid 20 ms", 960, Some(Mode::Hybrid)),
    ];
    eprintln!(
        "leak with both tones / with the other side silent, dB re the tone; noise case: / input's own level"
    );
    for (vname, frame, mode) in variants {
        eprint!("\n{vname:>12} kb/s");
        for p in pairs {
            eprint!(" | {:>13}", p.3);
        }
        eprintln!();
        for &r in &rates {
            if mode == Some(Mode::Hybrid) && r > 96_000 {
                continue;
            }
            eprint!("{:>17.1}", r as f32 / 1000.0);
            for (fl, fr, noise, _) in pairs {
                match leak(cfg(r, frame, mode), fl, fr, noise) {
                    Some(l) if noise => eprint!(" | {:>6.1} /{:>5.1}", l.both, l.input),
                    Some(l) => eprint!(" | {:>6.1} /{:>5.1}", l.both, l.alone),
                    None => eprint!(" | {:>13}", "range differs"),
                }
            }
            eprintln!();
        }
    }
}

/// Hybrid stereo (SILK below 8 kHz, CELT from band 17) at 16 to 64 kb/s.
#[test]
#[ignore = "a report; run with --ignored --nocapture"]
fn hybrid_separation_table() {
    let pairs: [(f32, f32, &str); 4] = [
        (440.0, 660.0, "440/660"),
        (400.0, 600.0, "400/600"),
        (1000.0, 1200.0, "1000/1200"),
        (3000.0, 5000.0, "3k/5k"),
    ];
    eprintln!("hybrid 20 ms: leak with both tones / with the other side silent, dB re the tone");
    eprint!("{:>10}", "kb/s");
    for p in pairs {
        eprint!(" | {:>13}", p.2);
    }
    eprintln!();
    for r in [
        16_000, 20_000, 24_000, 28_000, 32_000, 40_000, 48_000, 64_000,
    ] {
        eprint!("{:>10.1}", r as f32 / 1000.0);
        for (fl, fr, _) in pairs {
            match leak(cfg(r, 960, Some(Mode::Hybrid)), fl, fr, false) {
                Some(l) => eprint!(" | {:>6.1} /{:>5.1}", l.both, l.alone),
                None => eprint!(" | {:>13}", "range differs"),
            }
        }
        eprintln!();
    }
}

/// Hybrid stereo keeps a hard-panned pair 20 dB apart from 32 kb/s up.
/// With 55 % of the rate for the SILK layer the pair crossed over at -14
/// dB at 32 kb/s (and at -26 to -29 dB at 40-48 kb/s); with 80 % it is
/// -24 to -33 dB at 32 kb/s. Below that SILK's mid and side get too few
/// bits each for their coding errors to cancel in L and R.
#[test]
fn hybrid_hard_panned_tones_stay_apart() {
    let mut failures = Vec::new();
    for rate in [32_000, 40_000, 48_000, 64_000] {
        let l = leak(cfg(rate, 960, Some(Mode::Hybrid)), 440.0, 660.0, false)
            .expect("final ranges agree");
        if l.both > -20.0 {
            failures.push(format!("{rate} b/s: {:.1} dB", l.both));
        }
    }
    assert!(
        failures.is_empty(),
        "hybrid crosstalk above -20 dB: {failures:#?}"
    );
}

/// rivet's 5.1 test: FL FR FC LFE BL BR carry 400, 600, 800, 50, 1000 and
/// 1200 Hz, given in Vorbis order (FL FC FR BL BR LFE), at 320 kb/s.
/// Returns the worst level of another channel's tone in each channel, in
/// Vorbis order.
fn surround_51(report: bool) -> [f32; 6] {
    let tones = [400.0f32, 800.0, 600.0, 1000.0, 1200.0, 50.0];
    let names = ["FL", "FC", "FR", "BL", "BR", "LFE"];
    let cfg = EncoderConfig {
        channels: 6,
        bitrate: 320_000,
        ..EncoderConfig::default()
    };
    let mut enc = MultistreamEncoder::new(cfg).unwrap();
    let head = OpusHead::parse(&enc.head().to_bytes()).unwrap();
    let mut dec = MultistreamDecoder::from_head(&head, 48_000).unwrap();
    let n = enc.frame_samples();
    let mut out = Vec::new();
    for k in 0..(48_000.0 * SECONDS) as usize / n {
        let pcm: Vec<f32> = (0..n * 6)
            .map(|i| LEVEL * (2.0 * PI * tones[i % 6] * ((k * n + i / 6) as f32 / 48_000.0)).sin())
            .collect();
        out.extend(dec.decode(Some(&enc.encode(&pcm).unwrap())).unwrap());
        assert_eq!(dec.final_range(), enc.final_range());
    }
    let mut worst = [f32::NEG_INFINITY; 6];
    for c in 0..6 {
        for (o, &f) in tones.iter().enumerate() {
            if o != c {
                let l = db(amplitude(&out, 6, c, f));
                if report {
                    eprintln!(
                        "5.1: {} carries {}'s {f} Hz at {l:.1} dB",
                        names[c], names[o]
                    );
                }
                worst[c] = worst[c].max(l);
            }
        }
    }
    worst
}

#[test]
#[ignore = "a report; run with --ignored --nocapture"]
fn surround_51_report() {
    let w = surround_51(true);
    eprintln!("worst per channel (FL FC FR BL BR LFE): {w:.1?}");
}

/// Hard-panned tones sharing a band: the dual-stereo decision must see the
/// bands that hold the signal, or mid/side mixes each side's coding error
/// into the other (RFC 6716 §5.3.5 is the decision; the L1 norms are taken
/// on the spectrum, not the unit-norm bands). With mid/side these leaked
/// at -27 to -36 dB; dual stereo leaves each channel's own coding noise,
/// -42 dB and below.
#[test]
fn hard_panned_tones_stay_apart() {
    let mut failures = Vec::new();
    for rate in [96_000, COUPLED_51, 128_000, 192_000, 256_000] {
        for (fl, fr) in [(400.0, 600.0), (410.0, 610.0), (440.0, 660.0)] {
            let l = leak(cfg(rate, 960, None), fl, fr, false).expect("final ranges agree");
            if l.both > -38.0 {
                failures.push(format!("{rate} b/s, {fl}/{fr} Hz: {:.1} dB", l.both));
            }
        }
    }
    assert!(failures.is_empty(), "crosstalk above -38 dB: {failures:#?}");
}

/// The 5.1 layout at 320 kb/s keeps every channel's tone out of the others
/// by 30 dB (-27 dB before the dual-stereo decision looked at the spectrum;
/// what is left, some -34 dB, is each channel's own coding noise).
#[test]
fn surround_51_channels_stay_apart() {
    let w = surround_51(false);
    assert!(
        w.iter().all(|&l| l < -30.0),
        "worst leak per channel (FL FC FR BL BR LFE): {w:.1?}"
    );
}
