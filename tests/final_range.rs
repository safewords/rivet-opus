//! Encoder/decoder agreement: for every packet, the range coder state the
//! decoder ends a frame with (RFC 6716 §4.1, "final range") must be the
//! state the encoder ended it with. A mismatch means the two disagree on
//! how many bits some symbol took — the decoder has lost sync with what was
//! coded, whether or not the output sounds plausible.
//!
//! Every mode (SILK, hybrid, CELT), mono and stereo, CBR and VBR, a spread
//! of rates and frame sizes, on speech-like and music-like signals and on a
//! hard-panned pair (one channel silent: the case that drives SILK's
//! mid-only flag and the CELT layer's intensity / dual stereo choices).

mod common;

use common::{Lcg, music, speech};
use opus::{Application, Decoder, Encoder, EncoderConfig, Mode};
use std::f32::consts::PI;

/// Left: a tone sweep with noise; right: silence, then the left's signal
/// panned hard right, then both (half a second each, repeated).
fn panned(channels: usize, seconds: f32) -> Vec<f32> {
    let n = (48_000.0 * seconds) as usize;
    let mut g = Lcg(5);
    let mut out = Vec::with_capacity(n * channels);
    for i in 0..n {
        let t = i as f32 / 48_000.0;
        let s = 0.3 * (2.0 * PI * (300.0 + 2000.0 * (t * 0.25).fract()) * t).sin()
            + 0.02 * g.next_f32();
        let phase = ((t * 2.0) as usize) % 3;
        let (l, r) = match phase {
            0 => (s, 0.0),
            1 => (0.0, s),
            _ => (s, 0.7 * s),
        };
        out.push(l);
        if channels == 2 {
            out.push(r);
        }
    }
    out
}

/// The packets whose final range differs, as (index, encoder, decoder).
fn mismatches(cfg: EncoderConfig, pcm: &[f32]) -> Vec<(usize, u32, u32)> {
    let mut enc = Encoder::new(cfg).unwrap();
    let mut dec = Decoder::new(48_000, cfg.channels).unwrap();
    let fs = enc.frame_samples() * cfg.channels;
    let mut bad = Vec::new();
    for (k, chunk) in pcm.chunks_exact(fs).enumerate() {
        let p = enc.encode(chunk).unwrap();
        dec.decode(Some(&p)).unwrap();
        if dec.final_range() != enc.final_range() {
            bad.push((k, enc.final_range(), dec.final_range()));
        }
    }
    bad
}

fn check(mode: Mode, rates: &[u32], frames: &[usize]) {
    let mut failures = Vec::new();
    for channels in [1, 2] {
        let signals: [(&str, Vec<f32>); 3] = [
            ("speech", speech(48_000, channels, 1.5)),
            ("music", music(48_000, channels, 1.5)),
            ("panned", panned(channels, 3.0)),
        ];
        for &frame in frames {
            for &rate in rates {
                // VBR, CBR, and (SILK and hybrid) VBR with LBRR frames.
                let fec_too = mode != Mode::Celt;
                for (vbr, fec) in [(true, false), (false, false), (true, true)] {
                    if fec && !fec_too {
                        continue;
                    }
                    for (name, sig) in &signals {
                        let cfg = EncoderConfig {
                            channels,
                            bitrate: rate,
                            frame_size: frame,
                            vbr,
                            fec,
                            packet_loss_percent: if fec { 10 } else { 0 },
                            mode: Some(mode),
                            application: if mode == Mode::Celt {
                                Application::Audio
                            } else {
                                Application::Voip
                            },
                            ..EncoderConfig::default()
                        };
                        let bad = mismatches(cfg, sig);
                        if !bad.is_empty() {
                            failures.push(format!(
                                "{mode:?} {channels}ch {rate} b/s {frame} {}{} {name}: {} packets differ, first {:?}",
                                if vbr { "VBR" } else { "CBR" },
                                if fec { "+FEC" } else { "" },
                                bad.len(),
                                bad[0]
                            ));
                        }
                    }
                }
            }
        }
    }
    assert!(
        failures.is_empty(),
        "final range differs:\n{}",
        failures.join("\n")
    );
}

#[test]
fn silk_final_range_agrees() {
    check(
        Mode::Silk,
        &[6_000, 10_000, 16_000, 24_000, 40_000],
        &[480, 960, 1920, 2880],
    );
}

#[test]
fn hybrid_final_range_agrees() {
    check(
        Mode::Hybrid,
        &[
            16_000, 20_000, 24_000, 28_000, 32_000, 40_000, 48_000, 64_000,
        ],
        &[480, 960, 1920],
    );
}

#[test]
fn celt_final_range_agrees() {
    check(
        Mode::Celt,
        &[6_000, 16_000, 32_000, 64_000, 128_000, 256_000, 510_000],
        &[120, 240, 480, 960, 1920],
    );
}
