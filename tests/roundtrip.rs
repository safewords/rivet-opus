//! Encoder round trips: encode with this crate, check every packet against
//! RFC 6716 §3, decode with this crate, and measure SNR and bit rate.

mod common;

use common::*;
use opus::{Application, EncoderConfig, Mode};

struct Case {
    mode: Mode,
    bitrate: u32,
    frame: usize,
    channels: usize,
    vbr: bool,
    floor: f64,
}

fn run(cases: &[Case], signal: fn(u32, usize, f32) -> Vec<f32>, name: &str) {
    let mut failures = Vec::new();
    for c in cases {
        let sig = signal(48000, c.channels, 2.0);
        let cfg = EncoderConfig {
            channels: c.channels,
            bitrate: c.bitrate,
            frame_size: c.frame,
            vbr: c.vbr,
            mode: Some(c.mode),
            application: if c.mode == Mode::Celt {
                Application::Audio
            } else {
                Application::Voip
            },
            ..EncoderConfig::default()
        };
        let r = round_trip(cfg, &sig);
        for p in &r.packets {
            let toc = opus::packet::Toc::from_byte(p[0]);
            assert_eq!(toc.mode(), c.mode, "TOC mode");
            assert_eq!(toc.stereo, c.channels == 2, "TOC stereo bit");
        }
        let err = (r.bitrate - f64::from(c.bitrate)) / f64::from(c.bitrate) * 100.0;
        eprintln!(
            "{:?} {name} {}ch {:>6} b/s {:>4.1} ms {}: SNR {:6.2} dB, rate {:>8.0} b/s ({:+.2}%)",
            c.mode,
            c.channels,
            c.bitrate,
            c.frame as f64 / 48.0,
            if c.vbr { "VBR" } else { "CBR" },
            r.snr,
            r.bitrate,
            err
        );
        // CBR is exact; VBR may spend less on an easy signal but not more.
        let rate_ok = if c.vbr {
            err < 12.0 && err > -60.0
        } else {
            err.abs() < 1.0 || (c.frame == 480 && c.mode == Mode::Silk && err < 20.0)
        };
        if r.snr <= c.floor || !rate_ok {
            failures.push(format!(
                "{:?} {} b/s {} samples {}ch: SNR {:.2} (floor {}), rate {err:+.2}%",
                c.mode, c.bitrate, c.frame, c.channels, r.snr, c.floor
            ));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

#[test]
fn celt_round_trips() {
    let mut cases = Vec::new();
    for &frame in &[120usize, 240, 480, 960, 1920, 2880] {
        for &channels in &[1usize, 2] {
            for &(br, floor) in &[
                (16_000u32, 8.0),
                (32_000, 11.0),
                (64_000, 14.0),
                (128_000, 18.0),
                (256_000, 22.0),
                (510_000, 25.0),
            ] {
                let br = if channels == 2 { br.max(24_000) } else { br };
                // Short frames cost more overhead per second.
                // 2.5 and 5 ms frames at the lowest rates have only a few
                // bytes each: a floor only against breakage there.
                let floor = match frame {
                    120 => floor - 16.0,
                    240 => floor - 8.0,
                    _ => floor,
                };
                cases.push(Case {
                    mode: Mode::Celt,
                    bitrate: br,
                    frame,
                    channels,
                    vbr: false,
                    floor,
                });
            }
        }
    }
    for &channels in &[1usize, 2] {
        cases.push(Case {
            mode: Mode::Celt,
            bitrate: 96_000,
            frame: 960,
            channels,
            vbr: true,
            floor: 14.0,
        });
    }
    run(&cases, music, "music");
}

#[test]
fn celt_lowest_rates() {
    let cases = [
        Case {
            mode: Mode::Celt,
            bitrate: 8_000,
            frame: 960,
            channels: 1,
            vbr: false,
            floor: 2.0,
        },
        Case {
            mode: Mode::Celt,
            bitrate: 12_000,
            frame: 960,
            channels: 2,
            vbr: false,
            floor: 2.0,
        },
    ];
    run(&cases, music, "music");
}

#[test]
fn celt_transients() {
    let cases = [
        Case {
            mode: Mode::Celt,
            bitrate: 64_000,
            frame: 960,
            channels: 1,
            vbr: false,
            floor: 2.0,
        },
        Case {
            mode: Mode::Celt,
            bitrate: 128_000,
            frame: 480,
            channels: 2,
            vbr: false,
            floor: 2.0,
        },
    ];
    run(&cases, transients, "transients");
}

#[test]
fn silk_round_trips() {
    let mut cases = Vec::new();
    for &frame in &[480usize, 960, 1920, 2880] {
        for &channels in &[1usize, 2] {
            for &(br, floor) in &[
                (8_000u32, -5.0),
                (12_000, 1.0),
                (16_000, 6.0),
                (24_000, 12.0),
                (32_000, 15.0),
                (40_000, 15.0),
            ] {
                cases.push(Case {
                    mode: Mode::Silk,
                    bitrate: br * channels as u32,
                    frame,
                    channels,
                    vbr: false,
                    floor,
                });
            }
        }
    }
    for &channels in &[1usize, 2] {
        cases.push(Case {
            mode: Mode::Silk,
            bitrate: 20_000 * channels as u32,
            frame: 960,
            channels,
            vbr: true,
            floor: 10.0,
        });
    }
    run(&cases, speech, "speech");
}

#[test]
fn hybrid_round_trips() {
    let mut cases = Vec::new();
    for &frame in &[480usize, 960, 1920, 2880] {
        for &channels in &[1usize, 2] {
            for &(br, floor) in &[
                (24_000u32, 4.0),
                (32_000, 9.0),
                (48_000, 11.0),
                (64_000, 11.0),
            ] {
                cases.push(Case {
                    mode: Mode::Hybrid,
                    bitrate: br * channels as u32,
                    frame,
                    channels,
                    vbr: false,
                    floor,
                });
            }
        }
    }
    for &channels in &[1usize, 2] {
        cases.push(Case {
            mode: Mode::Hybrid,
            bitrate: 32_000 * channels as u32,
            frame: 960,
            channels,
            vbr: true,
            floor: 8.0,
        });
    }
    run(&cases, speech, "speech");
}

/// Input at every Opus rate, decoded at the same rate.
#[test]
fn other_input_rates() {
    for rate in [8000u32, 12000, 16000, 24000] {
        for channels in [1usize, 2] {
            let sig = music(rate, channels, 2.0);
            let cfg = EncoderConfig {
                sample_rate: rate,
                channels,
                bitrate: 64_000,
                frame_size: 960,
                vbr: false,
                ..EncoderConfig::default()
            };
            let r = round_trip(cfg, &sig);
            eprintln!(
                "input {rate} Hz {channels}ch: SNR {:.2} dB, rate {:.0} b/s",
                r.snr, r.bitrate
            );
            assert!(r.snr > 10.0, "{rate} Hz: {:.2}", r.snr);
        }
    }
}

/// Automatic mode choice per application and rate.
#[test]
fn automatic_modes() {
    for (app, br, want) in [
        (Application::Audio, 64_000u32, Mode::Celt),
        (Application::Voip, 16_000, Mode::Silk),
        (Application::Voip, 32_000, Mode::Hybrid),
        (Application::LowDelay, 24_000, Mode::Celt),
    ] {
        let sig = speech(48000, 1, 0.5);
        let cfg = EncoderConfig {
            channels: 1,
            bitrate: br,
            application: app,
            vbr: false,
            ..EncoderConfig::default()
        };
        let r = round_trip(cfg, &sig);
        let toc = opus::packet::Toc::from_byte(r.packets[0][0]);
        assert_eq!(toc.mode(), want, "{app:?} at {br}");
    }
}
