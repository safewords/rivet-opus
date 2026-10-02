//! Encoder round trips: encode, check the packets, decode, measure.

mod common;

use common::*;
use opus::{Application, EncoderConfig, Mode};

#[test]
fn celt_round_trip_report() {
    for &(bitrate, frame) in &[(64_000u32, 960usize), (128_000, 960), (32_000, 960), (96_000, 480), (96_000, 240), (96_000, 120), (64_000, 1920), (64_000, 2880)] {
        for channels in [1usize, 2] {
            for (name, sig) in [("music", music(48000, channels, 3.0)), ("transients", transients(48000, channels, 3.0))] {
                let cfg = EncoderConfig {
                    channels,
                    bitrate,
                    frame_size: frame,
                    vbr: false,
                    mode: Some(Mode::Celt),
                    application: Application::Audio,
                    ..EncoderConfig::default()
                };
                let r = round_trip(cfg, &sig);
                eprintln!("CELT {name} {channels}ch {bitrate} b/s {frame}: SNR {:.2} dB, rate {:.0} b/s", r.snr, r.bitrate);
            }
        }
    }
}

#[test]
fn silk_and_hybrid_report() {
    for &(mode, bitrate, frame) in &[
        (Mode::Silk, 8_000u32, 960usize),
        (Mode::Silk, 12_000, 960),
        (Mode::Silk, 20_000, 960),
        (Mode::Silk, 32_000, 960),
        (Mode::Silk, 16_000, 480),
        (Mode::Silk, 16_000, 1920),
        (Mode::Silk, 16_000, 2880),
        (Mode::Hybrid, 24_000, 960),
        (Mode::Hybrid, 32_000, 960),
        (Mode::Hybrid, 48_000, 480),
        (Mode::Hybrid, 32_000, 1920),
    ] {
        for channels in [1usize, 2] {
            let sig = speech(48000, channels, 3.0);
            let cfg = EncoderConfig {
                channels,
                bitrate: bitrate * channels as u32,
                frame_size: frame,
                vbr: false,
                mode: Some(mode),
                application: Application::Voip,
                ..EncoderConfig::default()
            };
            let r = round_trip(cfg, &sig);
            eprintln!("{mode:?} speech {channels}ch {} b/s {frame}: SNR {:.2} dB, rate {:.0} b/s", cfg.bitrate, r.snr, r.bitrate);
        }
    }
}
