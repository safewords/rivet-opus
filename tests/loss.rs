//! Packet loss: concealment (RFC 6716 §4.4) and LBRR forward error
//! correction (§4.2.5).

mod common;

use common::*;
use opus::{Application, Decoder, Encoder, EncoderConfig, Mode};

fn segment_snr(x: &[f32], y: &[f32]) -> f64 {
    let (mut s, mut e) = (0.0f64, 0.0f64);
    for (a, b) in x.iter().zip(y) {
        s += f64::from(*a) * f64::from(*a);
        e += f64::from(a - b) * f64::from(a - b);
    }
    10.0 * (s / e.max(1e-30)).log10()
}

/// Drops every seventh packet; the concealed stream keeps its length, stays
/// finite and bounded, and fades rather than clicks.
#[test]
fn concealment_in_every_mode() {
    for (mode, br, sig) in [
        (Mode::Celt, 64_000u32, music(48000, 2, 2.0)),
        (Mode::Silk, 24_000, speech(48000, 2, 2.0)),
        (Mode::Hybrid, 40_000, speech(48000, 2, 2.0)),
    ] {
        let cfg = EncoderConfig { channels: 2, bitrate: br, mode: Some(mode), vbr: false, application: Application::Voip, ..EncoderConfig::default() };
        let mut enc = Encoder::new(cfg).unwrap();
        let mut dec = Decoder::new(48000, 2).unwrap();
        let n = enc.frame_samples() * 2;
        let mut out = Vec::new();
        for (k, chunk) in sig.chunks_exact(n).enumerate() {
            let p = enc.encode(chunk).unwrap();
            let pcm = if k % 7 == 3 { dec.decode(None).unwrap() } else { dec.decode(Some(&p)).unwrap() };
            assert_eq!(pcm.len(), n, "{mode:?}: concealment keeps the packet duration");
            out.extend(pcm);
        }
        let peak_in = sig.iter().fold(0.0f32, |a, v| a.max(v.abs()));
        let peak_out = out.iter().fold(0.0f32, |a, v| a.max(v.abs()));
        assert!(out.iter().all(|v| v.is_finite()));
        assert!(peak_out < 2.0 * peak_in + 0.05, "{mode:?}: concealment overshoots ({peak_out} vs {peak_in})");
    }
}

/// With LBRR on, a lost SILK packet is recovered from the next one, and
/// the recovery is much closer to the input than concealment.
#[test]
fn lbrr_recovers_a_lost_packet() {
    for channels in [1usize, 2] {
        let sig = speech(48000, channels, 2.0);
        let cfg = EncoderConfig {
            channels,
            bitrate: 32_000 * channels as u32,
            mode: Some(Mode::Silk),
            vbr: false,
            fec: true,
            packet_loss_percent: 10,
            application: Application::Voip,
            ..EncoderConfig::default()
        };
        let mut enc = Encoder::new(cfg).unwrap();
        let n = enc.frame_samples() * channels;
        let packets: Vec<Vec<u8>> = sig.chunks_exact(n).map(|c| enc.encode(c).unwrap()).collect();
        let delay = enc.lookahead() * channels;
        let lost = 30;
        let reference = &sig[lost * n - delay..(lost + 1) * n - delay];
        let mut with_fec = Decoder::new(48000, channels).unwrap();
        let mut with_plc = Decoder::new(48000, channels).unwrap();
        let (mut fec_out, mut plc_out) = (Vec::new(), Vec::new());
        for (k, p) in packets.iter().enumerate().take(lost + 1) {
            if k == lost {
                fec_out = with_fec.decode_fec(&packets[lost + 1]).unwrap();
                plc_out = with_plc.decode(None).unwrap();
            } else {
                with_fec.decode(Some(p)).unwrap();
                with_plc.decode(Some(p)).unwrap();
            }
        }
        let fec_snr = segment_snr(reference, &fec_out);
        let plc_snr = segment_snr(reference, &plc_out);
        eprintln!("{channels}ch: lost packet recovered by LBRR at {fec_snr:.2} dB, concealed at {plc_snr:.2} dB");
        assert!(fec_snr > 5.0 && fec_snr > plc_snr + 3.0, "{channels}ch: FEC {fec_snr:.2} dB, PLC {plc_snr:.2} dB");
        // The stream carries on normally after the recovery.
        let after = with_fec.decode(Some(&packets[lost + 1])).unwrap();
        assert_eq!(after.len(), n);
    }
}
