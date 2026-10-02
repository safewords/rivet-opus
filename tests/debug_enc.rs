mod common;
use common::*;
use opus::{Application, EncoderConfig, Mode};
#[test]
#[ignore]
fn debug_silk_sizes() {
    let sig = speech(48000, 1, 1.0);
    let cfg = EncoderConfig { channels: 1, bitrate: 8000, frame_size: 960, vbr: false, mode: Some(Mode::Silk), application: Application::Voip, ..EncoderConfig::default() };
    let mut enc = opus::Encoder::new(cfg).unwrap();
    for (k, chunk) in sig.chunks(960).enumerate() {
        let p = enc.encode(chunk).unwrap();
        let e: f32 = chunk.iter().map(|v| v * v).sum::<f32>() / 960.0;
        eprintln!("packet {k}: {} bytes, rms {:.4}", p.len(), e.sqrt());
    }
}

#[test]
#[ignore]
fn debug_silk_delay() {
    for (mode, br) in [(Mode::Silk, 32000u32), (Mode::Hybrid, 32000), (Mode::Celt, 64000)] {
        let sig = speech(48000, 1, 2.0);
        let cfg = EncoderConfig { channels: 1, bitrate: br, frame_size: 960, vbr: false, mode: Some(mode), application: Application::Voip, ..EncoderConfig::default() };
        let r = round_trip(cfg, &sig);
        let d0 = r.lookahead as i64;
        let mut best = (f64::MIN, 0i64, 0.0);
        for lag in d0 - 60..d0 + 60 {
            let (mut xy, mut yy, mut xx) = (0.0, 0.0, 0.0);
            for i in 9600..80000usize {
                let j = (i as i64 + lag) as usize;
                if j >= r.decoded.len() { break; }
                let x = f64::from(sig[i]);
                let y = f64::from(r.decoded[j]);
                xy += x * y; yy += y * y; xx += x * x;
            }
            let c = xy / (xx * yy).sqrt();
            if c > best.0 { best = (c, lag, xy / yy); }
        }
        eprintln!("{mode:?}: lookahead {d0}, best lag {} corr {:.4} gain {:.3}", best.1, best.0, best.2);
    }
}

#[test]
#[ignore]
fn debug_silk_frames() {
    let sig = speech(48000, 1, 1.0);
    let cfg = EncoderConfig { channels: 1, bitrate: 32000, frame_size: 960, vbr: false, mode: Some(Mode::Silk), application: Application::Voip, ..EncoderConfig::default() };
    let r = round_trip(cfg, &sig);
    let d = r.lookahead;
    for f in 0..48 {
        let (mut ex, mut ey, mut ee) = (0.0f64, 0.0f64, 0.0f64);
        for i in f * 960..(f + 1) * 960 {
            if i + d >= r.decoded.len() { break; }
            let x = f64::from(sig[i]);
            let y = f64::from(r.decoded[i + d]);
            ex += x * x; ey += y * y; ee += (x - y).powi(2);
        }
        eprintln!("frame {f}: in {:.1} dB out {:.1} dB snr {:.1}", 10.0 * (ex / 960.0 + 1e-12).log10(), 10.0 * (ey / 960.0 + 1e-12).log10(), 10.0 * (ex / ee).log10());
    }
}
