//! Throughput benchmark: decodes the RFC 8251 test vectors and encodes the
//! decoded audio back in CELT, SILK and hybrid modes, printing how many
//! times faster than real time each runs (best of several passes).
//!
//! `cargo run --release --example opus_bench -- <vector dir> [passes] [mode filter]`

use std::path::Path;
use std::time::Instant;

fn read_bit(path: &Path) -> Vec<Vec<u8>> {
    let b = std::fs::read(path).unwrap();
    let (mut p, mut out) = (0, Vec::new());
    while p + 8 <= b.len() {
        let len = u32::from_be_bytes(b[p..p + 4].try_into().unwrap()) as usize;
        p += 8;
        out.push(b[p..p + len].to_vec());
        p += len;
    }
    out
}

fn best<F: FnMut()>(passes: usize, mut f: F) -> f64 {
    (0..passes)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64()
        })
        .fold(f64::INFINITY, f64::min)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = Path::new(args.get(1).map(String::as_str).unwrap_or("tests/vectors"));
    let passes: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(5);
    // Optional: only the encode modes whose name contains this.
    let only = args.get(3).cloned().unwrap_or_default();

    // Decode every vector at 48 kHz stereo.
    let vectors: Vec<Vec<Vec<u8>>> = (1..=12).map(|i| read_bit(&dir.join(format!("testvector{i:02}.bit")))).collect();
    let mut audio = Vec::new();
    let mut samples = 0usize;
    let t = best(passes, || {
        audio.clear();
        samples = 0;
        for v in &vectors {
            let mut dec = opus::Decoder::new(48_000, 2).unwrap();
            for p in v {
                let pcm = if p.is_empty() { dec.decode(None).unwrap() } else { dec.decode(Some(p)).unwrap() };
                samples += pcm.len() / 2;
                audio.extend_from_slice(&pcm);
            }
        }
    });
    let secs = samples as f64 / 48_000.0;
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for v in &audio {
        hash = (hash ^ u64::from(v.to_bits())).wrapping_mul(0x100_0000_01b3);
    }
    println!("decode  all vectors      {:8.1} x realtime ({secs:.1} s audio in {:.3} s, output hash {hash:016x})", secs / t, t);

    // Encode the first 60 s of the decoded audio.
    let n = audio.len().min(48_000 * 2 * 60);
    let audio = &audio[..n];
    let secs = n as f64 / 2.0 / 48_000.0;
    let modes: [(&str, Option<opus::Mode>, u32, opus::Application); 4] = [
        ("celt 128k", Some(opus::Mode::Celt), 128_000, opus::Application::Audio),
        ("celt 64k", Some(opus::Mode::Celt), 64_000, opus::Application::Audio),
        ("hybrid 32k", Some(opus::Mode::Hybrid), 32_000, opus::Application::Voip),
        ("silk 16k", Some(opus::Mode::Silk), 16_000, opus::Application::Voip),
    ];
    for (name, mode, bitrate, application) in modes {
        if !name.contains(only.as_str()) {
            continue;
        }
        let cfg = opus::EncoderConfig { channels: 2, bitrate, mode, application, ..opus::EncoderConfig::default() };
        let (mut bytes, mut hash) = (0usize, 0u64);
        let t = best(passes, || {
            let mut enc = opus::Encoder::new(cfg).unwrap();
            let f = enc.frame_samples() * 2;
            bytes = 0;
            hash = 0xcbf2_9ce4_8422_2325;
            for c in audio.chunks_exact(f) {
                let p = enc.encode(c).unwrap();
                bytes += p.len();
                for &b in &p {
                    hash = (hash ^ u64::from(b)).wrapping_mul(0x100_0000_01b3);
                }
            }
        });
        println!("encode  {name:<16} {:8.1} x realtime ({:.0} kb/s, stream hash {hash:016x})", secs / t, bytes as f64 * 8.0 / secs / 1000.0);
    }
}
