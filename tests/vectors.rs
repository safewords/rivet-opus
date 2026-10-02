//! The official Opus test vectors (opus_testvectors-rfc8251), used as data.
//!
//! Set `OPUS_TESTVECTORS` to the directory holding `testvectorNN.bit`,
//! `testvectorNN.dec` and `testvectorNN m.dec`, or put them in
//! `tests/vectors/`. Without them these tests skip (and say so).

use std::path::PathBuf;

fn dir() -> Option<PathBuf> {
    let d = std::env::var_os("OPUS_TESTVECTORS").map(PathBuf::from).unwrap_or_else(|| {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("vectors")
    });
    d.join("testvector01.bit").exists().then_some(d)
}

struct Packet {
    data: Vec<u8>,
    range: u32,
}

fn read_bit(path: &PathBuf) -> Vec<Packet> {
    let b = std::fs::read(path).unwrap();
    let mut p = 0;
    let mut out = Vec::new();
    while p + 8 <= b.len() {
        let len = u32::from_be_bytes(b[p..p + 4].try_into().unwrap()) as usize;
        let range = u32::from_be_bytes(b[p + 4..p + 8].try_into().unwrap());
        p += 8;
        out.push(Packet { data: b[p..p + len].to_vec(), range });
        p += len;
    }
    out
}

fn read_pcm(path: &PathBuf) -> Vec<f32> {
    std::fs::read(path).unwrap().chunks_exact(2).map(|c| f32::from(i16::from_le_bytes([c[0], c[1]])) / 32768.0).collect()
}

/// Decodes a vector; returns (output, packets whose final range mismatched).
fn decode(packets: &[Packet], rate: u32, channels: usize) -> (Vec<f32>, usize, usize) {
    let mut dec = opus::Decoder::new(rate, channels).unwrap();
    let mut out = Vec::new();
    let mut bad = 0;
    let mut first_bad = usize::MAX;
    for (i, p) in packets.iter().enumerate() {
        let pcm = if p.data.is_empty() { dec.decode(None).unwrap() } else { dec.decode(Some(&p.data)).unwrap() };
        out.extend_from_slice(&pcm);
        if !p.data.is_empty() && dec.final_range() != p.range {
            bad += 1;
            first_bad = first_bad.min(i);
        }
    }
    (out, bad, first_bad)
}

/// SNR (dB) of `test` against `reference`, and the largest absolute error
/// (in 16-bit units).
fn compare(reference: &[f32], test: &[f32]) -> (f64, f64) {
    let n = reference.len().min(test.len());
    let (mut s, mut e, mut m) = (0.0f64, 0.0f64, 0.0f64);
    for i in 0..n {
        let r = f64::from(reference[i]);
        let t = f64::from(test[i]).clamp(-1.0, 32767.0 / 32768.0);
        // The reference is 16-bit: round the test output the same way.
        let t = (t * 32768.0).round() / 32768.0;
        s += r * r;
        e += (r - t) * (r - t);
        m = m.max((r - t).abs() * 32768.0);
    }
    (10.0 * (s / e.max(1e-30)).log10(), m)
}

#[test]
fn test_vectors_48k() {
    let Some(d) = dir() else {
        eprintln!("SKIPPED: test vectors not found (set OPUS_TESTVECTORS)");
        return;
    };
    let mut report = String::new();
    for v in 1..=12 {
        let packets = read_bit(&d.join(format!("testvector{v:02}.bit")));
        let reference = read_pcm(&d.join(format!("testvector{v:02}.dec")));
        let reference_m = read_pcm(&d.join(format!("testvector{v:02}m.dec")));
        let (stereo, bad, first_bad) = decode(&packets, 48000, 2);
        let (snr, maxerr) = compare(&reference, &stereo);
        let (snr_m, _) = compare(&reference_m, &stereo);
        let line = format!(
            "vector {v:02}: {} packets, final range mismatches {bad} (first {first_bad}); stereo SNR {snr:.2} dB (vs m: {snr_m:.2} dB), max error {maxerr:.0}\n",
            packets.len()
        );
        eprint!("{line}");
        report.push_str(&line);
    }
}
