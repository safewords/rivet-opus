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

/// Decodes a vector; returns (output, packets whose final range mismatched,
/// first mismatching packet).
fn decode(packets: &[Packet], rate: u32, channels: usize, no_inversion: bool) -> (Vec<f32>, usize, Option<usize>) {
    let mut dec = opus::Decoder::new(rate, channels).unwrap();
    dec.set_phase_inversion_disabled(no_inversion);
    let mut out = Vec::new();
    let mut bad = 0;
    let mut first_bad = None;
    for (i, p) in packets.iter().enumerate() {
        let pcm = if p.data.is_empty() { dec.decode(None).unwrap() } else { dec.decode(Some(&p.data)).unwrap() };
        out.extend_from_slice(&pcm);
        if !p.data.is_empty() && dec.final_range() != p.range {
            bad += 1;
            first_bad.get_or_insert(i);
        }
    }
    (out, bad, first_bad)
}

/// SNR (dB) of `test` against `reference`, and the largest absolute error
/// in 16-bit units. The test output is rounded to 16 bits like the
/// reference.
fn compare(reference: &[f32], test: &[f32]) -> (f64, f64) {
    let n = reference.len().min(test.len());
    let (mut s, mut e, mut m) = (0.0f64, 0.0f64, 0.0f64);
    for i in 0..n {
        let r = f64::from(reference[i]);
        let t = (f64::from(test[i]) * 32768.0).round().clamp(-32768.0, 32767.0) / 32768.0;
        s += r * r;
        e += (r - t) * (r - t);
        m = m.max((r - t).abs() * 32768.0);
    }
    (10.0 * (s / e.max(1e-30)).log10(), m)
}

/// The lowest SNR each vector is allowed (dB): the CELT-only vectors are
/// float rounding away from the reference; the SILK and hybrid ones carry
/// the difference of the (non-normative, float versus fixed-point) SILK
/// synthesis and resampling.
const FLOOR: [f64; 12] = [100.0, 45.0, 45.0, 42.0, 40.0, 40.0, 95.0, 80.0, 80.0, 55.0, 100.0, 43.0];

#[test]
fn test_vectors_48k() {
    let Some(d) = dir() else {
        eprintln!("SKIPPED: test vectors not found (set OPUS_TESTVECTORS)");
        return;
    };
    for v in 1..=12 {
        let packets = read_bit(&d.join(format!("testvector{v:02}.bit")));
        let reference = read_pcm(&d.join(format!("testvector{v:02}.dec")));
        let reference_m = read_pcm(&d.join(format!("testvector{v:02}m.dec")));
        let (stereo, bad, first_bad) = decode(&packets, 48000, 2, false);
        let (snr, maxerr) = compare(&reference, &stereo);
        let (stereo_m, bad_m, _) = decode(&packets, 48000, 2, true);
        let (snr_m, maxerr_m) = compare(&reference_m, &stereo_m);
        let (mono, bad_mono, _) = decode(&packets, 48000, 1, true);
        let mono_ref: Vec<f32> = reference_m.chunks_exact(2).map(|c| 0.5 * (c[0] + c[1])).collect();
        let (snr_mono, maxerr_mono) = compare(&mono_ref, &mono);
        eprintln!(
            "vector {v:02}: {} packets, final-range mismatches {bad}/{bad_m}/{bad_mono};              stereo SNR {snr:.2} dB max err {maxerr:.0}; stereo (no inversion) vs m {snr_m:.2} dB max err {maxerr_m:.0};              mono vs m downmix {snr_mono:.2} dB max err {maxerr_mono:.0}",
            packets.len()
        );
        assert_eq!(bad + bad_m + bad_mono, 0, "vector {v}: range coder state differs from the reference (first at packet {first_bad:?})");
        assert_eq!(stereo.len(), reference.len(), "vector {v}: length");
        let floor = FLOOR[v - 1];
        assert!(snr > floor && snr_m > floor, "vector {v}: SNR {snr:.2} / {snr_m:.2} below {floor} dB");
    }
}

/// Every vector also decodes at the other rates and channel counts with the
/// reference's final range on every packet and the right length (the
/// reference outputs exist only at 48 kHz).
#[test]
fn test_vectors_other_rates() {
    let Some(d) = dir() else {
        eprintln!("SKIPPED: test vectors not found (set OPUS_TESTVECTORS)");
        return;
    };
    for v in 1..=12 {
        let packets = read_bit(&d.join(format!("testvector{v:02}.bit")));
        let n48 = read_pcm(&d.join(format!("testvector{v:02}.dec"))).len() / 2;
        for rate in [8000u32, 12000, 16000, 24000] {
            for channels in [1usize, 2] {
                let (pcm, bad, first) = decode(&packets, rate, channels, false);
                assert_eq!(bad, 0, "vector {v} at {rate} Hz x{channels}: first mismatch {first:?}");
                assert_eq!(pcm.len(), n48 * rate as usize / 48000 * channels, "vector {v} at {rate} Hz");
                assert!(pcm.iter().all(|x| x.is_finite() && x.abs() <= 2.0));
            }
        }
    }
}
