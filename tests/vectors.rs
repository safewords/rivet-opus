//! The official Opus test vectors (opus_testvectors-rfc8251), used as data.
//!
//! Set `OPUS_TESTVECTORS` to the directory holding `testvectorNN.bit`,
//! `testvectorNN.dec` and `testvectorNN m.dec`, or put them in
//! `tests/vectors/`. Without them these tests skip (and say so).

#![allow(clippy::needless_range_loop, clippy::chunks_exact_to_as_chunks)]

use std::path::PathBuf;

mod opus_compare;

fn dir() -> Option<PathBuf> {
    let d = std::env::var_os("OPUS_TESTVECTORS")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("vectors")
        });
    let found = d.join("testvector01.bit").exists();
    assert!(
        found || std::env::var_os("OPUS_REQUIRE_VECTORS").is_none(),
        "OPUS_REQUIRE_VECTORS is set but no test vectors are in {}",
        d.display()
    );
    found.then_some(d)
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
        out.push(Packet {
            data: b[p..p + len].to_vec(),
            range,
        });
        p += len;
    }
    out
}

fn read_pcm(path: &PathBuf) -> Vec<f32> {
    std::fs::read(path)
        .unwrap()
        .chunks_exact(2)
        .map(|c| f32::from(i16::from_le_bytes([c[0], c[1]])) / 32768.0)
        .collect()
}

/// Decodes a vector; returns (output, packets whose final range mismatched,
/// first mismatching packet).
fn decode(
    packets: &[Packet],
    rate: u32,
    channels: usize,
    no_inversion: bool,
) -> (Vec<f32>, usize, Option<usize>) {
    let mut dec = opus::Decoder::new(rate, channels).unwrap();
    dec.set_phase_inversion_disabled(no_inversion);
    let mut out = Vec::new();
    let mut bad = 0;
    let mut first_bad = None;
    for (i, p) in packets.iter().enumerate() {
        let pcm = if p.data.is_empty() {
            dec.decode(None).unwrap()
        } else {
            dec.decode(Some(&p.data)).unwrap()
        };
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
        let t = (f64::from(test[i]) * 32768.0)
            .round()
            .clamp(-32768.0, 32767.0)
            / 32768.0;
        s += r * r;
        e += (r - t) * (r - t);
        m = m.max((r - t).abs() * 32768.0);
    }
    (10.0 * (s / e.max(1e-30)).log10(), m)
}

/// The lowest SNR allowed for the CELT-only vectors (dB), which are float
/// rounding away from the reference. Vectors with SILK content carry the
/// difference of the non-normative SILK synthesis and resampler (which
/// RFC 6716 §4.2.9 lets a decoder choose freely): their waveform SNR is
/// reported but not a criterion; the RFC's criterion,
/// [`test_vectors_conformance`], applies to them.
const CELT_FLOOR: [Option<f64>; 12] = [
    Some(100.0),
    None,
    None,
    None,
    None,
    None,
    Some(95.0),
    None,
    None,
    None,
    Some(100.0),
    None,
];

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
        let mono_ref: Vec<f32> = reference_m
            .chunks_exact(2)
            .map(|c| 0.5 * (c[0] + c[1]))
            .collect();
        let (snr_mono, maxerr_mono) = compare(&mono_ref, &mono);
        eprintln!(
            "vector {v:02}: {} packets, final-range mismatches {bad}/{bad_m}/{bad_mono}; stereo SNR {snr:.2} dB max err {maxerr:.0}; stereo (no inversion) vs m {snr_m:.2} dB max err {maxerr_m:.0}; mono vs m downmix {snr_mono:.2} dB max err {maxerr_mono:.0}",
            packets.len()
        );
        assert_eq!(
            bad + bad_m + bad_mono,
            0,
            "vector {v}: range coder state differs from the reference (first at packet {first_bad:?})"
        );
        assert_eq!(stereo.len(), reference.len(), "vector {v}: length");
        if let Some(floor) = CELT_FLOOR[v - 1] {
            assert!(
                snr > floor && snr_m > floor,
                "vector {v}: SNR {snr:.2} / {snr_m:.2} below {floor} dB"
            );
        }
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
                assert_eq!(
                    bad, 0,
                    "vector {v} at {rate} Hz x{channels}: first mismatch {first:?}"
                );
                assert_eq!(
                    pcm.len(),
                    n48 * rate as usize / 48000 * channels,
                    "vector {v} at {rate} Hz"
                );
                assert!(pcm.iter().all(|x| x.is_finite() && x.abs() <= 2.0));
            }
        }
    }
}

/// The lowest opus_compare quality accepted by this crate's own tests: a
/// regression guard well above the RFC's threshold of 0.
const Q_GUARD: f64 = 30.0;

/// The conformance criterion of RFC 6716 §6, with the test vectors and
/// two reference sets of RFC 8251 §11: for every vector, every output rate
/// (8, 12, 16, 24, 48 kHz) and both channel counts, the decoder's final
/// range must match on every packet and its output must reach opus_compare
/// quality Q >= 0 against the reference output. RFC 8251 accepts either
/// set; the normal set (`.dec`) is decoded with phase inversion, the `m`
/// set (`m.dec`) without, and both are required here. Every Q must also
/// clear [`Q_GUARD`].
#[test]
fn test_vectors_conformance() {
    let Some(d) = dir() else {
        eprintln!("SKIPPED: test vectors not found (set OPUS_TESTVECTORS)");
        return;
    };
    let rates = [8000u32, 12000, 16000, 24000, 48000];
    let results: Vec<(Vec<String>, Vec<String>)> = std::thread::scope(|s| {
        let handles: Vec<_> = (1..=12usize)
            .map(|v| {
                let d = d.clone();
                s.spawn(move || {
                    let packets = read_bit(&d.join(format!("testvector{v:02}.bit")));
                    let mut lines = Vec::new();
                    let mut failures = Vec::new();
                    for (set, file, no_inv) in [("normal", format!("testvector{v:02}.dec"), false), ("m", format!("testvector{v:02}m.dec"), true)] {
                        let reference: Vec<f32> = read_pcm(&d.join(&file)).iter().map(|x| x * 32768.0).collect();
                        for channels in [1usize, 2] {
                            let spec = opus_compare::reference_spectrum(&reference, channels);
                            let mut qs = Vec::new();
                            for rate in rates {
                                let (pcm, bad, first) = decode(&packets, rate, channels, no_inv);
                                assert_eq!(bad, 0, "vector {v} at {rate} Hz x{channels}: final range differs (first at packet {first:?})");
                                assert_eq!(pcm.len() / channels * (48000 / rate) as usize, reference.len() / 2, "vector {v} at {rate} Hz: length");
                                let q = opus_compare::quality(&spec, &opus_compare::to_pcm16_units(&pcm), rate, channels).expect("length");
                                if q < Q_GUARD {
                                    failures.push(format!("vector {v:02} set {set} {channels} ch {rate} Hz: Q = {q:.1}"));
                                }
                                qs.push(format!("{q:>6.1}"));
                            }
                            lines.push(format!("vector {v:02} {set:>6} {channels} ch  Q at 8/12/16/24/48 kHz: {}", qs.join(" ")));
                        }
                    }
                    (lines, failures)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut failures = Vec::new();
    for (lines, f) in results {
        for l in lines {
            eprintln!("{l}");
        }
        failures.extend(f);
    }
    assert!(
        failures.is_empty(),
        "opus_compare quality below {Q_GUARD}:\n{}",
        failures.join("\n")
    );
}
