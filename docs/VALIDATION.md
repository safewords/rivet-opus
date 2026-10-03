# Validation

How this crate is checked, and the measured results. Every figure here is
produced by the tests named; rerun them to reproduce.

## The conformance criterion of RFC 6716 §6

RFC 6716 §6 defines a compliant decoder as one that, for each official test
vector, **reproduces the reference decoder's final range-coder state on
every packet** and produces output **within the thresholds of the
`opus_compare` quality metric** against the reference output, "for each
output sampling rate and channel count supported". RFC 8251 §11 replaced the
vectors and provides two reference outputs (`.dec` with the CELT intensity
phase inversion, `m.dec` without); passing either set is compliant.

The metric is given in the RFC only by `opus_compare.c` (Appendix A); it is
re-implemented in `tests/opus_compare/mod.rs` from a description of that
specification (module documentation). It scores 100 for identical output
and 0 at the pass threshold.

`tests/vectors.rs::test_vectors_conformance` decodes all 12 vectors of the
RFC 8251 set at 8, 12, 16, 24 and 48 kHz, mono and stereo, against both
reference sets (the normal set with phase inversion, the `m` set without),
and requires the reference final range on every packet and Q ≥ 30 (a
regression guard well above the RFC's 0). Run it with
`OPUS_TESTVECTORS=<dir> cargo test --release --test vectors -- --nocapture`.

**Result: every packet of every vector has the reference final range at
every rate and channel count, and every one of the 240 comparisons passes;
the lowest quality is 37.9** (vector 02, SILK NB, at 48 kHz). RFC 6716 §6.1
recommends Q above 90 at 48 kHz and notes that "as low as 50" is normal at
other rates "because of harmless mismatch with the delay and phase of the
internal sampling rate conversion". This decoder exceeds 90 at 48 kHz on
the CELT-only and CELT-dominated vectors (01, 07, 08 and 11, at 94.8–99.9) but
not on the SILK and hybrid ones, where it uses its own resampler and float
SILK synthesis (below).

Quality Q per vector, reference set, output channels and rate:

| vector | set | ch | 8 kHz | 12 kHz | 16 kHz | 24 kHz | 48 kHz |
|---|---|---|---|---|---|---|---|
| 01 | normal | 1 | 91.4 | 95.6 | 95.2 | 96.7 | 98.4 |
| 01 | normal | 2 | 93.3 | 96.0 | 94.7 | 96.8 | 99.9 |
| 01 | m | 1 | 91.4 | 95.6 | 95.3 | 97.4 | 98.2 |
| 01 | m | 2 | 93.3 | 96.0 | 94.7 | 96.8 | 99.9 |
| 02 | normal | 1 | 59.7 | 42.2 | 50.3 | 42.0 | 37.9 |
| 02 | normal | 2 | 57.6 | 42.6 | 48.9 | 42.2 | 38.2 |
| 02 | m | 1 | 59.7 | 42.2 | 50.3 | 42.0 | 37.9 |
| 02 | m | 2 | 57.6 | 42.6 | 48.9 | 42.2 | 38.2 |
| 03 | normal | 1 | 74.3 | 79.3 | 61.5 | 61.4 | 61.3 |
| 03 | normal | 2 | 73.8 | 74.7 | 60.7 | 60.6 | 60.5 |
| 03 | m | 1 | 74.3 | 79.3 | 61.5 | 61.4 | 61.3 |
| 03 | m | 2 | 73.8 | 74.7 | 60.7 | 60.6 | 60.5 |
| 04 | normal | 1 | 75.9 | 84.3 | 85.1 | 75.0 | 74.9 |
| 04 | normal | 2 | 76.0 | 84.9 | 76.2 | 74.0 | 73.9 |
| 04 | m | 1 | 75.9 | 84.3 | 85.1 | 75.0 | 74.9 |
| 04 | m | 2 | 76.0 | 84.9 | 76.2 | 74.0 | 73.9 |
| 05 | normal | 1 | 75.7 | 77.5 | 77.7 | 56.2 | 56.1 |
| 05 | normal | 2 | 75.8 | 77.9 | 77.9 | 56.0 | 55.9 |
| 05 | m | 1 | 75.7 | 77.5 | 77.6 | 44.9 | 44.8 |
| 05 | m | 2 | 75.8 | 77.8 | 77.8 | 43.0 | 43.1 |
| 06 | normal | 1 | 68.1 | 69.1 | 68.5 | 66.2 | 66.2 |
| 06 | normal | 2 | 68.2 | 69.3 | 68.7 | 66.3 | 66.2 |
| 06 | m | 1 | 68.1 | 69.1 | 68.6 | 53.5 | 53.6 |
| 06 | m | 2 | 68.2 | 69.3 | 68.7 | 66.6 | 66.5 |
| 07 | normal | 1 | 76.0 | 85.1 | 89.9 | 94.0 | 98.1 |
| 07 | normal | 2 | 76.0 | 85.1 | 89.8 | 94.0 | 99.9 |
| 07 | m | 1 | 76.0 | 85.1 | 89.7 | 94.0 | 98.1 |
| 07 | m | 2 | 76.0 | 85.1 | 89.8 | 94.0 | 99.9 |
| 08 | normal | 1 | 70.0 | 81.6 | 86.2 | 92.4 | 96.0 |
| 08 | normal | 2 | 70.4 | 81.8 | 87.6 | 93.6 | 94.8 |
| 08 | m | 1 | 70.0 | 81.6 | 87.5 | 93.6 | 96.0 |
| 08 | m | 2 | 70.3 | 81.8 | 87.6 | 93.6 | 94.8 |
| 09 | normal | 1 | 75.5 | 84.8 | 78.8 | 87.8 | 88.1 |
| 09 | normal | 2 | 76.2 | 85.1 | 78.9 | 87.4 | 87.5 |
| 09 | m | 1 | 75.5 | 84.8 | 78.8 | 87.8 | 88.1 |
| 09 | m | 2 | 76.2 | 85.1 | 78.9 | 87.4 | 87.5 |
| 10 | normal | 1 | 86.2 | 92.4 | 88.1 | 89.7 | 89.8 |
| 10 | normal | 2 | 86.2 | 87.6 | 87.6 | 87.7 | 87.7 |
| 10 | m | 1 | 86.2 | 92.4 | 88.1 | 89.7 | 89.8 |
| 10 | m | 2 | 86.2 | 87.6 | 87.6 | 87.7 | 87.7 |
| 11 | normal | 1 | 91.0 | 94.9 | 96.3 | 98.0 | 98.4 |
| 11 | normal | 2 | 93.0 | 95.8 | 96.3 | 98.1 | 99.8 |
| 11 | m | 1 | 91.0 | 94.9 | 96.4 | 98.0 | 98.5 |
| 11 | m | 2 | 92.9 | 95.8 | 96.3 | 98.0 | 99.8 |
| 12 | normal | 1 | 43.3 | 47.2 | 51.3 | 47.0 | 43.6 |
| 12 | normal | 2 | 43.5 | 47.4 | 51.5 | 47.2 | 43.8 |
| 12 | m | 1 | 43.3 | 47.2 | 51.3 | 47.0 | 43.6 |
| 12 | m | 2 | 43.5 | 47.4 | 51.5 | 47.2 | 43.8 |

### Waveform SNR at 48 kHz

`tests/vectors.rs::test_vectors_48k` also reports the SNR of the 48 kHz
output against the reference (16-bit units for the maximum error). The
CELT-only vectors are float rounding away from the reference and are held
to SNR floors (01 and 11 above 100 dB, 07 above 95 dB). Vectors with SILK
content are not: the RFC leaves the resampler free (§4.2.9), so their
waveform differs from the reference by the two resamplers' different phase
responses while the conformance metric above is met.

| vector | content | final range | stereo SNR | max err | no-inversion stereo vs `m` | mono vs `m` downmix |
|---|---|---|---|---|---|---|
| 01 | CELT stereo | 2147/2147 | 106.27 dB | 1 | 106.29 dB | 76.83 dB |
| 02 | SILK NB | 1185/1185 | 15.89 dB | 2735 | 15.89 dB | 16.17 dB |
| 03 | SILK MB | 998/998 | 19.37 dB | 1317 | 19.37 dB | 20.92 dB |
| 04 | SILK WB | 1265/1265 | 18.78 dB | 3444 | 18.78 dB | 20.53 dB |
| 05 | hybrid SWB | 2037/2037 | 20.87 dB | 1522 | 20.87 dB | 23.71 dB |
| 06 | hybrid FB | 1876/1876 | 21.41 dB | 1308 | 21.41 dB | 23.63 dB |
| 07 | CELT, all sizes and bandwidths | 4186/4186 | 99.87 dB | 1 | 99.82 dB | 67.42 dB |
| 08 | mixed SILK/CELT | 1247/1247 | 81.60 dB | 4 | 81.60 dB | 65.19 dB |
| 09 | mixed SILK/CELT | 1337/1337 | 76.87 dB | 10 | 76.87 dB | 66.30 dB |
| 10 | hybrid/CELT switching | 1912/1912 | 52.67 dB | 102 | 52.67 dB | 53.38 dB |
| 11 | CELT stereo, all packet codes | 553/553 | 106.17 dB | 1 | 106.27 dB | 79.79 dB |
| 12 | SILK bandwidth switching | 1332/1332 | 19.78 dB | 2665 | 19.78 dB | 19.78 dB |

(The earlier fitted resampler gave 44–48 dB on vectors 02–06 because it was
fitted to these same outputs; that is not a conformance criterion.)

## The SILK resampler

`src/resample.rs` (design in its module documentation). Its unit tests
measure the response of every converter the decoder and encoder build, from
the prototype filter reassembled out of the polyphase branches:
`frequency_response` asserts passband ripple < 0.8 dB, group-delay error
< 0.04 ms over 90 % of the passband, −6 ± 1.5 dB at the lower Nyquist
frequency, > 30 dB rejection in the near stopband and > 60 dB beyond;
`streaming_tone` checks block-wise streaming and the delay of a tone. The
table below is printed by
`cargo test --lib print_response_table -- --ignored --nocapture`.
"Delay error" is the largest deviation of the group delay from the target
over 90 % of the passband; the near stopband is `[fc + Δ, 0.75·lo]`, the far
stopband beyond (`lo` the lower rate, `fc = lo/2`, `Δ = 0.05·lo`). Equal rates
are a delay line of the nearest whole number of samples (§4.2.9 allows the
rounding), hence their delay error of up to half a sample.

| in → out | delay (ms) | taps (per branch) | passband ripple | delay error | at fc | near stopband | far stopband |
|---|---|---|---|---|---|---|---|
| 8 → 8 kHz | 0.538 | 5 (1) | ±0.00 dB | 0.038 ms | 0.0 dB | — | — |
| 8 → 12 kHz | 0.538 | 150 (50) | ±0.58 dB | 0.029 ms | -5.3 dB | 34.0 dB | 66.7 dB |
| 8 → 16 kHz | 0.538 | 100 (50) | ±0.44 dB | 0.033 ms | -5.8 dB | 35.6 dB | 73.7 dB |
| 8 → 24 kHz | 0.538 | 150 (50) | ±0.58 dB | 0.029 ms | -5.3 dB | 34.0 dB | 66.7 dB |
| 8 → 48 kHz | 0.538 | 300 (50) | ±0.76 dB | 0.036 ms | -5.2 dB | 33.4 dB | 65.3 dB |
| 12 → 8 kHz | 0.692 | 150 (75) | ±0.69 dB | 0.031 ms | -6.2 dB | 37.0 dB | 72.5 dB |
| 12 → 12 kHz | 0.692 | 9 (1) | ±0.00 dB | 0.025 ms | 0.0 dB | — | — |
| 12 → 16 kHz | 0.692 | 200 (50) | ±0.46 dB | 0.017 ms | -6.3 dB | 40.3 dB | 75.1 dB |
| 12 → 24 kHz | 0.692 | 100 (50) | ±0.35 dB | 0.014 ms | -6.4 dB | 41.9 dB | 80.7 dB |
| 12 → 48 kHz | 0.692 | 200 (50) | ±0.46 dB | 0.017 ms | -6.3 dB | 40.3 dB | 75.1 dB |
| 16 → 8 kHz | 0.706 | 100 (100) | ±0.73 dB | 0.032 ms | -6.3 dB | 37.3 dB | 73.7 dB |
| 16 → 12 kHz | 0.706 | 200 (67) | ±0.39 dB | 0.015 ms | -6.3 dB | 40.8 dB | 76.2 dB |
| 16 → 16 kHz | 0.706 | 12 (1) | ±0.00 dB | 0.019 ms | 0.0 dB | — | — |
| 16 → 24 kHz | 0.706 | 150 (50) | ±0.13 dB | 0.005 ms | -6.3 dB | 45.4 dB | 79.4 dB |
| 16 → 48 kHz | 0.706 | 150 (50) | ±0.13 dB | 0.005 ms | -6.3 dB | 45.4 dB | 79.4 dB |
| 8 → 48 kHz | 1 | 300 (50) | ±0.59 dB | 0.029 ms | -6.6 dB | 40.5 dB | 77.1 dB |
| 12 → 48 kHz | 1 | 200 (50) | ±0.12 dB | 0.008 ms | -6.4 dB | 46.6 dB | 84.8 dB |
| 16 → 48 kHz | 1 | 150 (50) | ±0.02 dB | 0.002 ms | -6.2 dB | 52.3 dB | 93.7 dB |
| 24 → 48 kHz | 1 | 100 (50) | ±0.02 dB | 0.001 ms | -5.9 dB | 55.2 dB | 100.5 dB |
| 48 → 8 kHz | 1.337 | 300 (300) | ±0.16 dB | 0.013 ms | -6.4 dB | 44.5 dB | 82.3 dB |
| 48 → 12 kHz | 1.2247 | 200 (200) | ±0.04 dB | 0.004 ms | -6.2 dB | 51.0 dB | 94.1 dB |
| 48 → 16 kHz | 1.2315 | 150 (150) | ±0.04 dB | 0.001 ms | -6.0 dB | 54.3 dB | 104.2 dB |

The first 15 rows are the decoder's (SILK internal rate to each output
rate, Table 54 delay); the last 7 the encoder's (input to 48 kHz with 1 ms
delay, 48 kHz to the SILK rate with the rest of its lookahead). The NB rows
are the hardest: the Table 54 budget is 0.538 ms, 4.3 samples at 8 kHz,
which bounds how sharp a filter with that delay can be.

## The CELT layer

`src/celt/` is a clean-room implementation of `docs/CELT_SPEC.md` (see
`docs/PROVENANCE.md`). Beyond the vectors above, its unit tests check every
table derived in code against the data extracted from RFC 6716 Appendix A
(`log2_frac_tables_match_appendix_a`, `pulse_cache_matches_appendix_a`,
`caps_match_appendix_a`, `small_tables_match_appendix_a`,
`float_tables_are_scaled_fixed_point`, `window_is_the_rfc_formula`), plus
Laplace round trips, bits-to-pulses, the conservative log2, the TF
interleave inverse, the spreading-rotation inverse, the integer square root
and the bit-exact cosine/log-tangent ranges.

Points where the implementer had to choose a reading of the specification
(none affects the range coder; all were settled by the vectors):

- CELT_SPEC §10.4 states the reference inverse MDCT as the textbook one
  times `1 + (π/8n)²`. Applying that factor with this crate's textbook IMDCT
  lowers vector 01 from 106.3 to 96.0 dB, so it is not applied (gain 1); the
  stated relation evidently already includes it in the reference's own
  convention.
- §8.1 step 7, the collapse-mask range of the folding source, is read as
  inclusive (output samples only; no measurable difference).
- The encoder-side rotation is the exact inverse of the decoder's (passes
  undone in reverse order) rather than the spec's description of the
  reference encoder; the encoder is non-normative.

## Round trips, loss, unit tests

`tests/roundtrip.rs` encodes and decodes synthetic music and speech across
2.5–60 ms, 8–510 kb/s, mono and stereo, CBR and VBR, all three modes and
8–24 kHz input; every packet is parsed against RFC 6716 §3 and encoder and
decoder final ranges must agree on every packet. `tests/loss.rs` drops
packets in every mode and recovers a lost SILK packet from LBRR data.
Unit tests cover the range coder, the PVQ codebook (exhaustive for
N, K ≤ 6), the FFT and MDCT against their definitions, the SILK tables and
filters, packet rules [R1]–[R7] and self-delimited framing, the resampler and
the conformance metric (`identity_and_monotonicity`).
