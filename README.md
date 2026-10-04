# rivet-opus

[![CI](https://github.com/safewords/rivet-opus/actions/workflows/ci.yml/badge.svg)](https://github.com/safewords/rivet-opus/actions/workflows/ci.yml)

An **Opus** encoder and decoder in Rust: no C, no system libraries, no
build script, nothing to install on a build host. Written from RFC 6716 as
updated by RFC 8251 (and RFC 7845 for the identification header and channel
mappings), not translated from any implementation. The decoder reproduces
the reference decoder's final range-coder state on every packet of all
twelve official test vectors (the figures are [below](#how-it-is-checked)).

Written for the **[rivet](https://github.com/safewords/rivet)**
transcoder, where it replaces libopus on both sides: the encoder behind
`audio=opus`, and the decoder for Opus sources in MP4, Matroska and Ogg.

Published as `rivet-opus`; **imported as `opus`** (`use opus::…`). One
dependency (`thiserror`), no features, no build script.

```toml
[dependencies]
opus = { package = "rivet-opus", git = "https://github.com/safewords/rivet-opus", branch = "develop" }
```

## What it decodes

Everything RFC 6716 defines:

| | |
|---|---|
| **Modes** | SILK (NB, MB, WB), hybrid (SWB, FB), CELT (NB, WB, SWB, FB) |
| **Frames** | 2.5, 5, 10, 20, 40 and 60 ms; packet codes 0–3 with padding; every framing rule [R1]–[R7] enforced |
| **Channels** | mono and stereo streams, decoded to mono or stereo; multistream with mapping families 0, 1 (mono to 7.1, Vorbis order) and 255 |
| **Output** | interleaved `f32` at 8, 12, 16, 24 or 48 kHz |
| **Transitions** | SILK⇄hybrid⇄CELT with and without redundancy (§4.5), bandwidth changes, mono⇄stereo |
| **Loss** | concealment (CELT: pitch-periodic extension through the MDCT overlap; SILK: LPC extrapolation of the last excitation), LBRR forward error correction (`Decoder::decode_fec`) |
| **RFC 8251** | all of it, including the optional no-phase-inversion decode (`set_phase_inversion_disabled`, on by default for mono output) |

Pre-skip is the caller's business (the decoder returns the stream from its
first sample, like libopus); `OpusHead::pre_skip` carries it, and
`MultistreamDecoder::from_head` applies the head's output gain.

## What it encodes

| | |
|---|---|
| **Modes** | CELT (music; the default for `Application::Audio`), SILK and hybrid (speech; `Application::Voip`), or forced |
| **Rates** | 6 to 510 kb/s, CBR (exact, with Opus padding) or VBR |
| **Frames** | 2.5 to 60 ms (CELT and hybrid above 20 ms as multi-frame packets) |
| **Input** | interleaved `f32` at 8, 12, 16, 24 or 48 kHz; mono and stereo, 3–8 channels through `MultistreamEncoder` (family 1) |
| **FEC** | LBRR frames for the next packet (`EncoderConfig::fec`) |
| **Delay** | 312 samples at 48 kHz (6.5 ms) at 48 kHz input — the `OpusHead` pre-skip |

The CELT encoder: pre-emphasis, MDCT with energy-rise transient detection
and short blocks, coarse/fine energy quantization with intra decisions, band
boost (§5.3.4.1), allocation trim, band skipping, dual/intensity stereo
decisions (§5.3.5), tonality-based spreading, greedy PVQ search. The SILK
encoder: autocorrelation LPC with lookahead, LSF conversion and a two-stage
quantizer that reconstructs exactly as the decoder does, residual-domain
pitch with contour refinement, 5-tap LTP with rate-distortion codebook
choice, adaptive mid/side with the decoder's prediction, and a
noise-feedback quantizer that runs the decoder's own synthesis, so encoder
and decoder cannot drift; a rate loop scales the gains to the budget.

**Not implemented** (all optional for an encoder): the CELT pitch
pre-filter, RDO TF analysis, the SILK delayed-decision quantizer and
pre-filter, DTX, and redundancy frames at mode switches (a stream whose mode
changes is concealed across the switch by the decoder, as RFC 6716 Figure 19
allows). Opus Custom is not implemented.

## How it is checked

**The official test vectors** (opus_testvectors-rfc8251, used as data;
`tests/vectors.rs`, run by CI), checked by RFC 6716 §6's own conformance
criterion: on every one of the 12 vectors, at 8, 12, 16, 24 and 48 kHz, mono
and stereo, against both RFC 8251 reference sets, the decoder reproduces the
reference's **final range-coder state on every packet** and passes the
**`opus_compare` quality metric** (re-implemented in `tests/opus_compare/`
from its specification in RFC 6716 Appendix A). All 240 comparisons pass;
the lowest quality is 37.9 (SILK NB at 48 kHz; the threshold is 0, and 100
means identical output). The CELT-only vectors are float rounding away from
the reference (99.9–106 dB SNR, quality 98–99.9 at 48 kHz). The full tables,
waveform SNRs and the resampler's measured response are in
[docs/VALIDATION.md](docs/VALIDATION.md).

| vector | content | final range | Q at 48 kHz (stereo) | lowest Q, any rate/channels |
|---|---|---|---|---|
| 01 | CELT stereo | 2147/2147 | 99.9 | 91.4 |
| 02 | SILK NB | 1185/1185 | 38.2 | 37.9 |
| 03 | SILK MB | 998/998 | 60.5 | 60.5 |
| 04 | SILK WB | 1265/1265 | 73.9 | 73.9 |
| 05 | hybrid SWB | 2037/2037 | 55.9 | 43.0 |
| 06 | hybrid FB | 1876/1876 | 66.2 | 53.5 |
| 07 | CELT, all sizes and bandwidths | 4186/4186 | 99.9 | 76.0 |
| 08 | mixed SILK/CELT | 1247/1247 | 94.8 | 70.0 |
| 09 | mixed SILK/CELT | 1337/1337 | 87.5 | 75.5 |
| 10 | hybrid/CELT switching | 1912/1912 | 87.7 | 86.2 |
| 11 | CELT stereo, all packet codes | 553/553 | 99.8 | 91.0 |
| 12 | SILK bandwidth switching | 1332/1332 | 43.8 | 43.3 |

The SILK resampler is this crate's own design (a delay-constrained
least-squares FIR with exactly the RFC's Table 54 delay; RFC 6716 §4.2.9
leaves the resampler to the decoder), so SILK and hybrid output differs from
the reference in phase and fine detail while meeting the criterion.

**Round trips** (`tests/roundtrip.rs`): this encoder, this decoder, 2 s of
synthetic music (harmonic tones with vibrato) or speech (a glottal pulse
train through formant resonators). Every packet is parsed against RFC 6716
§3, encoder and decoder final ranges must agree on every packet, and SNR is
measured after the encoder's lookahead. CBR rates are exact. A selection at
20 ms:

| mode | signal | 1 ch | 2 ch |
|---|---|---|---|
| CELT 16 / 24 kb/s | music | 14.7 dB | 10.5 dB (24 kb/s) |
| CELT 64 kb/s | music | 22.7 dB | 15.6 dB |
| CELT 128 kb/s | music | 27.5 dB | 20.8 dB |
| CELT 510 kb/s | music | 32.0 dB | 32.0 dB |
| SILK 8 kb/s per channel | speech | 7.7 dB | 17.0 dB |
| SILK 16 kb/s per channel | speech | 15.4 dB | 26.2 dB |
| SILK 32 kb/s per channel | speech | 27.2 dB | 33.2 dB |
| hybrid 24 kb/s per channel | speech | 18.4 dB | 29.2 dB |
| hybrid 32 kb/s per channel | speech | 21.3 dB | 31.3 dB |

(The stereo speech signal is one voice at two levels, so its side channel is nearly free. The 32 dB ceiling at 510 kb/s is the test signal's noise above 20 kHz, which
Opus does not code.) The full matrix covers 2.5–60 ms, 8–510 kb/s, mono and
stereo, CBR and VBR, in all three modes, plus 8–24 kHz input. `tests/loss.rs`
drops packets in every mode (concealment stays bounded) and recovers a lost
SILK packet from the next one's LBRR data (26 dB mono, 33 dB stereo, against
6 dB for concealment). The multistream round trip keeps each of 1–8
channels' tones in its own channel. `tests/final_range.rs` checks encoder
and decoder final ranges on every packet across SILK, hybrid and CELT, mono
and stereo, CBR, VBR and LBRR, 6–510 kb/s and every frame size, on speech,
music and a hard-panned pair. Hybrid stereo keeps hard-panned tone pairs
(440/660 Hz to 3/5 kHz) 24–33 dB apart at 32 kb/s, 38–43 dB at 40 kb/s and 57–61 dB at
64 kb/s (`tests/stereo_separation.rs`, `hybrid_separation_table`).

**Unit tests**: range coder round trips with matching encoder/decoder state
after every symbol; Laplace coding; the PVQ codebook enumerated exhaustively
for N, K ≤ 6 (a bijection with indices `0..V(N,K)`); the FFT and MDCT
against their defining sums and TDAC reconstruction; the pulse cache;
power-complementary window; SILK tables, NLSF weights, stabilisation and
filter stability; packet rules [R1]–[R7] and self-delimited framing.

### NEON on ARM hardware

CI runs on x86-64 Linux only, so the NEON (aarch64) code paths are not tested
there. They are verified by hand on ARM hardware (an aarch64 Linux machine,
or Apple silicon) after a change to them and before a release:

```sh
cargo test --release
cargo test --release --features force-scalar
```

The second run compiles the vector paths out; both must pass unchanged.

## Performance

On a Ryzen 9 9950X (one thread, `cargo run --release --example bench --
<vector dir>`): decoding all twelve RFC 8251 vectors at 48 kHz stereo runs
about 600x real time; encoding stereo music about 290x (CELT 128 kb/s),
420x (CELT 64 kb/s), 45x (hybrid 32 kb/s) and 90x (SILK 16 kb/s).

The few vector kernels (`src/simd.rs`: AVX2 on x86-64, NEON on AArch64,
chosen at run time) multiply and add separately in eight lanes and reduce
in a fixed order, exactly as their portable versions do, so the decoded
output and the encoded packets are the same on every CPU. The
`force-scalar` feature compiles the vector paths out; CI runs the tests
both ways.

## Provenance and licensing

Written from the RFCs' text. libopus, FFmpeg and other implementations were
not read. Where RFC 6716 makes its attached reference code (Appendix A) the
normative definition, that code was used **as specification only**: one
author described the CELT layer's normative behaviour in prose
([docs/CELT_SPEC.md](docs/CELT_SPEC.md)) and another, who never saw that
code, implemented it from the description; the conformance metric was
likewise re-implemented from its description. No filter or table is fitted
to any decoder's output. [docs/PROVENANCE.md](docs/PROVENANCE.md) records
every source, every table and how each was obtained.

**Patents.** Opus is covered by patents whose holders have made royalty-free
licensing commitments to the IETF. Nothing here is a licence to any patent.

## Using it

```rust
use opus::{Decoder, Encoder, EncoderConfig, MultistreamDecoder, OpusHead};

// Encoding: interleaved f32 at 48 kHz, one packet per frame.
let mut enc = Encoder::new(EncoderConfig { channels: 2, bitrate: 96_000, ..EncoderConfig::default() })?;
let head = OpusHead::new(2, enc.lookahead() as u16, 48_000)?; // dOps / OpusHead
let packet = enc.encode(&vec![0.0f32; enc.frame_samples() * 2])?;

// Decoding a stream described by its OpusHead (any family).
let mut dec = MultistreamDecoder::from_head(&head, 48_000)?;
let pcm = dec.decode(Some(&packet))?; // None conceals a lost packet
```

## License

Open Encoding Attribution License v1.0 — a source-available (not OSI open-source)
license, royalty-free, with a commercial-attribution requirement. See
[LICENSE.md](LICENSE.md) and [NOTICE](NOTICE).
