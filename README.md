# rivet-opus

[![CI](https://github.com/rivet-transcoder/rivet-opus/actions/workflows/ci.yml/badge.svg)](https://github.com/rivet-transcoder/rivet-opus/actions/workflows/ci.yml)

An **Opus** encoder and decoder in Rust: no C, no system libraries, no
build script, nothing to install on a build host. Written from RFC 6716 as
updated by RFC 8251 (and RFC 7845 for the identification header and channel
mappings), not translated from any implementation. The decoder reproduces
the reference decoder's final range-coder state on every packet of all
twelve official test vectors (the figures are [below](#how-it-is-checked)).

Written for the **[rivet](https://github.com/rivet-transcoder/rivet)**
transcoder, where it replaces libopus on both sides: the encoder behind
`audio=opus`, and the decoder for Opus sources in MP4, Matroska and Ogg.

Published as `rivet-opus`; **imported as `opus`** (`use opus::…`). One
dependency (`thiserror`), no features, no build script.

```toml
[dependencies]
opus = { package = "rivet-opus", git = "https://github.com/rivet-transcoder/rivet-opus", branch = "develop" }
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
`tests/vectors.rs`, run by CI). RFC 6716 §6 defines compliance through
`opus_compare`, a quality metric given only as source code, so it is not
used here; instead each vector is reported by its **final range-coder
state** (which a compliant decoder must reproduce exactly) and by **SNR and
largest error** against the reference output. 48 kHz:

| vector | content | final range | stereo SNR | max err | no-inversion stereo vs `m` | mono vs `m` downmix |
|---|---|---|---|---|---|---|
| 01 | CELT stereo | 2147/2147 | 106.30 dB | 1 | 106.36 dB | 76.83 dB |
| 02 | SILK NB | 1185/1185 | 48.22 dB | 31 | 48.22 dB | 48.54 dB |
| 03 | SILK MB | 998/998 | 47.39 dB | 43 | 47.39 dB | 47.99 dB |
| 04 | SILK WB | 1265/1265 | 44.18 dB | 66 | 44.18 dB | 44.91 dB |
| 05 | hybrid SWB | 2037/2037 | 43.13 dB | 54 | 43.13 dB | 43.65 dB |
| 06 | hybrid FB | 1876/1876 | 42.49 dB | 72 | 42.49 dB | 42.87 dB |
| 07 | CELT, all sizes and bandwidths, mono/stereo | 4186/4186 | 99.82 dB | 1 | 99.66 dB | 67.42 dB |
| 08 | mixed SILK/CELT | 1247/1247 | 84.15 dB | 2 | 84.15 dB | 65.24 dB |
| 09 | mixed SILK/CELT | 1337/1337 | 83.96 dB | 2 | 83.96 dB | 66.63 dB |
| 10 | hybrid/CELT switching | 1912/1912 | 58.94 dB | 67 | 58.94 dB | 62.23 dB |
| 11 | CELT stereo, all packet codes | 553/553 | 106.17 dB | 1 | 106.22 dB | 79.79 dB |
| 12 | SILK bandwidth switching | 1332/1332 | 45.56 dB | 45 | 45.56 dB | 45.56 dB |

Max error is in 16-bit units. The CELT vectors differ from the reference by
float rounding. The SILK and hybrid ones carry the difference between this
float synthesis and the reference's fixed-point one and its resampler: the
SILK-to-48 kHz resamplers are least-squares fits to the reference's output on
vectors 02–04 (the RFC leaves the resampler non-normative), so vectors 02–04
are in-sample for that fit and 05, 06, 10 and 12 are not. The vectors also
decode at 8, 12, 16 and 24 kHz, mono and stereo, with the reference's final
range on every packet (no reference PCM exists at those rates).

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
| hybrid 24 kb/s per channel | speech | 12.4 dB | 23.8 dB |
| hybrid 32 kb/s per channel | speech | 16.8 dB | 27.4 dB |

(The stereo speech signal is one voice at two levels, so its side channel is nearly free. The 32 dB ceiling at 510 kb/s is the test signal's noise above 20 kHz, which
Opus does not code.) The full matrix covers 2.5–60 ms, 8–510 kb/s, mono and
stereo, CBR and VBR, in all three modes, plus 8–24 kHz input. `tests/loss.rs`
drops packets in every mode (concealment stays bounded) and recovers a lost
SILK packet from the next one's LBRR data (26 dB mono, 33 dB stereo, against
6 dB for concealment). The multistream round trip keeps each of 1–8
channels' tones in its own channel.

**Unit tests**: range coder round trips with matching encoder/decoder state
after every symbol; Laplace coding; the PVQ codebook enumerated exhaustively
for N, K ≤ 6 (a bijection with indices `0..V(N,K)`); the FFT and MDCT
against their defining sums and TDAC reconstruction; the pulse cache;
power-complementary window; SILK tables, NLSF weights, stabilisation and
filter stability; packet rules [R1]–[R7] and self-delimited framing.

## Provenance and licensing

Written from the RFCs' text; **no Opus implementation's source was read** —
not the reference code attached to RFC 6716, libopus, FFmpeg's or any other.
[docs/PROVENANCE.md](docs/PROVENANCE.md) records every source and every
table.

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
