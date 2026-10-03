# Provenance

Every source this crate was written from, and what came from where.

## Specifications

| source | used for |
|---|---|
| RFC 6716 §3 | packet framing, TOC, rules [R1]–[R7] (`src/packet.rs`) |
| RFC 6716 §4.1, §5.1 | range decoder and encoder, `tell`/`tell_frac` (`src/range.rs`) |
| RFC 6716 §4.2 | SILK decoder: header and LBRR flags, every parameter, the NLSF fixed-point reconstruction, stabilisation, NLSF→LPC with range and prediction-gain limiting, excitation, LTP/LPC synthesis (as the RFC's floating-point description), stereo unmixing (`src/silk/`) |
| RFC 6716 §4.3 | CELT decoder (`src/celt/`) |
| RFC 6716 §4.4, §4.5 | concealment, transitions, redundancy (`src/decoder.rs`) |
| RFC 6716 §5.2, §5.3 | the encoders' structure (the encoders themselves are an original design) |
| RFC 6716 Appendix B | self-delimited framing (`src/packet.rs`) |
| RFC 8251 | stereo state reset (§3), padding parsing (§4), inverse-gain and LSF overflow fixes (§6, §7), the band-energy cap (§8), hybrid folding (§9), the no-phase-inversion option (§10) |
| RFC 7845 §5.1 | `OpusHead`, mapping families 0, 1 and 255 (`src/multistream.rs`) |

The base64-encoded reference source of RFC 6716 Appendix A was not read.

## Tables

- **SILK** (Tables 4–53 of RFC 6716): extracted from the RFC's plain text by
  a script (column parsing of each table, PDFs converted to inverse CDFs)
  into `src/silk/tables.rs`. Unit tests check that every PDF sums to 256,
  that the codebooks rise and that their weights fall in the RFC's stated
  range.
- **CELT**, from the RFC text: band edges (Table 55), static allocation
  (Table 57, extracted by script), trim/spread/tapset PDFs (Tables 56, 58,
  59), TF adjustments (Tables 60–63), post-filter taps and de-emphasis
  (§4.3.7).
- **CELT**, named by the RFC's prose but given there only as data of the
  reference implementation (`e_prob_model`, `cache_caps50`, the energy
  prediction coefficients, the band means, `LOG2_FRAC_TABLE`, the pulse
  cache, the band-split constants): written from the author's knowledge of
  the specification, or derived from their mathematical definition (the
  pulse cache from V(N,K), `LOG2_FRAC_TABLE` as `ceil(8·log2(i+1))`, the
  window from §4.3.7). Any error in these desynchronises the range coder,
  which the test vectors' per-packet final-range check would catch; all
  12 vectors match on every packet.
- **The SILK-to-48 kHz resamplers** (`src/resample_fit.rs`): the RFC makes
  the resampler non-normative, fixing only its delay (Table 54). These
  polyphase filters were identified by least squares
  (`tests/fit_resampler.rs`) from test vectors 02, 03 and 04 — this crate's
  SILK output at the internal rate against the reference decoder's 48 kHz
  output — so the SILK path follows the reference's delay and response.
  Every other rate pair uses a windowed-sinc design with the Table 54 delay.

## Validation data

The official test vectors, `opus_testvectors-rfc8251.tar.gz` from
opus-codec.org (SHA-1s listed in RFC 8251 §11), are downloaded at test time
(CI does this) and used only as data.
