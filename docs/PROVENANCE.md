# Provenance

Every source this crate was written from, and what came from where.

## Specifications

| source | used for |
|---|---|
| RFC 6716 §3 | packet framing, TOC, rules [R1]–[R7] (`src/packet.rs`) |
| RFC 6716 §4.1, §5.1 | range decoder and encoder, `tell`/`tell_frac` (`src/range.rs`) |
| RFC 6716 §4.2 | SILK decoder: header and LBRR flags, every parameter, the NLSF fixed-point reconstruction, stabilisation, NLSF→LPC with range and prediction-gain limiting, excitation, LTP/LPC synthesis (as the RFC's floating-point description), stereo unmixing (`src/silk/`) |
| RFC 6716 §4.2.9 | the SILK output delay (Table 54); the resampler itself is this crate's own design (below) |
| RFC 6716 §4.3, Appendix A (as specification), RFC 8251 | CELT layer (`src/celt/`), through `docs/CELT_SPEC.md` (below) |
| RFC 6716 §4.4, §4.5 | concealment, transitions, redundancy (`src/decoder.rs`, `src/celt/decoder.rs`) |
| RFC 6716 §5.2, §5.3 | the encoders' structure; their analysis heuristics are this crate's own designs following that prose |
| RFC 6716 §6, Appendix A `opus_compare.c` (as specification) | the conformance quality metric, re-implemented for the tests (`tests/opus_compare/`) |
| RFC 6716 Appendix B | self-delimited framing (`src/packet.rs`) |
| RFC 8251 | stereo state reset (§3), padding parsing (§4), inverse-gain and LSF overflow fixes (§6, §7), the band-energy cap (§8), hybrid folding (§9), the no-phase-inversion option (§10), the two reference output sets (§11) |
| RFC 7845 §5.1 | `OpusHead`, mapping families 0, 1 and 255 (`src/multistream.rs`) |

No implementation of Opus other than RFC 6716 Appendix A was consulted, and
Appendix A only as described here. libopus and FFmpeg were not read or used,
and nothing in this repository (code, tests, data or tools) uses FFmpeg.

## How RFC 6716 Appendix A was used

RFC 6716 §6 makes the reference code of Appendix A normative: "should the
description contradict the source code of the reference implementation, the
latter shall take precedence", and it defines compliance through
`opus_compare.c`. Where the RFC prose defers to that code, it was used **as
specification only**, through a separation of roles:

1. **Describer.** One author read RFC 6716 §4.3, RFC 8251 and the Appendix A
   CELT sources and wrote `docs/CELT_SPEC.md`: a prose and mathematical
   description, with no code, that cites for every item whether it comes
   from the RFC prose, from Appendix A (file and function), from RFC 8251,
   or is derived mathematically, and lists the 12 places where prose and
   code disagree (§14 of that document; the code governs). The same
   author wrote `tools/appendix_a_tables.py`, which recovers Appendix A from
   the RFC text exactly as RFC 6716 §A.1 describes and prints its constant
   data tables (data only).
2. **Implementer.** A second author, who did not read Appendix A, any other
   Opus implementation, or the previous versions of the files being
   replaced, wrote `src/celt/{tables,mode,energy,rate,bands,decoder}.rs`
   from `docs/CELT_SPEC.md`, the RFC prose and the extracted data tables.
   Tables with a mathematical definition are derived in code and unit-tested
   against the extracted data; the rest are copied as data with their
   source named.

The conformance metric (`tests/opus_compare/mod.rs`) is described in its
module documentation from `opus_compare.c` and re-implemented, not
translated.

## Items previously written from memory, and their status now

The first version of the CELT layer was written partly from the author's
memory of the reference implementation rather than from the RFC text. Every
such piece and its current basis:

| item | basis now | how |
|---|---|---|
| `e_prob_model` (Laplace energy model) | Appendix A `quant_bands.c` (data) | extracted by script; copied. The remembered values were identical |
| energy prediction alpha/beta per LM | Appendix A `quant_bands.c` (data); intra beta = 4915/32768 is RFC §4.3.2.1 prose | extracted; copied. Identical |
| band means `eMeans` | Appendix A `quant_bands.c` (data) | extracted; copied; tested equal to the fixed-point table / 16. Identical |
| `cache_caps50` (band caps) | Appendix A `rate.c` procedure (RFC §4.3.3 names it) | derived in code from the pulse-cache computation; test: equals the extracted table. Remembered values were identical |
| pulse cache (`cache_index50`, `cache_bits50`) | V(N,K) (RFC §4.3.4.2) + Appendix A `rate.c` construction | derived; test: all 392 bytes equal the extracted table |
| `LOG2_FRAC_TABLE`, logN | Appendix A `rate.c` / `entcode.c` (conservative log2) | derived; test: equal to the extracted tables |
| band-split constants (`exp2_table8`, theta offsets, `compute_qn`) | Appendix A `bands.c` | derived (`floor(16384·2^(k/8))`), tested; procedure re-implemented from CELT_SPEC §8.2.4 |
| sequency order (`ordery`), bit (de)interleave masks | Appendix A `bands.c` | derived, tested equal |
| low-overlap window | RFC §4.3.7 formula; Appendix A `window120` | table copied (5 of 120 entries differ from the formula by 1 ulp); test: within 1 ulp of the formula |
| Laplace decoding/encoding | Appendix A `laplace.c` | re-implemented from CELT_SPEC §2.1 |
| coarse/fine/final energy procedures | RFC §4.3.2 prose + Appendix A `quant_bands.c` | re-implemented from CELT_SPEC §2–3 |
| TF decoding, spread, boosts, trim | RFC §4.3.1, §4.3.3, §4.3.4.5 prose (corrected by Appendix A `celt.c`) | re-implemented from CELT_SPEC §4–5 |
| bit allocation (reservations, search, interpolation, skipping, intensity/dual, fine/shape split) | RFC §4.3.3 prose outline + Appendix A `rate.c` | re-implemented from CELT_SPEC §6–7 |
| band shapes (split recursion, theta, PVQ, spreading, TF Haar/Hadamard, folding, stereo merge, collapse masks) | RFC §4.3.4 prose + Appendix A `bands.c`, `vq.c`; RFC 8251 §9, §10 | re-implemented from CELT_SPEC §8 |
| anti-collapse | RFC §4.3.5 prose + Appendix A `bands.c` | re-implemented from CELT_SPEC §9 |
| denormalisation, IMDCT conventions, post-filter, de-emphasis, decoder state | RFC §4.3.6–4.3.7 prose + Appendix A `celt.c`, `mdct.c`; RFC 8251 §8 | re-implemented from CELT_SPEC §10–11 |
| CELT packet loss concealment (non-normative) | RFC §4.4 prose | new design (pitch-periodic continuation through the MDCT synthesis, noise after repeated losses) |
| PVQ index coding | RFC §4.3.4.2, §5.3.8.2 prose | unchanged (`src/celt/cwrs.rs` was written from the prose) |
| SILK-to-48 kHz resampler filters | none: least-squares fitted to vectors 02–04 | **removed**, with the fitting tool; replaced by the design below |

Where the RFC prose itself specifies an item it was used directly. Where it
defers to the reference code, the item rests on the normative Appendix A
read as specification through `docs/CELT_SPEC.md`.

## The resampler

RFC 6716 §4.2.9 makes the resampler non-normative and only its delay
(Table 54) normative. `src/resample.rs` designs every converter at
construction: a delay-constrained weighted least-squares FIR whose target is
a pure delay of exactly the Table 54 amount with a raised-cosine roll-off at
the lower Nyquist frequency, solved from its normal equations (the
derivation is in the module documentation). Equal rates are a whole-sample
delay line. No coefficient was fitted to any decoder's output; the handful
of design constants (band fractions and weights) were chosen on a small
grid against the RFC §6 conformance metric and passband flatness. Measured
responses are in `docs/VALIDATION.md`.

## Tables

- **SILK** (Tables 4–53 of RFC 6716): extracted from the RFC's plain text by
  a script (column parsing of each table, PDFs converted to inverse CDFs)
  into `src/silk/tables.rs`. Unit tests check that every PDF sums to 256,
  that the codebooks rise and that their weights fall in the RFC's stated
  range.
- **CELT**: band edges (Table 55), static allocation (Table 57), trim,
  spread, tapset and TF PDFs and adjustments (Tables 56, 58–63), post-filter
  taps and de-emphasis (§4.3.7) are in the RFC text; the remaining tables
  are listed above.

## Validation data

The official test vectors, `opus_testvectors-rfc8251.tar.gz` from
opus-codec.org (SHA-1s listed in RFC 8251 §11), are downloaded at test time
(CI does this) and used only as data: as the reference for the final range
and for the conformance metric. No filter or table is fitted to them.
