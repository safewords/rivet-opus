# CELT layer: normative specification for re-implementation

This document specifies the CELT layer of Opus (decoder, plus the parts an
encoder must mirror) precisely enough to reproduce the reference decoder's
range-coder state on every packet and its output to floating-point precision.
It was written from:

* RFC 6716 prose (cited as "RFC §x.y (prose)"),
* the normative reference implementation of RFC 6716 Appendix A, used only as
  a specification (cited as "App. A, `file`, `function`"); RFC 6716 §6 makes
  this code take precedence over the prose,
* RFC 8251, which patches that code (cited as "RFC 8251 §n").

No code is reproduced; procedures are described in prose and numbered steps.
Constant tables are data: they are produced exactly by
`tools/appendix_a_tables.py` (run it on `rfc6716.txt`; output names are given
as `TABLE_NAME` below), or derived as explained.

Where the RFC prose and the code disagree, **the code wins**; every such case is
flagged with **[PROSE ≠ CODE]** and collected in §14. The range-decoder
primitives used are listed in the appendix at the end.

Conventions used throughout:

* "BITRES" = 3: allocation quantities are in 1/8 bit ("eighth bits").
* `a >> k` on a signed integer is an **arithmetic (flooring) shift**; `a / b` on
  integers is C division, **truncating toward zero**. Where the distinction
  matters it is said again.
* `ilog(x)` = number of significant bits of an unsigned x (0 for 0, 1 for 1,
  2 for 2..3, ...), i.e. floor(log2 x)+1 (RFC §1.1.10).
* Bands, LM, M: the mode is the 48 kHz, 960-sample-maximum mode. A frame of
  `120·2^LM` samples at 48 kHz has LM ∈ {0,1,2,3} (2.5/5/10/20 ms), M = 2^LM.
  Band i covers bins `M·eBands[i] .. M·eBands[i+1]-1` of each channel, where
  `eBands = EBAND5MS = {0,1,2,3,4,5,6,7,8,10,12,14,16,20,24,28,34,40,48,60,78,100}`
  (RFC §4.3 Table 55, App. A `modes.c` `eband5ms`). nbEBands = 21; the number
  of "effective" bands effEBands = 21 (all band edges ≤ 120). "Width" of band
  i means `eBands[i+1]-eBands[i]` (the LM=0 width); the band size at the
  current LM is `width<<LM`.
* "Float build": the reference decoder whose output matches the test vectors
  is the floating-point build. All quantities that influence the range coder
  or the allocation are integers and must be bit-exact; quantities that only
  affect output samples are 32-bit floats and only need float precision.
  Each section says which is which.
* C = number of coded channels in the stream (1 or 2, from the TOC);
  CC = number of output channels of the decoder (1 or 2).

---

## 1. Frame-level decoding

Source: App. A, `celt.c`, `celt_decode_with_ec`; RFC §4.3 Table 56 (prose) for
the order and PDFs; the gating conditions are in the code only.

### 1.1 Inputs and parameters

* `len`: the number of bytes of the CELT frame. In CELT-only mode, the whole
  Opus frame. In Hybrid mode CELT continues decoding with the **same** range
  decoder that SILK used (no re-initialisation); `len` is the Opus frame length
  minus any redundancy bytes, and the decoder's buffer size (`storage`) is
  reduced by the redundancy bytes too (App. A `src/opus_decoder.c`,
  `opus_decode_frame`; RFC §4.5.1). For CELT-only frames a fresh range decoder
  is initialised on the frame.
* `start` (first coded band): 0 in CELT-only mode, 17 in Hybrid mode (RFC §4.3
  prose: "In Hybrid mode, the first 17 bands (up to 8 kHz) are not coded";
  App. A `opus_decoder.c`). Redundant 5 ms CELT frames inside SILK/Hybrid
  packets are decoded with start = 0.
* `end` (one past the last coded band), from the TOC bandwidth (App. A
  `opus_decoder.c`): narrowband 13, mediumband 17, wideband 17,
  super-wideband 19, fullband 21. `effEnd = min(end, effEBands) = end`.
* LM from the frame duration: 2.5 ms → 0, 5 → 1, 10 → 2, 20 → 3. N = 120·M
  is the number of MDCT bins (and output samples at 48 kHz) per channel.
* `total_bits = 8·len` (whole bits).
* `downsample` = 48000 / output rate ∈ {1,2,3,4,6} (§10.6).

If the frame is missing or `len ≤ 1`, the frame is concealed (§11.4) instead of
decoded. (`len` > 1275 is invalid.)

Before decoding, every coefficient of the normalised spectrum X (C channels ×
N bins) below `M·eBands[start]` and at or above `M·eBands[effEnd]` is set to 0.

If C = 1, first set, for every band i (all 21): `oldBandE[0][i] =
max(oldBandE[0][i], oldBandE[1][i])` (the decoder always keeps two channels of
energy state; see §11).

### 1.2 Symbol order and gating

`tell` below always means `ec_tell()` (whole bits, RFC §4.1.6.1) at that
point, and `tell_frac` means `ec_tell_frac()` (eighth bits, RFC §4.1.6.2).
Symbols marked "raw" are raw bits read from the end of the frame with
`ec_dec_bits` (RFC §4.1.4); all others are range-coded.

1. **silence**. If `tell ≥ total_bits`: silence = 1 (nothing read). Else if
   `tell == 1` exactly: silence = `ec_dec_bit_logp(15)`. Otherwise silence = 0.
   (So silence is never read in Hybrid mode, where SILK has used bits.)
   If silence = 1, the decoder's whole-bit counter is increased so that
   `ec_tell()` now returns exactly `8·len` ("pretend all bits were read"): add
   `8·len − ec_tell()` to `nbits_total`. Everything after that is decoded by
   the normal procedure; every gated symbol is then skipped by its gate and the
   ungated coarse energies take their "no bits" value (§2.3).
2. **post-filter** (only if `start == 0` and `tell + 16 ≤ total_bits`):
   one flag `ec_dec_bit_logp(1)`. If set:
   * octave = `ec_dec_uint(6)` (value 0..5);
   * period = `(16 << octave) + ec_dec_bits(4 + octave) − 1` (raw), range
     15..1022;
   * qg = `ec_dec_bits(3)` (raw); gain = 0.09375·(qg+1) (= 3(qg+1)/32);
   * tapset: if `ec_tell() + 2 ≤ total_bits`, `ec_dec_icdf(TAPSET_ICDF, 2)`
     (PDF {2,1,1}/4), else 0.
   If not set (or not read): period = 0, gain = 0, tapset = 0.
   **[PROSE ≠ CODE]** RFC §4.3.7.1 says the octave is "between 0 and 6"; it is
   `ec_dec_uint(6)`, i.e. 0..5 (consistent with the stated 1022 maximum).
3. **transient**: if `LM > 0` and `tell + 3 ≤ total_bits`,
   `ec_dec_bit_logp(3)`, else 0. shortBlocks = M if transient else 0 (in
   §8, B = M for transient frames, else 1).
4. **intra**: if `tell + 3 ≤ total_bits`, `ec_dec_bit_logp(3)`, else 0.
5. **coarse energy** for bands start..end−1 (§2).
6. **tf_change / tf_select** (§4).
7. **spread**: if `tell + 4 ≤ total_bits`, `ec_dec_icdf(SPREAD_ICDF, 5)`
   (PDF {7,2,21,2}/32), else 2 (SPREAD_NORMAL). Values: 0 none, 1 light,
   2 normal, 3 aggressive.
8. **caps** (§6.1), then **dynalloc boosts** (§5.2).
9. **allocation trim** (§5.3).
10. **anti-collapse reservation and allocation** (§6), which decodes the
    **skip** flags, **intensity** and **dual stereo**.
11. **fine energy** (§3.1).
12. **band shapes** (§8), passing the shape budget
    `8·8·len − anti_collapse_rsv` (eighth bits).
13. **anti-collapse flag**: if `anti_collapse_rsv > 0`, read 1 raw bit
    (`ec_dec_bits(1)`), else 0. **[PROSE ≠ CODE]** Table 56 lists PDF {1,1}/2;
    it is a raw bit.
14. **final fine energy bits** (§3.2) with `bits_left = 8·len − ec_tell()`.

Then synthesis: anti-collapse (§9) if the flag is set, energy to amplitude
(§10.1), silence override, denormalisation, inverse MDCT, post-filter,
de-emphasis (§10), and the state updates of §11.

**Silence override** (after the energy finalisation and anti-collapse): if
silence, then for all `i < C·21` (i.e. only the C coded channels): linear band
amplitude = 0 and `oldBandE = −28`.

### 1.3 After decoding (range decoder)

The decoder's RNG seed for the next frame is set to the range decoder's final
`rng` value (§11.2). The Opus layer's "final range" is the range decoder's
`rng` (XORed with the redundant frame's final `rng` if any); that is what the
test vectors check per packet.

---

## 2. Coarse energy

### 2.1 Laplace decoder

Source: App. A, `laplace.c`, `ec_laplace_decode` (RFC §4.3.2.1 prose names it
but gives no procedure). All arithmetic is unsigned 32-bit integer and exact.

Parameters: `fs` = probability of 0 in 1/32768 units, `decay` (Q14).
Constants: minimum probability MINP = 1, NMIN = 16.

1. `fm = ec_decode_bin(15)`; `fl = 0`; `val = 0`.
2. If `fm ≥ fs`:
   1. `val = 1`; `fl = fs`;
      `fs = (((32768 − 2·NMIN·MINP − fs) · (16384 − decay)) >> 15) + MINP`
      (with the old fs; i.e. `((32736 − fs)(16384 − decay)) >> 15`, plus 1).
   2. While `fs > MINP` and `fm ≥ fl + 2·fs`:
      `fs = 2·fs`; `fl = fl + fs`; `fs = (((fs − 2·MINP) · decay) >> 15) + MINP`;
      `val = val + 1`.
   3. If `fs ≤ MINP` (the tail of equiprobable values): `di = (fm − fl) >> 1`;
      `val = val + di`; `fl = fl + 2·di`.
   4. If `fm < fl + fs`, the value is negative: `val = −val`. Otherwise
      `fl = fl + fs` (positive value).
3. `ec_dec_update(fl, min(fl + fs, 32768), 32768)`; return val.

Each magnitude ≥ 1 thus occupies two adjacent intervals of width fs
(negative first, then positive).

### 2.2 Probability model

`E_PROB_MODEL[LM][intra][2·k]` and `[2·k+1]` (App. A `quant_bands.c`
`e_prob_model`, 4×2×42 bytes, Q8) are the zero probability and decay for band
index k = min(i, 20). The Laplace parameters are `fs = model[2k] << 7` and
`decay = model[2k+1] << 6`.

### 2.3 Decoding procedure

Source: App. A, `quant_bands.c`, `unquant_coarse_energy`; RFC §4.3.2.1 (prose)
for the predictor.

Coefficients (floats in the float build; App. A `quant_bands.c`):
* inter frames: α = `PRED_COEF[LM]` = {29440, 26112, 21248, 16384}/32768,
  β = `BETA_COEF[LM]` = {30147, 22282, 12124, 6554}/32768;
* intra frames: α = 0, β = `BETA_INTRA` = 4915/32768.

`budget = 8·storage` (the range decoder's buffer size in bits, which equals
`8·len`). For each channel keep `prev[c]`, initialised to 0.

For i = start..end−1, and for c = 0..C−1 (channel loop innermost):

1. `t = ec_tell()`.
2. Decode the integer residual qi:
   * if `budget − t ≥ 15`: qi = Laplace decode (§2.1) with the parameters of
     §2.2;
   * else if `budget − t ≥ 2`: s = `ec_dec_icdf(SMALL_ENERGY_ICDF, 2)` (PDF
     {2,1,1}/4), qi = 0 for s=0, −1 for s=1, +1 for s=2 (i.e.
     `(s>>1) XOR −(s&1)`);
   * else if `budget − t ≥ 1`: qi = −`ec_dec_bit_logp(1)`;
   * else qi = −1 (nothing read).
3. Prediction (float): `old = max(−9, oldBandE[c][i])`;
   `oldBandE[c][i] = α·old + prev[c] + qi`;
   `prev[c] = prev[c] + qi − β·qi`.

Energies are in log2 units (1.0 = 6.02 dB) relative to the band mean
(§10.1). The time-domain prediction uses the previous frame's final energies
(including fine and final-fine refinements). The float build has no further
clamping. (The fixed-point build clamps the predicted value at −28; this only
matters for fixed point and is not part of the float behaviour.) Energies only
influence output samples, never the bitstream parsing, so the float
arithmetic here need not be bit-exact.

---

## 3. Fine energy

### 3.1 Fine energy bits

Source: RFC §4.3.2.2 (prose: correction `(f+1/2)/2^B − 1/2`); App. A
`quant_bands.c`, `unquant_fine_energy`.

For i = start..end−1 with `fine_quant[i] > 0`, for c = 0..C−1:
q2 = `ec_dec_bits(fine_quant[i])` (raw), and
`oldBandE[c][i] += (q2 + 0.5)/2^fine_quant[i] − 0.5`.
(In the float build this is computed as `(q2+0.5)·2^(14−fq)·(1/16384) − 0.5`,
which is exact.)

### 3.2 Final fine bits

Source: RFC §4.3.2.2 (prose); App. A `quant_bands.c`,
`unquant_energy_finalise`.

`bits_left = 8·len − ec_tell()` (after the anti-collapse bit). For priority
p = 0, then p = 1: for i = start..end−1, **checking `bits_left ≥ C` before each
band** (stop the pass for this priority once it fails): skip the band if
`fine_quant[i] ≥ 8` (MAX_FINE_BITS) or `fine_priority[i] ≠ p`; otherwise for
c = 0..C−1: q2 = `ec_dec_bits(1)` (raw);
`oldBandE[c][i] += (q2 − 0.5)/2^(fine_quant[i]+1)`; `bits_left −= 1`.

The RFC prose ("starting from band 0") omits the per-band `bits_left ≥ C`
test; for stereo both channels of a band always get their bit together.

---

## 4. Time-frequency resolution (tf_change, tf_select)

Source: RFC §4.3.1 and §4.3.4.5 (prose, PDFs and Tables 60–63); App. A
`celt.c`, `tf_decode` (budget rules). Integer, exact.

1. `budget = 8·storage`, `t = ec_tell()`. logp = 2 if transient else 4.
2. `tf_select_rsv = (LM > 0 and t + logp + 1 ≤ budget)`; if so,
   `budget −= 1`.
3. `curr = 0`, `tf_changed = 0`. For i = start..end−1:
   * if `t + logp ≤ budget`: `curr = curr XOR ec_dec_bit_logp(logp)`;
     `t = ec_tell()`; `tf_changed = tf_changed OR curr`;
   * `tf_res[i] = curr`;
   * logp = 4 if transient else 5 (from the second band on, whether or not
     the first was read).
4. tf_select = 0. If `tf_select_rsv` and
   `TF_SELECT_TABLE[LM][4·transient + tf_changed] ≠
   TF_SELECT_TABLE[LM][4·transient + 2 + tf_changed]`, then tf_select =
   `ec_dec_bit_logp(1)`.
5. For i = start..end−1: `tf_res[i] = TF_SELECT_TABLE[LM][4·transient +
   2·tf_select + tf_res[i]]`.

`TF_SELECT_TABLE` (App. A `celt.c` `tf_select_table`, rows LM=0..3, columns
[non-transient sel0: flag0, flag1, non-transient sel1: 0,1, transient sel0:
0,1, transient sel1: 0,1]):
LM0 {0,−1,0,−1, 0,−1,0,−1}; LM1 {0,−1,0,−2, 1,0,1,−1};
LM2 {0,−2,0,−3, 2,0,1,−1}; LM3 {0,−2,0,−3, 3,0,1,−1}.
This equals RFC Tables 60–63 exactly.

The resulting `tf_res[i]` is the per-band tf_change used in §8.4: positive =
more frequency resolution (only meaningful for transients), negative = more
time resolution.

Encoder (App. A `tf_encode`): identical gating; when a flag cannot be coded the
encoder must replace its decision by the running value `curr`, and it codes
tf_select only under the same condition (otherwise it is 0).

---

## 5. Spread, dynalloc boosts, trim

### 5.1 Spread

See §1.2 step 7. `SPREAD_ICDF = {25,23,2,0}` (ftb 5).

### 5.2 Band boosts (dynalloc)

Source: RFC §4.3.3 (prose, detailed); App. A `celt.c`. Integer, exact.

1. `dynalloc_logp = 6`; `total = 8·8·len` (eighth bits); `t = tell_frac`.
2. For i = start..end−1:
   1. `width = C · (width_i << LM)` (bins of all coded channels);
      `quanta = min(8·width, max(48, width))`.
   2. `loop_logp = dynalloc_logp`; `boost = 0`.
   3. While `t + 8·loop_logp < total` and `boost < cap[i]`:
      flag = `ec_dec_bit_logp(loop_logp)`; `t = tell_frac`; if flag = 0 stop;
      else `boost += quanta`, `total −= quanta`, `loop_logp = 1`.
   4. `offsets[i] = boost`; if `boost > 0`,
      `dynalloc_logp = max(2, dynalloc_logp − 1)`.

`total` keeps its reduced value for the trim gate.

**[PROSE ≠ CODE]** RFC §4.3.3 states the loop condition as "tell plus the
cost is less than total_bits plus total_boost", which (since the prose also
subtracts each quanta from total_bits) compares against the *original* frame
size. The code compares against the *reduced* total (`total` above). Code
wins. The prose's "N" in the quanta formula is the bin count of all coded
channels (C·(width<<LM)).

### 5.3 Allocation trim

Source: RFC §4.3.3 (prose, Table 58); App. A `celt.c`.
`alloc_trim = ec_dec_icdf(TRIM_ICDF, 7)` if `t + 48 ≤ total` (t = the latest
tell_frac, `total` = the reduced total of §5.2), else 5.
`TRIM_ICDF = {126,124,119,109,87,41,19,9,4,2,0}` = PDF
{2,2,5,10,22,46,22,10,5,2,2}/128. This matches the prose exactly.

---

## 6. Bit allocation

Source: RFC §4.3.3 (prose, a substantial outline); App. A `celt.c` (`init_caps`,
the reservation), `rate.c` (`compute_allocation`, `interp_bits2pulses`). The
whole of §6 is integer arithmetic and must be bit-exact.

### 6.1 Caps

Source: RFC §4.3.3 (prose) and App. A `celt.c` `init_caps`:
`cap[i] = ((CACHE_CAPS50[2·LM + C − 1][i] + 64) · C · (width_i << LM)) >> 2`
for all 21 bands. (`CACHE_CAPS50` is 8 rows × 21 bands; row index
`2·LM + C − 1`, as in the prose's `nbBands·(2·LM+stereo)`.)

`CACHE_CAPS50` is derivable from the pulse cache; see §7.3 (derivation verified
to reproduce the table exactly).

### 6.2 Reservations

Source: App. A `celt.c` and `rate.c` `compute_allocation`; RFC §4.3.3 (prose).

1. `bits = 8·8·len − tell_frac − 1` (after the trim).
2. `anti_collapse_rsv = 8` if transient and `LM ≥ 2` and `bits ≥ 8·(LM+2)`,
   else 0. `bits −= anti_collapse_rsv`.
3. Inside the allocation: `total = max(bits, 0)`.
4. `skip_rsv = 8` if `total ≥ 8`, else 0; `total −= skip_rsv`.
   **[PROSE ≠ CODE]** prose says "greater than 8"; the code uses ≥.
5. If C = 2: `intensity_rsv = LOG2_FRAC_TABLE[end − start]`; if
   `intensity_rsv > total`, `intensity_rsv = 0`; otherwise
   `total −= intensity_rsv` and `dual_stereo_rsv = 8` if `total ≥ 8` else 0,
   `total −= dual_stereo_rsv`. If C = 1 both are 0.
   **[PROSE ≠ CODE]** prose says dual is reserved if total "is still greater
   than 8"; the code uses ≥. The prose's "ilog2(end−start) bits" is really
   `LOG2_FRAC_TABLE[end−start]` eighth bits (a conservative log2 of the number
   of intensity values).

`LOG2_FRAC_TABLE` (24 entries, App. A `rate.c`) =
{0, 8,13, 16,19,21,23, 24,26,27,28,29,30,31,32, 32,33,34,34,35,36,36,37,37}.
Derivation: entry i = ⌈8·log2(i+1)⌉ = `log2_frac(i+1, 3)` (§7.1); both give the
table exactly.

### 6.3 Per-band thresholds and trim offsets

For j = start..end−1 (N_j = width_j):
* `thresh[j] = max(8·C, (3·(N_j << LM) << 3) >> 4)` (one bit per channel, or
  3/16 bit per bin — the prose's "24·N/16").
* `trim_offset[j] = (C · N_j · (alloc_trim − 5 − LM) · (end − j − 1) ·
  2^(LM+3)) >> 6`, an **arithmetic shift of a possibly negative product
  (flooring)**. **[PROSE ≠ CODE]** the prose says "divide by 64"; truncating
  division gives different results for negative values. Then, if
  `N_j << LM == 1`, `trim_offset[j] −= 8·C`.

### 6.4 Search over the static table rows

`BAND_ALLOCATION` (11 rows × 21 bands, units 1/32 bit per bin; App. A
`modes.c` `band_allocation`; equal to RFC Table 57 transposed).

Define for row q and band j:
`raw(q, j) = (C · N_j · BAND_ALLOCATION[q][j] << LM) >> 2`, and
`adj(q, j) = max(0, raw + trim_offset[j])` if `raw > 0`, else `raw` (= 0).

Bisection: `lo = 1`, `hi = 10`. Repeat while `lo ≤ hi`:
1. `mid = (lo + hi) >> 1`; `psum = 0`; `done = false`.
2. For j = end−1 down to start: `v = adj(mid, j) + offsets[j]`;
   if `v ≥ thresh[j]` or done: `done = true`, `psum += min(v, cap[j])`;
   else if `v ≥ 8·C`: `psum += 8·C`.
3. If `psum > total`: `hi = mid − 1`, else `lo = mid + 1`.

Then `hi = lo`, `lo = lo − 1` (so 0 ≤ lo ≤ 10, hi = lo + 1 ≤ 11).

For j = start..end−1:
* `bits1 = adj(lo, j)`; `bits2 = adj(hi, j)` if hi ≤ 10, else
  `max(0, cap[j] + trim_offset[j])` (the cap is treated as row 11, and the
  "if > 0, add trim and clamp" rule is applied to it too);
* if `lo > 0`: `bits1 += offsets[j]`; always `bits2 += offsets[j]`;
* if `offsets[j] > 0`: `skip_start = j` (initially skip_start = start; the
  last boosted band wins);
* `bits2 = max(0, bits2 − bits1)`.

### 6.5 Interpolation (1/64 steps)

`alloc_floor = 8·C`. Bisection over 6 steps: `lo = 0`, `hi = 64`; six times:
`mid = (lo + hi) >> 1`; compute psum as in §6.4 step 2 with
`v = bits1[j] + ((mid · bits2[j]) >> 6)` (thresholds and cap as there, floor
contribution `alloc_floor`); if `psum > total`, `hi = mid`, else `lo = mid`.

Then, with this `lo`, for j = end−1 down to start (with `done = false`
initially): `v = bits1[j] + ((lo · bits2[j]) >> 6)`; if `v < thresh[j]` and
not done: `v = alloc_floor` if `v ≥ alloc_floor`, else 0; otherwise
`done = true`. `bits[j] = min(v, cap[j])`; `psum += bits[j]` (psum restarted
at 0).

### 6.6 Skipping (decoder), and the encoder's choice point

Loop with `codedBands` starting at `end`, decrementing each iteration;
`j = codedBands − 1`:

1. If `j ≤ skip_start`: `total += skip_rsv`; stop. (The first band, and the
   last boosted band and those below it, are never skipped.)
2. `left = total − psum`; `D = eBands[codedBands] − eBands[start]`;
   `percoeff = left / D`; `left −= D·percoeff`;
   `rem = max(left − (eBands[j] − eBands[start]), 0)`;
   `band_width = eBands[codedBands] − eBands[j]`;
   `band_bits = bits[j] + percoeff·band_width + rem`.
   (Widths here are LM=0 widths, not scaled by M.)
3. If `band_bits ≥ max(thresh[j], alloc_floor + 8)`:
   * decoder: if `ec_dec_bit_logp(1)` = 1, stop (keep codedBands);
   * encoder: it codes 1 (and stops) when it decides to keep the band. The
     reference encoder keeps the band if
     `band_bits > ((j < prev ? 7 : 9) · band_width << LM << 3) >> 4`, where
     `prev` is the previous frame's codedBands — this is the only
     non-normative decision in the allocation;
   * after a 0: `psum += 8`; `band_bits −= 8`.
   Otherwise the band is skipped without signalling.
4. `psum −= bits[j] + intensity_rsv`; if `intensity_rsv > 0`,
   `intensity_rsv = LOG2_FRAC_TABLE[j − start]`; `psum += intensity_rsv`.
5. If `band_bits ≥ alloc_floor`: `psum += alloc_floor`, `bits[j] =
   alloc_floor` (one fine-energy bit per channel); else `bits[j] = 0`.
6. Continue with `codedBands − 1`.

At the end codedBands > start.

### 6.7 Intensity and dual stereo

* If `intensity_rsv > 0`: `intensity = start + ec_dec_uint(codedBands + 1 −
  start)`; else `intensity = 0`. (Encoder: `intensity = min(intensity,
  codedBands)` before coding.)
* If `intensity ≤ start`: `total += dual_stereo_rsv`, `dual_stereo_rsv = 0`.
* If `dual_stereo_rsv > 0`: `dual_stereo = ec_dec_bit_logp(1)`, else 0.

### 6.8 Distributing the remainder

`left = total − psum`; `D = eBands[codedBands] − eBands[start]`;
`percoeff = left / D`; `left −= D·percoeff`. For j = start..codedBands−1:
`bits[j] += percoeff · width_j`. Then for j = start..codedBands−1:
`t = min(left, width_j)`; `bits[j] += t`; `left −= t`.

### 6.9 Fine energy / shape split

`balance = 0`; `logM = 8·LM`. For j = start..codedBands−1:

1. `N = width_j << LM`; `bits[j] += balance`.
2. If `N > 1`:
   1. `excess = max(bits[j] − cap[j], 0)`; `bits[j] −= excess`.
   2. `den = C·N + 1` if (C = 2 and N > 2 and not dual_stereo and
      j < intensity), else `C·N`.
   3. `NClogN = den · (logN[j] + logM)`, where `logN[j] = LOGN400[j]` =
      `log2_frac(width_j, 3)` = ⌈8·log2(width_j)⌉ =
      {0,0,0,0,0,0,0,0,8,8,8,8,16,16,16,21,21,24,29,34,36}.
   4. `offset = (NClogN >> 1) − 21·den` (FINE_OFFSET = 21); if `N == 2`,
      `offset += 2·den` (`den << 3 >> 2`).
   5. If `bits[j] + offset < 16·den`: `offset += NClogN >> 2`; else if
      `bits[j] + offset < 24·den`: `offset += NClogN >> 3`.
   6. `ebits[j] = max(0, (bits[j] + offset + 4·den) / (8·den))`.
   7. If `C·ebits[j] > bits[j] >> 3`: `ebits[j] = (bits[j] >> (C−1)) >> 3`.
   8. `ebits[j] = min(ebits[j], 8)`.
   9. `fine_priority[j] = (ebits[j] · 8·den ≥ bits[j] + offset)`.
   10. `bits[j] −= 8·C·ebits[j]`.
3. Else (N = 1): `excess = max(0, bits[j] − 8·C)`; `bits[j] −= excess`;
   `ebits[j] = 0`; `fine_priority[j] = 1`.
4. If `excess > 0`: `extra_fine = min(excess >> (C−1+3), 8 − ebits[j])`;
   `ebits[j] += extra_fine`; `extra_bits = 8·C·extra_fine`;
   `fine_priority[j] = (extra_bits ≥ excess − balance)` (balance = the value
   carried *into* this band); `excess −= extra_bits`.
5. `balance = excess`.

The final `balance` is returned (the "balance" passed to §8). For the skipped
bands j = codedBands..end−1: `ebits[j] = (bits[j] >> (C−1)) >> 3`,
`bits[j] = 0`, `fine_priority[j] = (ebits[j] < 1)`.

Outputs: `pulses[j] = bits[j]` (shape budget in eighth bits),
`fine_quant[j] = ebits[j]`, `fine_priority[j]`, `codedBands`, `balance`,
`intensity`, `dual_stereo`.

---

## 7. Pulse cache, bits ↔ pulses

Source: RFC §4.3.4.1 (prose: nearest K, ties round down, balance); RFC
§4.3.4.2 (prose: V(N,K)); App. A `rate.h` (`get_pulses`, `bits2pulses`,
`pulses2bits`), `rate.c` `compute_pulse_cache` (the generator, compiled only
for custom modes), `cwrs.c` (`log2_frac`, `get_required_bits`),
`static_modes_float.h` (`cache_index50`, `cache_bits50`, `cache_caps50`).
Integer, exact.

### 7.1 Conservative log2: `log2_frac(val, frac)`

Returns an integer fixed-point log2 (frac fractional bits) of an unsigned
32-bit `val` > 0 that is never smaller than the true value (App. A `cwrs.c`
`log2_frac`):

1. `l = ilog(val)`.
2. If val is a power of two, return `(l − 1) << frac`.
3. Otherwise normalise val to 16 fractional bits, **rounding up**: if
   `l > 16`, `val = ⌈val / 2^(l−16)⌉` (computed without overflow as
   `val >> (l−16)`, plus 1 if any of the discarded low bits is set); else
   `val = val << (16 − l)`. (Now 2^15 < val ≤ 2^16.)
4. `l = (l − 1) << frac`.
5. For f = frac, frac−1, …, 0 (frac+1 iterations): `b = val >> 16` (0 or 1);
   `l += b << f`; `val = (val + b) >> b`; `val = (val·val + 0x7FFF) >> 15`
   (32-bit unsigned).
6. Return `l + 1` if `val > 0x8000`, else `l`.

This is NOT always equal to ⌈2^frac·log2(val)⌉: in the cache below there is
exactly one difference (N = 11, K = 9: 8·log2 V(11,9) = 176.99997, the ceiling
is 177 but log2_frac gives 178). Use this procedure.

### 7.2 The pulse cache

`V(N,K)` = number of N-dimensional integer vectors with sum of absolute values
K: V(N,0) = 1, V(0,K>0) = 0, V(N,K) = V(N−1,K) + V(N,K−1) + V(N−1,K−1)
(RFC §4.3.4.2 prose).

Pseudo-pulse mapping: `get_pulses(q) = q` for q < 8, else
`(8 + (q AND 7)) << ((q >> 3) − 1)` (q = 8..40 → K = 8,9,…,15, 16,18,…,30,
32,36,…,60, 64,72,…,120, 128). MAX_PSEUDO = 40.

Layout (`CACHE_INDEX50`, 5 rows × 21; `CACHE_BITS50`, 392 bytes): row r ∈
0..4 corresponds to "LM = r − 1" and to band size `N = (width_j << r) >> 1`
(r = 0 gives N = width/2, which is 0 for width-1 bands → index −1). Each
distinct N (scanning r = 0..4 outer, j = 0..20 inner) gets an entry at the
next free offset `o` in CACHE_BITS50; repeated sizes reuse the first entry's
offset. Entry content:
* `bits[o] = Kmax`, the largest q ≤ 40 such that V(N, get_pulses(q)) < 2^32
  (equivalent to App. A `fits_in32`; verified);
* `bits[o + q] = log2_frac(V(N, get_pulses(q)), 3) − 1` for q = 1..Kmax
  (for N = 1, V = 2, so every entry is 7);
* entry length Kmax + 1.

This construction reproduces `cache_index50` and `cache_bits50` exactly
(verified). Entries in order, as (N, Kmax): (1,40) (2,40) (3,40) (4,40) (6,35)
(9,21) (11,17) (8,25) (12,16) (18,11) (22,9) (16,12) (24,9) (36,7) (44,6)
(32,7) (48,6) (72,5) (88,5) (64,5) (96,5) (144,4) (176,4); total 392 bytes.

For band i at "LM" l (l may be −1 after splits) the entry is
`CACHE_BITS50[CACHE_INDEX50[l+1][i] ..]`, called `cache` below
(`cache[0]` = Kmax).

### 7.3 Derivation of the caps (`CACHE_CAPS50`)

Source: App. A `rate.c` `compute_pulse_cache` (second half). For LM i = 0..3,
C = 1..2, band j (row `2·i + C − 1`), with width N0 = width_j and logN as in
§6.9 (integers; all divisions here have non-negative operands):

1. If `N0 << i == 1`: `max_bits = 8·C·(1 + 8)`.
2. Else:
   1. `LM0 = 0`. If N0 > 2: `N0 = N0/2`, `LM0 = −1`. Else if N0 ≤ 1:
      `LM0 = min(i, 1)`, `N0 = N0 << LM0`.
   2. `p` = the cache entry of band j at row LM0+1;
      `max_bits = p[p[0]] + 1` (cost of the largest codebook).
   3. `N = N0`. For k = 0 .. i − LM0 − 1: `max_bits = 2·max_bits`;
      `offset = ((logN[j] + 8·(LM0 + k)) >> 1) − 4`;
      `num = 459·((2N − 1)·offset + max_bits)`; `den = ((2N − 1) << 9) − 459`;
      `qb = min((num + den/2) / den, 57)`; `max_bits += qb`; `N = 2N`.
   4. If C = 2: `max_bits = 2·max_bits`;
      `offset = ((logN[j] + 8·i) >> 1) − (16 if N == 2 else 4)`;
      `ndof = 2N − 1 − (1 if N == 2 else 0)`; `f = 512 if N == 2 else 487`;
      `num = f·(max_bits + ndof·offset)`; `den = (ndof << 9) − f`;
      `qb = min((num + den/2)/den, 64 if N == 2 else 61)`; `max_bits += qb`.
   5. `ndof = C·N + (1 if C == 2 and N > 2 else 0)`;
      `offset = ((logN[j] + 8·i) >> 1) − 21`, plus 2 if N == 2;
      `num = max_bits + ndof·offset`; `den = (ndof − 1) << 3`;
      `qb = min((num + den/2)/den, 8)`; `max_bits += 8·C·qb`.
3. `cap = 4·max_bits / (C·(width_j << i)) − 64` (fits a byte).

This reproduces `cache_caps50` exactly (verified).

### 7.4 bits2pulses and pulses2bits

Source: App. A `rate.h`. For band i at LM l and budget `b` (eighth bits),
`cache` as in §7.2:

* `bits2pulses`: `lo = 0`, `hi = cache[0]`, `b' = b − 1`. Six times:
  `mid = (lo + hi + 1) >> 1`; if `cache[mid] ≥ b'`, `hi = mid`, else
  `lo = mid`. Return lo if `b' − (−1 if lo == 0 else cache[lo]) ≤
  cache[hi] − b'`, else hi. (Nearest cost, ties to the smaller q — RFC
  §4.3.4.1.)
* `pulses2bits(q) = 0` if q = 0, else `cache[q] + 1`.

---

## 8. Band shapes

Source: App. A `bands.c` (`quant_all_bands`, `quant_band`, `compute_qn`,
helpers), `vq.c`, `cwrs.c`; RFC §4.3.4.x (prose outline); RFC 8251 §9 and §10.
Everything that determines what is read from the range decoder (budgets,
splits, qn, theta, pulses, the shared `remaining_bits`) is integer and must
be exact; the float operations (rotations, normalisation, Haar steps, stereo
merge, folding values) only affect output samples.

Notation: in the decoder, the normalised spectrum `X` (channel 0) and `Y`
(channel 1, absent for mono) are float arrays of N = 120·M bins. A per-frame
scratch buffer `norm` (and `norm2` for channel 1) of `M·eBands[21]` = 100·M
floats holds, for each decoded band, a scaled copy of its normalised output
used as the folding source of later bands. Its content before being written
in a frame is irrelevant (only already-written parts are ever read).

### 8.1 Per-band loop (`quant_all_bands`)

Inputs: `pulses[]`, `tf_res[]`, `spread`, `dual_stereo`, `intensity`,
`codedBands`, `balance` (from §6.9), `shape_total = 8·8·len −
anti_collapse_rsv`, and the 32-bit `seed` (decoder state, §11.2), updated in
place by the folding generator.

`B = M` if transient else 1. `lowband_offset = 0`, `update_lowband = true`.
For i = start..end−1 (`N = M·width_i`, band offset `o_i = M·eBands[i]`):

1. `tell = tell_frac`. If `i ≠ start`: `balance −= tell`.
2. `remaining_bits = shape_total − tell − 1`.
3. Budget: if `i ≤ codedBands − 1`:
   `curr_balance = balance / min(3, codedBands − i)` (**truncating**; balance
   may be negative), `b = max(0, min(16383, min(remaining_bits + 1,
   pulses[i] + curr_balance)))`; else `b = 0`.
4. Folding position (RFC 8251 §9 form): if (`o_i − N ≥ M·eBands[start]` or
   `i == start + 1`) and (`update_lowband` or `lowband_offset == 0`):
   `lowband_offset = i`.
5. RFC 8251 §9: if `i == start + 1`: let `n1 = M·width_start`,
   `n2 = M·width_{start+1}`, `s = M·eBands[start]`; copy the `n2 − n1`
   values `norm[s + 2·n1 − n2 .. s + n1)` to `norm[s + n1 .. s + n2)` (and the
   same in norm2 if C = 2). In CELT-only mode (start = 0) n1 = n2 and nothing
   is copied; in Hybrid (start = 17) this repeats the last 4M values of band
   17's folding data so band 18 (12M bins) can fold from it.
6. `tf_change = tf_res[i]`.
7. Folding source and collapse estimate: if `lowband_offset ≠ 0` and
   (`spread ≠ 3` or `B > 1` or `tf_change < 0`):
   * `eff = max(M·eBands[start], M·eBands[lowband_offset] − N)`; the folding
     source is `norm[eff .. eff + N)` (and `norm2[...]` for dual stereo);
   * `fold_start` = the largest index f < lowband_offset with
     `M·eBands[f] ≤ eff`;
   * `fold_end` = the first index f ≥ lowband_offset with `f ≥ i` or
     `M·eBands[f] ≥ eff + N` (RFC 8251 §9 added the `f ≥ i` bound);
   * `x_cm` = OR of `collapse_masks[f][0]`, `y_cm` = OR of
     `collapse_masks[f][C−1]`, over f = fold_start..max(fold_start,
     fold_end−1).

   Otherwise there is no folding source (noise is used) and
   `x_cm = y_cm = 2^B − 1`.
8. If `dual_stereo` and `i == intensity`: dual stereo is switched off for the
   rest of the frame, and `norm[k] = 0.5·(norm[k] + norm2[k])` for
   k = M·eBands[start] .. o_i − 1.
9. Decode:
   * dual stereo: `x_cm = band(X, mono, b >> 1, lowband = norm + eff (or
     none), lowband_out = norm + o_i, fill = x_cm)`, then
     `y_cm = band(Y, mono, b >> 1, lowband = norm2 + eff (or none),
     lowband_out = norm2 + o_i, fill = y_cm)`;
   * otherwise: `x_cm = band(X [and Y if C = 2], b, lowband = norm + eff (or
     none), lowband_out = norm + o_i, fill = x_cm OR y_cm)`, `y_cm = x_cm`.
   All with level 0, gain 1.0, the frame's B and LM, tf_change, spread,
   intensity.
10. `collapse_masks[i][0] = x_cm`, `collapse_masks[i][C−1] = y_cm` (low 8
    bits).
11. `balance += pulses[i] + tell`.
12. `update_lowband = (b > 8·N)`.

`remaining_bits` is a single counter shared by all nested calls of a band
and decremented as bits are consumed. Note that in joint (non-dual) stereo
only `norm` is written (from the mid), never `norm2`.

The encoder runs the same loop. It need not compute folding sources, collapse
masks or any resynthesis: they never influence the bitstream (the reference
encoder without resynthesis keeps `lowband_offset = 0`); but every budget,
split and theta resolution must match.

### 8.2 Band decoding: `band(...)`

Arguments: X, optional Y (stereo), N, b, spread, B, intensity, tf_change,
lowband (may be absent), `remaining_bits` (shared), LM, lowband_out (may be
absent), level, gain, fill. Returns the collapse mask `cm` (initially 0).
Local copies: `N0 = N`, `B0 = B`, `longBlocks = (B0 == 1)`, `N_B = N / B`,
stereo = (Y present), `split = stereo`, `inv = 0`, `mid = side = 0`.

#### 8.2.1 N = 1

For each of the 1 or 2 channels (X then Y): if `remaining_bits ≥ 8`: read one
raw bit `sign`, `remaining_bits −= 8`, `b −= 8`; else sign = 0. Set the
coefficient to −1.0 if sign else +1.0. If lowband_out is present,
`lowband_out[0] = X[0]`. Return 1. (No theta, no stereo merge, no TF.)

#### 8.2.2 Time-frequency changes on entry (mono only, level 0)

If not stereo and level = 0:

1. `recombine = max(tf_change, 0)`.
2. If a lowband is present and (`recombine > 0` or (`N_B` even and
   `tf_change < 0`) or `B0 > 1`), work on a private copy of it (the
   transforms below must not modify `norm`). Always copying is equivalent.
3. For k = 0..recombine−1: if lowband present, `haar(lowband, N >> k,
   1 << k)`; `fill = BIT_INTERLEAVE_TABLE[fill AND 15] OR
   (BIT_INTERLEAVE_TABLE[fill >> 4] << 2)`.
4. `B = B >> recombine`; `N_B = N_B << recombine`.
5. `time_divide = 0`. While `N_B` is even and `tf_change < 0`: if lowband
   present, `haar(lowband, N_B, B)`; `fill = fill OR (fill << B)`;
   `B = 2B`; `N_B = N_B / 2`; `time_divide += 1`; `tf_change += 1`.
6. `B0 = B`; `N_B0 = N_B`.
7. If `B0 > 1` and lowband present: `deinterleave(lowband, N_B >> recombine,
   B0 << recombine, hadamard = longBlocks)`.

`haar(v, n, stride)`: for each i < stride and j < n/2, with
`a = v[stride·2j + i]`, `c = v[stride·(2j+1) + i]`:
`v[stride·2j + i] = 0.70710678·a + 0.70710678·c`,
`v[stride·(2j+1) + i] = 0.70710678·a − 0.70710678·c` (float: each product
first, then the sum/difference).

`deinterleave(v, n0, stride, hadamard)` builds `t` of length n0·stride with
`t[ord(i)·n0 + j] = v[j·stride + i]` for i < stride, j < n0, where
`ord(i) = ORDERY_TABLE[stride − 2 + i]` if hadamard, else `ord(i) = i`; then
`v = t`. `interleave` is its inverse: `t[j·stride + i] = v[ord(i)·n0 + j]`.
(Coefficients interleaved by block become contiguous per block; for
non-transient frames the blocks come out in "sequency" order.)

`BIT_INTERLEAVE_TABLE[x]` (4-bit x): bit 0 = ((x AND 3) ≠ 0), bit 1 =
((x AND 12) ≠ 0), i.e. {0,1,1,1,2,3,3,3,2,3,3,3,2,3,3,3}.
`BIT_DEINTERLEAVE_TABLE[x]` duplicates each bit k of x into bits 2k and 2k+1:
{0x00,0x03,0x0C,0x0F,0x30,0x33,0x3C,0x3F,0xC0,0xC3,0xCC,0xCF,0xF0,0xF3,0xFC,0xFF}.
`ORDERY_TABLE` (strides 2, 4, 8, 16 concatenated):
{1,0, 3,0,2,1, 7,0,4,3,6,1,5,2, 15,0,8,7,12,3,11,4,14,1,9,6,13,2,10,5};
derivation: for stride S = 2^s, entry i = (S − 1) − G⁻¹(bitrev_s(i)), where
bitrev_s reverses the s-bit index and G⁻¹ is the inverse Gray code
(x XOR x>>1 XOR x>>2 XOR …). All three derivations reproduce the tables
exactly (verified).

#### 8.2.3 Split decision (mono)

`cache` = pulse-cache entry of band i at the current LM (§7.2). If not
stereo, `LM ≠ −1`, `b > cache[cache[0]] + 12` and `N > 2`:
`N = N/2`; `Y = X + N` (the second half becomes the "side"); `split = true`;
`LM = LM − 1`; if `B == 1`, `fill = (fill AND 1) OR (fill << 1)`;
`B = (B + 1) >> 1`. (A band can thus be split at most LM+1 times — RFC
§4.3.4.4.)

#### 8.2.4 Theta (split parameter)

When `split` (mono split or stereo):

1. `pulse_cap = logN[i] + 8·LM` (current LM);
   `offset = (pulse_cap >> 1) − (16 if stereo and N == 2 else 4)`
   (QTHETA_OFFSET_TWOPHASE = 16, QTHETA_OFFSET = 4).
2. `qn` (App. A `compute_qn`): `N2 = 2N − 1`, minus 1 if stereo and N == 2;
   `qb = min(b − pulse_cap − 32, (b + N2·offset) / N2)` (truncating);
   `qb = min(64, qb)`; if `qb < 4`, `qn = 1`; else
   `qn = EXP2_TABLE8[qb AND 7] >> (14 − (qb >> 3))`, then rounded up to even:
   `qn = ((qn + 1) >> 1) << 1`. `EXP2_TABLE8 =
   {16384,17866,19483,21247,23170,25267,27554,30048}` = ⌊16384·2^(k/8)⌋
   (exact; note it is the floor, not the rounding). qn ∈ {1, 2, 4, …, 256}.
3. If stereo and `i ≥ intensity`: `qn = 1`.
4. `t0 = tell_frac`.
5. If `qn ≠ 1`, decode `itheta` ∈ 0..qn:
   * stereo and N > 2 — **step PDF**: `p0 = 3`, `x0 = qn/2`,
     `ft = p0·(x0+1) + x0`. `fs = ec_decode(ft)`. If `fs < (x0+1)·p0`,
     `x = fs / p0`, else `x = x0 + 1 + (fs − (x0+1)·p0)`. Its interval is
     [p0·x, p0·(x+1)) if `x ≤ x0`, else [(x−1−x0) + (x0+1)·p0,
     (x−x0) + (x0+1)·p0). `ec_dec_update` with it; itheta = x.
     (Values ≤ qn/2 have weight 3, the others weight 1.)
   * else if `B0 > 1` or stereo — **uniform**: itheta = `ec_dec_uint(qn+1)`.
   * else — **triangular PDF** (weights 1,2,…,qn/2+1,…,2,1):
     `h = qn >> 1`, `ft = (h+1)²`, `fm = ec_decode(ft)`.
     If `fm < h(h+1)/2`: `itheta = (isqrt(8·fm + 1) − 1) >> 1`,
     `fs = itheta + 1`, `fl = itheta(itheta+1)/2`.
     Else: `itheta = (2(qn+1) − isqrt(8(ft − fm − 1) + 1)) >> 1`,
     `fs = qn + 1 − itheta`, `fl = ft − (qn+1−itheta)(qn+2−itheta)/2`.
     `ec_dec_update(fl, fl + fs, ft)`. isqrt = exact ⌊√x⌋ of a 32-bit
     unsigned (App. A `mathops.c` `isqrt32`).
   Then `itheta = itheta·16384 / qn` (integer division).
   (B0 is the B0 of this call: after the level-0 TF changes, or the B received
   from the parent at deeper levels and in stereo.)
6. Else if stereo (qn = 1): `inv = ec_dec_bit_logp(2)` if `b > 16` and
   `remaining_bits > 16`, else inv = 0. itheta = 0. RFC 8251 §10: a decoder
   MAY choose not to apply the phase inversion; the bit is still decoded
   (the bitstream is unchanged). The RFC 8251 "m" reference outputs are made
   that way.
7. Else (mono, qn = 1): itheta = 0, nothing read.
8. `qalloc = tell_frac − t0`; `b −= qalloc`.
9. Gains: `orig_fill = fill`. If itheta = 0: `imid = 32767`, `iside = 0`,
   `fill = fill AND (2^B − 1)`, `delta = −16384`. If itheta = 16384:
   `imid = 0`, `iside = 32767`, `fill = fill AND ((2^B − 1) << B)`,
   `delta = 16384`. Otherwise `imid = cosx(itheta)`,
   `iside = cosx(16384 − itheta)`,
   `delta = frac_mul16((N − 1) << 7, log2tan(iside, imid))`.
   `mid = imid/32768`, `side = iside/32768` (floats; note mid = 32767/32768
   when itheta = 0).

Bit-exact helpers (App. A `bands.c`, `mathops.h`; 32-bit signed integers):
* `frac_mul16(a, c) = (16384 + a·c) >> 15` (arithmetic shift; a and c fit in
  16 bits in every use).
* `cosx(x)` for x ∈ [64, 16320] (the reachable range, since itheta is a
  nonzero multiple of 16384/qn rounded down): `t = (4096 + x²) >> 13`;
  `r = (32767 − t) + frac_mul16(t, −7651 + frac_mul16(t, 8277 +
  frac_mul16(−626, t)))`; return `1 + r` (∈ [200, 32767]).
* `log2tan(s, c)`: `lc = ilog(c)`, `ls = ilog(s)`; `c = c << (15 − lc)`;
  `s = s << (15 − ls)`; return `(ls − lc)·2048 + frac_mul16(s,
  frac_mul16(s, −2597) + 7932) − frac_mul16(c, frac_mul16(c, −2597) + 7932)`
  (∈ [−15059, 15059]).

#### 8.2.5 Stereo with N = 2

If stereo and N = 2 (after theta):
1. `mbits = b`; `sbits = 8` if itheta ∉ {0, 16384}, else 0; `mbits −= sbits`.
2. `c = (itheta > 8192)`; `remaining_bits −= qalloc + sbits`.
3. `x2` = Y if c else X; `y2` = the other one.
4. If sbits > 0: `sign` = one raw bit, else 0; `s = 1 − 2·sign`.
5. `cm = band(x2, mono, N = 2, mbits, spread, B, intensity, tf_change,
   lowband, LM, lowband_out, level, gain, fill = orig_fill)` (same level and
   gain as this call; this writes lowband_out from the unit-norm x2).
6. `y2[0] = −s·x2[1]`, `y2[1] = s·x2[0]`.
7. `X[k] = mid·X[k]`, `Y[k] = side·Y[k]` for k = 0, 1; then for each k,
   `(X[k], Y[k]) = (X[k] − Y[k], X[k] + Y[k])`.
8. If inv: `Y = −Y` (unless the RFC 8251 §10 option is used). No stereo merge.

#### 8.2.6 Normal split

1. Mono only, if `B0 > 1` and `(itheta AND 0x3FFF) ≠ 0`: if itheta > 8192,
   `delta −= delta >> (4 − LM)`; else `delta = min(0, delta + ((N << 3) >>
   (5 − LM)))` (LM and N are the already-decremented/halved values).
2. `mbits = max(0, min(b, (b − delta) / 2))` (**truncating** division; b can
   be negative here); `sbits = b − mbits`; `remaining_bits −= qalloc`.
3. Sub-call parameters:
   * mono: level+1 for both; the mid (first half, X) gets lowband (if any),
     gain·mid and fill; the side (second half) gets `lowband + N` (if a
     lowband is present), gain·side and `fill >> B`; neither gets
     lowband_out; the side's cm is shifted left by `B0 >> 1`;
   * stereo: level 0 for both; the mid (X) gets lowband, lowband_out, gain
     1.0 (the mid must stay unit-norm because it is the folding source) and
     fill; the side (Y) gets no lowband, no lowband_out, gain·side and
     `fill >> B` (always 0, so the side is never folded and is zero if it has
     no pulses); the side's cm is not shifted.
   Both use the current B, LM, spread, intensity, tf_change.
4. `r0 = remaining_bits`. If `mbits ≥ sbits`: decode the mid with mbits;
   `rebalance = mbits − (r0 − remaining_bits)`; if `rebalance > 24` and
   itheta ≠ 0, `sbits += rebalance − 24`; decode the side with sbits.
   Else: decode the side first with sbits; `rebalance = sbits − (r0 −
   remaining_bits)`; if `rebalance > 24` and itheta ≠ 16384,
   `mbits += rebalance − 24`; decode the mid with mbits.
   `cm = cm_mid OR (cm_side << shift)`.

#### 8.2.7 No split: PVQ

1. `q = bits2pulses(i, LM, b)`; `cost = pulses2bits(q)`;
   `remaining_bits −= cost`. While `remaining_bits < 0` and `q > 0`:
   `remaining_bits += cost`; `q −= 1`; `cost = pulses2bits(q)`;
   `remaining_bits −= cost`.
2. If `q > 0`: `K = get_pulses(q)`; decode (§8.3) with N, K, spread, the
   current B, gain; `cm` = its collapse mask.
3. If `q = 0`: `fill = fill AND (2^B − 1)`, then
   * fill = 0: X = 0, cm = 0;
   * no lowband: for each coefficient in order, `seed = lcg(seed)`, value =
     (seed reinterpreted as signed 32-bit) >> 20 (arithmetic; −2048..2047);
     then renormalise (§8.3.3) to `gain`; `cm = 2^B − 1`;
   * lowband: for each coefficient j in order, `seed = lcg(seed)`,
     `X[j] = lowband[j] + (1/256 if (seed AND 0x8000) ≠ 0 else −1/256)`;
     renormalise to `gain`; `cm = fill`.
   `lcg(s) = (1664525·s + 1013904223) mod 2^32` (App. A `bands.c`
   `celt_lcg_rand`).

#### 8.2.8 Output stage

* Stereo: if N ≠ 2, stereo merge (§8.5) with `mid`; then if inv, `Y = −Y`
  (unless the RFC 8251 §10 option is used).
* Mono at level 0:
  1. If `B0 > 1`: `interleave(X, N_B0 >> recombine, B0 << recombine,
     hadamard = longBlocks)`.
  2. `N_B = N_B0`, `B = B0`. Repeat time_divide times: `B = B/2`;
     `N_B = 2·N_B`; `cm = cm OR (cm >> B)`; `haar(X, N_B, B)`.
  3. For k = 0..recombine−1: `cm = BIT_DEINTERLEAVE_TABLE[cm]`;
     `haar(X, N0 >> k, 1 << k)`.
  4. `B = B << recombine`.
  5. If lowband_out present: `lowband_out[j] = √N0 · X[j]` for j < N0.
  6. `cm = cm AND (2^B − 1)`.
* Mono at level > 0: nothing.

Return cm.

### 8.3 PVQ decoding, normalisation and spreading

Source: RFC §4.3.4.2–4.3.4.3 (prose); App. A `vq.c` (`alg_unquant`,
`exp_rotation`, `normalise_residual`, `extract_collapse_mask`,
`renormalise_vector`); `cwrs.c` `decode_pulses`.

#### 8.3.1 Index decoding (integer, exact)

`idx = ec_dec_uint(V(N,K))` (V < 2^32 is guaranteed by the cache). Convert to
the pulse vector y exactly as RFC §4.3.4.2 (prose) describes: k = K; for
j = 0..N−1: `p = (V(N−j−1,k) + V(N−j,k))/2`; if `idx < p` sign = +, else
sign = − and `idx −= p`; `k0 = k`; `p −= V(N−j−1,k)`; while `p > idx`:
`k −= 1`, `p −= V(N−j−1,k)`; `y[j] = sign·(k0 − k)`; `idx −= p`. (The
reference computes rows with recurrences and closed forms for N ≤ 4; the RFC
allows any exact method.)

#### 8.3.2 Normalisation (float)

`Ryy = Σ y[j]²` (exact integer); `X[j] = (gain / √Ryy) · y[j]`.

#### 8.3.3 Renormalisation (float)

`E = 1e−15 + Σ X[j]²`; `X[j] = (gain / √E) · X[j]`.

#### 8.3.4 Spreading rotation (float)

Applied after 8.3.2 with len = N (the current partition size), B (current),
K, spread:
1. If `2K ≥ N` or spread = 0: no rotation. **[PROSE ≠ CODE]** the prose does
   not mention the `2K ≥ N` exemption.
2. `f = {15, 10, 5}[spread − 1]`; `g = N / (N + f·K)`; `θ = g²/2`;
   `c = cos(π/2·θ)`, `s = cos(π/2·(1 − θ))` (= sin(π/2·θ)). The angle is
   πg²/4 as in the prose.
3. `stride2 = 0`; if `N ≥ 8·B`: `stride2 = 1`, then increment it while
   `(stride2² + stride2)·B + (B >> 2) < N` (≈ round(√(N/B))). **[PROSE ≠
   CODE]** the prose says "round(sqrt(N/nb_blocks))"; use this exact rule.
4. For each block b = 0..B−1, on the contiguous sub-vector
   `X[b·N/B .. (b+1)·N/B)` of length L = N/B, the decoder applies:
   if stride2 > 0: `rot(L, stride2, c' = s, s' = c)`; then `rot(L, 1, c, s)`.
   `rot(L, d, c, s)`: for k = 0, 1, …, L−d−1 (increasing), then for
   k = L−2d−1 down to 0: with `a = x[k]`, `e = x[k+d]`:
   `x[k+d] = c·e + s·a`, `x[k] = c·a − s·e`.
   **[PROSE ≠ CODE]** RFC §4.3.4.3 gives `x_i' = cos·x_i + sin·x_j`,
   `x_j' = −sin·x_i + cos·x_j`: the opposite sign of what the decoder applies
   (the prose's is the encoder's direction). The order of the pairs and "the
   extra rotation first" agree with the prose.
   The encoder applies the inverse before its pulse search: first
   `rot(L, 1, c, −s)`, then, if stride2 > 0, `rot(L, stride2, s, −c)`.

#### 8.3.5 Collapse mask (integer)

If B ≤ 1: cm = 1. Else with `n0 = N/B`: bit b of cm is set iff
`y[b·n0 + j] ≠ 0` for some j < n0.

### 8.4 Where tf_change acts

All TF processing happens inside `band()` at level 0 on mono calls (§8.2.2,
§8.2.8). In joint stereo the outer call is "stereo" and does none; each of its
mid and side sub-calls is a level-0 mono call that does its own. N = 1 bands
do none.

### 8.5 Stereo merge (float)

Source: App. A `bands.c` `stereo_merge`. Given the unit-norm mid X, the side
Y (already scaled by `side`), and `mid`:
1. `xp = mid · Σ X[j]·Y[j]`; `sy = Σ Y[j]²`.
2. `El = mid² + sy − 2·xp`; `Er = mid² + sy + 2·xp`.
3. If `Er < 6e−4` or `El < 6e−4`: `Y = X` (copy) and stop.
4. `X[j] = (mid·X[j] − Y[j]) / √El`; `Y[j] = (mid·X[j] + Y[j]) / √Er`
   (both from the old X[j], Y[j]).

So channel 0 = normalised(M − S), channel 1 = normalised(M + S), each of unit
norm. For intensity bands (itheta = 0, side all zero) both equal the mid.

### 8.6 RFC 8251 changes in this section

* §9 (hybrid folding): §8.1 steps 4, 5 and the `f ≥ i` bound of step 7.
  Applies in all modes but only changes anything when band start+1 is wider
  than band start (i.e. in Hybrid).
* §10 (phase inversion): optional; when chosen, the decoded `inv` is ignored
  (§8.2.5 step 8, §8.2.8).

---

## 9. Anti-collapse

Source: RFC §4.3.5 (prose, qualitative); App. A `bands.c` `anti_collapse`
(exact procedure); `celt.c` (when it runs). Float except where noted; it only
affects output samples, but its seed sequence must be followed exactly to get
float-precision output.

Runs after the final fine bits (§3.2), only if the anti-collapse flag
(§1.2 step 13) is 1 (so only for transient frames with LM ≥ 2 and enough
bits). Inputs: X/Y, collapse masks (§8.1 step 10), `pulses[]` (shape bits of
§6.9), the current energies `oldBandE` (after all fine refinements), the
histories `oldLogE` (prev1) and `oldLogE2` (prev2) from earlier frames, and a
**local copy** of the seed as it is after the band decoding of §8 (not written
back).

For i = start..end−1:
1. `N0 = width_i`; `depth = (1 + pulses[i]) / (N0 << LM)` (integer division,
   exact).
2. `thresh = 0.5·2^(−depth/8)`; `sqrt_1 = 1/√(N0 << LM)`.
3. For c = 0..C−1:
   1. `p1 = oldLogE[c][i]`, `p2 = oldLogE2[c][i]`; if C = 1:
      `p1 = max(p1, oldLogE[1][i])`, `p2 = max(p2, oldLogE2[1][i])`.
   2. `Ediff = max(0, oldBandE[c][i] − min(p1, p2))`.
   3. `r = 2·2^(−Ediff)`; if LM = 3, `r = r·1.41421356`;
      `r = min(thresh, r)·sqrt_1`.
   4. For each short block k = 0..M−1 whose bit k in `collapse_masks[i][c]`
      is 0: for j = 0..N0−1: `seed = lcg(seed)`;
      `X_c[o_i + j·M + k] = r` if `(seed AND 0x8000) ≠ 0`, else `−r`.
      (Short blocks are interleaved: coefficient j of block k is at offset
      `j·M + k` within the band.)
   5. If any block was filled, renormalise the whole band
      (`N0·M` coefficients) to unit norm (§8.3.3 with gain 1).

The seed carries over from band to band and channel to channel in this order
(i outer, c inner, blocks k, then j).

---

## 10. Synthesis

### 10.1 Energy to amplitude (with the RFC 8251 §8 cap)

Source: App. A `quant_bands.c` `log2Amp`; RFC 8251 §8. For c < C and bands
start..end−1: `lg = oldBandE[c][i] + eMeans[i]`; **RFC 8251 §8:**
`lg = min(lg, 32)`; amplitude `A[c][i] = 2^lg` (float; the reference uses
`exp(ln2·lg)` in double, rounded to float). Bands outside start..end−1 have
A = 0. `eMeans = E_MEANS` = {6.4375, 6.25, 5.75, 5.3125, 5.0625, 4.8125, 4.5,
4.375, 4.875, 4.6875, 4.5625, 4.4375, 4.875, 4.625, 4.3125, 4.5, 4.375, 4.625,
4.75, 4.4375, 3.75, 3.75, 3.75, 3.75, 3.75} (25 entries, 21 used). These are
exactly `E_MEANS_Q4[i]/16` (the fixed-point table {103,100,92,85,81,77,72,70,
78,75,73,71,78,74,69,72,70,74,76,71,60,…}; the source comment calls it "Q6"
but the scale relative to log2 units is 1/16).

Then the silence override of §1.2 is applied if needed.

### 10.2 Denormalisation

Source: RFC §4.3.6 (prose); App. A `bands.c` `denormalise_bands`, `celt.c`.
For c < C, every bin j of band i < effEnd: `F_c[j] = X_c[j] · A[c][i]`;
bins ≥ `M·eBands[effEnd]` are 0. Then bins below `M·eBands[start]` are set to
0, and so are bins ≥ `bound`, where `bound = M·eBands[effEnd]`, reduced to
`min(bound, N/downsample)` when downsample ≠ 1 (this band-limits the signal
for decimation, §10.6).

The signal scale is that of 16-bit PCM (full scale ≈ 32768).

### 10.3 Channel mapping

If CC = 2 and C = 1: `F_1 = F_0` (copy). If CC = 1 and C = 2:
`F_0[j] = 0.5·(F_0[j] + F_1[j])`. (Done on the MDCT coefficients, before the
inverse MDCT.)

### 10.4 Inverse MDCT, windowing and overlap-add

Source: RFC §4.3.7 (prose); App. A `mdct.c` `clt_mdct_backward`, `celt.c`
`compute_inv_mdcts`. Float; output only. The exact relation below was
verified numerically against the reference procedure.

Window: overlap `L = 120`; `W[i] = WINDOW120[i]` ≈
`sin(π/2 · sin²(π/2 · (i + 1/2)/120))`, i < 120 (RFC §4.3.7 formula with
L = 120; App. A `window120`). The float table equals the float rounding of
this formula except for 5 entries that differ by one float ulp; use the table
for the closest match.

Inverse transform of one block of n coefficients `Z[0..n)` (n = 120·M for a
long block, 120 for a short block): define, for m = 0 .. n + L − 1,
`m' = m + (n − L)/2` and

  `y[m] = g · w(m) · Σ_{k=0}^{n−1} Z[k] · cos(π/n · (m' + 1/2 + n/2) · (k + 1/2))`

with `w(m) = W[m]` for m < L, 1 for L ≤ m < n, `W[n + L − 1 − m]` for
n ≤ m < n + L, and `g = 1 + s²`, `s = π/(8n)`. (This is the textbook
2n-point IMDCT, without any 1/n factor, of which only the n + L samples under
the low-overlap window are kept. The factor g = 1 + s², about 1 + 10⁻⁵ to
1 + 2·10⁻⁷, comes from the reference's first-order approximation of the
π/(8n) twiddle rotation; include it for float-precision agreement.)

Per output channel c, with a work buffer `x` of N + L samples, all zero:
* long block (not transient): one transform of `Z[k] = F_c[k]`, k < N; add
  y[m] into x[m].
* transient: for b = 0..M−1, one transform of n = 120 coefficients
  `Z[k] = F_c[b + k·M]`, adding y[m] into `x[120·b + m]` (consecutive short
  blocks overlap by L and add).

Then, with the per-channel persistent overlap memory `ov_c[0..L)`:
`out[j] = x[j] + ov_c[j]` for j < L, `out[j] = x[j]` for L ≤ j < N;
`ov_c[j] = x[N + j]` for j < L. `out[0..N)` is the frame's time signal before
the post-filter.

**Note on the prose:** RFC §4.3.7 says the IMDCT "scales by 1/2"; the
normative behaviour is the formula above. An encoder's forward MDCT only has
to be the transform this inverse reconstructs (time-domain aliasing
cancellation with the same window, matching scale); it never affects parsing.

### 10.5 Post-filter

Source: RFC §4.3.7.1 (prose: parameters, taps, window-squared cross-fade);
App. A `celt.c` `comb_filter` and its calls. Float; output only.

Tap gains `COMB_FILTER_GAINS[tapset]` = {0.3066406250, 0.2170410156,
0.1296386719}, {0.4638671875, 0.2680664062, 0}, {0.7998046875, 0.1000976562,
0} (equal to the prose).

`comb(seg, n, (T0, g0, tap0) → (T1, g1, tap1), ov)` acts **in place** on a
segment of the channel's continuous post-filtered history: let
`a_k = g0·GAINS[tap0][k]`, `b_k = g1·GAINS[tap1][k]`. For i < ov,
`f = W[i]²`:
`seg[i] += (1−f)·a0·seg[i−T0] + (1−f)·a1·(seg[i−T0−1] + seg[i−T0+1])
          + (1−f)·a2·(seg[i−T0−2] + seg[i−T0+2])
          + f·b0·seg[i−T1] + f·b1·(seg[i−T1−1] + seg[i−T1+1])
          + f·b2·(seg[i−T1−2] + seg[i−T1+2])`;
for ov ≤ i < n: `seg[i] += b0·seg[i−T1] + b1·(seg[i−T1−1] + seg[i−T1+1]) +
b2·(seg[i−T1−2] + seg[i−T1+2])`.
Because it is in place and T ≥ 15, every `seg[i−T±k]` read is an
already-filtered sample (of this frame or of earlier frames): the filter is
recursive, `y(n) = x(n) + G·(g0·y(n−T) + g1·(y(n−T−1) + y(n−T+1)) +
g2·(y(n−T−2) + y(n−T+2)))`. **[PROSE ≠ CODE]** RFC §4.3.7.1's formula
writes `y(n−T+1)+y(n−T+1)` and `y(n−T+2)+y(n−T+2)`: a typo for the ±1, ±2
taps. In the reference each weighted product `(1−f)·a_k` (or `f·b_k`) is
formed first and then multiplied by the sample. (If T1 = 0 with g1 = 0, the
b terms read not-yet-filtered samples but contribute exactly 0.)

Per frame, with decoder state (period, gain, tapset) "current" `P` and "old"
`P_old`, and the newly decoded parameters `P_new` (§1.2 step 2; period 0 and
gain 0 when absent), for each output channel:
1. `P.period = max(P.period, 15)`, `P_old.period = max(P_old.period, 15)`.
2. `comb(out[0..120), 120, P_old → P, ov = 120)`.
3. If LM ≠ 0: `comb(out[120..N), N − 120, P → P_new, ov = 120)`.
Then `P_old = P`, `P = P_new`; and if LM ≠ 0, also `P_old = P_new`.
(So for 2.5 ms frames the new parameters take effect one frame later; for
longer frames the first 2.5 ms use the previous frame's filter and the next
2.5 ms cross-fade to the new one. The filter is always run, also with zero
gains.)

### 10.6 De-emphasis, decimation, output scaling

Source: RFC §4.3.7.2 (prose: 1/(1 − 0.8500061035·z⁻¹)); App. A `celt.c`
`deemphasis`. Float; output only.

Per output channel with persistent state `m` (float): for j = 0..N−1:
`t = out[j] + m`; `m = 0.8500061·t` (`PREEMPH[0]`, float 0.85000610);
if `j mod downsample == 0`, output sample `j/downsample` = `t / 32768`.
Output rates below 48 kHz are produced by this plain decimation (keep every
downsample-th sample), relying on the spectral band limit of §10.2; there is
no other filter. Float output is in [−1, 1) scale; 16-bit output would be the
rounded and saturated `t`.

---

## 11. Decoder state

Source: App. A `celt.c` (`OpusCustomDecoder`, `celt_decode_with_ec`,
`OPUS_RESET_STATE`, `celt_decode_lost`), `src/opus_decoder.c`.

### 11.1 Contents and reset values

* `oldBandE[2][21]` (log2 energies relative to eMeans): reset 0.
* `oldLogE[2][21]`, `oldLogE2[2][21]` (energy history for anti-collapse):
  reset −28.
* `backgroundLogE[2][21]` (concealment only): reset 0.
* `rng` (32-bit seed): reset 0.
* post-filter: `P = P_old = (period 0, gain 0, tapset 0)`: reset 0.
* de-emphasis memory per channel: reset 0.
* per channel: time history ≥ 1024 + 2 + N samples of post-filtered output (the
  reference keeps 2048), and the IMDCT overlap memory `ov` (120): reset 0.
* concealment: loss count, last pitch, LPC: reset 0.

Two channels of energy state are always kept, whatever C is.

The Opus layer resets the CELT state (App. A `opus_decoder.c`): at
initialisation; before decoding the CELT part of a packet whose mode differs
from the previous packet's mode (unless the previous packet carried a
SILK→CELT redundant frame); and before decoding a SILK→CELT redundant frame.

### 11.2 Updates by a decoded frame (in this order)

1. At the start: if C = 1, `oldBandE[0][i] = max(oldBandE[0][i],
   oldBandE[1][i])` for all 21 bands.
2. Coarse, fine, final-fine energies update `oldBandE[c][start..end)` for
   c < C (§2, §3).
3. Folding updates the seed (§8.2.7); anti-collapse uses a copy (§9).
4. Silence: `oldBandE[c][i] = −28` for c < C, all i.
5. History shift: the per-channel time history advances by N; the new frame
   occupies its last N samples (and is post-filtered in place, §10.5).
6. Post-filter parameters as in §10.5.
7. If C = 1: `oldBandE[1][i] = oldBandE[0][i]` for all i.
8. If not transient: `oldLogE2 = oldLogE`, `oldLogE = oldBandE`,
   `backgroundLogE = min(backgroundLogE + M·0.001, oldBandE)` (all 2×21
   entries). If transient: `oldLogE = min(oldLogE, oldBandE)` (entrywise; the
   older history is kept).
9. For both channels c = 0, 1 (regardless of C), for bands i < start and
   i ≥ end: `oldBandE[c][i] = 0`, `oldLogE[c][i] = oldLogE2[c][i] = −28`.
10. `rng` = the range decoder's final `rng` (this overwrites the folding
    updates of step 3; the next frame's folding starts from it).
11. Loss count = 0; de-emphasis memories as in §10.6; overlap memories as in
    §10.4.

### 11.3 Mono/stereo switches

A mono stream decoded to stereo output duplicates the spectrum (§10.3); a
stereo stream decoded to mono averages it. Energy state: steps 1 and 7 keep
the two channels' energies consistent when C changes between frames;
anti-collapse with C = 1 uses the maximum over both channels' histories
(§9 step 3.1).

### 11.4 Packet loss (non-normative; summary)

The reference conceals a lost frame (also used for `len ≤ 1`) by:
* after 5 or more consecutive losses, or in Hybrid (start ≠ 0): noise
  shaped by energies — the energies are decayed in `oldBandE` (by 1.5 on the
  first loss, 0.5 afterwards, log2 units, bands start..end−1) or taken from
  `backgroundLogE` after 5 losses; each band gets LCG noise (seed from `rng`,
  `rng` updated) normalised to unit norm, then denormalised and inverse
  transformed as a long block;
* otherwise: pitch-based extrapolation of the time-domain history (pitch
  search, 24th-order LPC, periodic excitation with decay, energy guard,
  windowed TDAC blend with the overlap memory, post-filter/pre-filter of the
  overlap with the current post-filter parameters);
* de-emphasis as usual; loss count incremented.

What matters for later normal decoding: `oldBandE` decay (noise branch),
`rng`, the time history and overlap memory, the de-emphasis memory, and the
loss count (0 again after a good frame). Post-filter parameters and
`oldLogE`/`oldLogE2` are not changed by concealment.

---

## 12. Encoder obligations

What an encoder must mirror so that this decoder reads its stream as
intended (App. A `celt.c` `celt_encode_with_ec`, `quant_bands.c`, `bands.c`
with `encode = 1`, `rate.c`, `laplace.c` `ec_laplace_encode`, `vq.c`
`alg_quant`, `cwrs.c` `encode_pulses`):

* The same symbol order, the same gates (computed from the encoder's
  `ec_tell`/`ec_tell_frac`, which equal the decoder's), and the same
  allocation computation (§6) with the same reservations, using the same
  boosts, trim, tf, spread values it codes. The only free choices in the
  allocation are the skip decisions (§6.6) and intensity ≤ codedBands /
  dual (§6.7).
* Coarse energy: the coded qi must be representable with the chosen symbol:
  with ≥ 15 bits left, Laplace; the Laplace encoder clamps very large |qi| to
  the largest value representable in the tail (and the encoder must use the
  clamped value for its prediction state); with 2..14 bits, qi is clamped to
  [−1, 1]; with 1 bit, to min(qi, 0) and coded as a bit (1 = −1); with 0 bits
  qi = −1. The reference additionally (for i ≠ start, with
  `bits_left = budget − tell − 3·C·(end − i)`) forces qi ≤ 1 if
  bits_left < 24 and qi ≥ −1 if bits_left < 16, and limits how fast energy
  may drop — encoder policy, not normative.
* Fine energy: `q2 = clamp(⌊(error + 0.5)·2^fq⌋, 0, 2^fq − 1)` (reference);
  the decoder reconstructs per §3.
* Final fine bits: same priority loop; the bit is 1 if the remaining error
  is ≥ 0.
* tf: §4 (encoder replaces uncodable flags by the running value).
* Theta (encoder): `itheta = ⌊(itheta_exact·qn + 8192) / 16384⌋` (integer, the
  reference computes itheta_exact from atan2 of the band norms); for stereo
  with qn = 1 the encoder decides `inv` and codes it only under the decoder's
  condition. For N = 2 stereo the sign bit is the sign of the cross product.
  The rebalancing and the order of mid/side coding must follow §8.2.6.
* PVQ: any search giving a vector with exactly K pulses; the encoder applies
  the inverse spreading rotation before the search (§8.3.4).
* Encoders need not resynthesise: folding, collapse masks and anti-collapse
  never affect parsing. The anti-collapse flag is a free choice (when
  reserved).
* Post-filter: period 15..1022 must be coded as octave/fine pitch per §1.2;
  gain index and tapset are free.

---

## 13. Constants

| constant | value | source |
|---|---|---|
| BITRES | 3 | App. A `entcode.h` |
| frame sizes / LM | 120·2^LM, LM 0..3 | RFC §4.3 Table 55 |
| nbEBands / effEBands | 21 / 21 | App. A `static_modes_float.h` |
| start band (Hybrid) | 17 | RFC §4.3 prose; App. A `opus_decoder.c` |
| end band per bandwidth | NB 13, MB/WB 17, SWB 19, FB 21 | App. A `opus_decoder.c` |
| silence logp | 15 | RFC Table 56; App. A `celt.c` |
| post-filter gate | tell + 16 ≤ total | App. A `celt.c` |
| octave | `ec_dec_uint(6)` (0..5) | App. A `celt.c` |
| period | (16<<oct) + raw(4+oct) − 1 | RFC §4.3.7.1; App. A `celt.c` |
| post-filter gain | 0.09375·(qg+1) | RFC §4.3.7.1; App. A `celt.c` |
| tapset gate / icdf | tell + 2 ≤ total; {2,1,0}, ftb 2 | App. A `celt.c` |
| COMBFILTER_MINPERIOD / MAXPERIOD | 15 / 1024 | App. A `celt.c` |
| comb gains | §10.5 | RFC §4.3.7.1; App. A `comb_filter` |
| transient gate / logp | LM > 0, tell + 3 ≤ total / 3 | App. A `celt.c`; RFC Table 56 |
| intra gate / logp | tell + 3 ≤ total / 3 | same |
| Laplace MINP / NMIN | 1 / 16 | App. A `laplace.c` |
| Laplace model scaling | fs = p<<7, decay = d<<6, k = min(i,20) | App. A `quant_bands.c` |
| coarse thresholds | 15 / 2 / 1 bits left | App. A `quant_bands.c` |
| energy floor in prediction | −9 | App. A `quant_bands.c` |
| α (inter), β (inter) | {29440,26112,21248,16384}/32768, {30147,22282,12124,6554}/32768 | App. A `quant_bands.c`; RFC §4.3.2.1 (β intra) |
| α, β (intra) | 0, 4915/32768 | RFC §4.3.2.1; App. A |
| MAX_FINE_BITS | 8 | App. A `rate.h` |
| FINE_OFFSET | 21 | App. A `rate.h` |
| QTHETA_OFFSET / _TWOPHASE | 4 / 16 | App. A `rate.h` |
| MAX_PSEUDO / LOG_MAX_PSEUDO | 40 / 6 | App. A `rate.h` |
| ALLOC_STEPS | 6 | App. A `rate.c` |
| tf logp | 2/4 first band, 4/5 others (transient/not) | RFC §4.3.4.5; App. A `tf_decode` |
| spread gate, default | tell + 4 ≤ total, 2 | App. A `celt.c` |
| dynalloc logp start / min | 6 / 2 | RFC §4.3.3; App. A |
| boost quanta | min(8w, max(48, w)) | RFC §4.3.3; App. A |
| trim gate / default | tell_frac + 48 ≤ total / 5 | RFC §4.3.3; App. A |
| anti-collapse rsv | 8 if transient, LM ≥ 2, bits ≥ 8(LM+2) | App. A `celt.c` |
| skip / dual rsv | 8 if total ≥ 8 | App. A `rate.c` |
| thresh | max(8C, (24·N·M) >> 4) | RFC §4.3.3; App. A |
| split threshold | b > cache[cache[0]] + 12 | App. A `bands.c` |
| theta qb limit | min(64, …), qn = 1 if qb < 4 | App. A `compute_qn` |
| rebalance threshold | 24 (3 bits) | App. A `bands.c` |
| per-band b cap | 16383 | App. A `quant_all_bands` |
| inv bit gate / logp | b > 16 and remaining > 16 / 2 | App. A `bands.c` |
| mid/side gains for 0/16384 | 32767/0, delta ∓16384 | App. A `bands.c` |
| cosx coefficients | 4096, 13, 32767, −7651, 8277, −626 | App. A `bitexact_cos` |
| log2tan coefficients | 2048, −2597, 7932 | App. A `bitexact_log2tan` |
| LCG | 1664525, 1013904223 | App. A `celt_lcg_rand` |
| noise fill value | (int32)seed >> 20 | App. A `quant_band` |
| fold dither | ±1/256, bit 0x8000 | App. A `quant_band` |
| SPREAD_FACTOR | {15, 10, 5} | RFC Table 59; App. A `vq.c` |
| stride2 rule | (s²+s)·B + B/4 < N | App. A `exp_rotation` |
| normalisation epsilon | 1e−15 | App. A `arch.h` EPSILON |
| stereo merge floor | 6e−4 | App. A `stereo_merge` |
| Haar coefficient | 0.70710678 | App. A `haar1` |
| anti-collapse | 0.5·2^(−depth/8), 2·2^(−Ediff), ×1.41421356 at LM = 3 | App. A `anti_collapse` |
| energy history reset | −28 | App. A `celt.c` |
| background increment | 0.001·M per frame | App. A `celt.c` |
| energy cap | lg ≤ 32 | RFC 8251 §8 |
| overlap | 120 | App. A mode |
| IMDCT gain factor | 1 + (π/(8n))² | derived from App. A `mdct.c` (§10.4) |
| de-emphasis | 0.85000610 (float of 0.8500061035) | RFC §4.3.7.2; App. A `PREEMPH` |
| output scale | 1/32768 | App. A `arch.h` (SCALEOUT) |
| downsample factors | 48k:1, 24k:2, 16k:3, 12k:4, 8k:6 | App. A `resampling_factor` |
| max len | 1275 bytes | App. A `celt.c` |

---

## 14. RFC prose vs. normative code (code wins)

1. §4.3.3 dynalloc loop bound: prose compares with the original frame size
   (total_bits + total_boost); code compares with the boost-reduced total
   (§5.2).
2. §4.3.3 `skip_rsv` and `dual_stereo_rsv`: prose "greater than 8"; code
   "≥ 8" (§6.2).
3. §4.3.3 intensity reservation: prose "ilog2(end−start) bits"; code
   `LOG2_FRAC_TABLE[end−start]` eighth bits (§6.2).
4. §4.3.3 trim offsets: prose "divide by 64"; code arithmetic right shift by
   6 (floor for negative values) (§6.3).
5. §4.3.4.3 rotation sign: prose gives the encoder-direction rotation; the
   decoder rotates the opposite way (§8.3.4).
6. §4.3.4.3 extra-rotation stride: prose round(√(N/B)); code's exact
   incremental rule (§8.3.4); and the prose omits the "no rotation when
   2K ≥ N" rule.
7. §4.3.7.1 octave range: prose "0 to 6"; code `ec_dec_uint(6)` = 0..5.
8. §4.3.7.1 post-filter formula: prose taps `y(n−T+1)+y(n−T+1)`,
   `y(n−T+2)+y(n−T+2)`; code uses n−T∓1 and n−T∓2 (§10.5).
9. Table 56 lists the anti-collapse flag with PDF {1,1}/2; it is a raw bit
   (§1.2 step 13). (Likewise the fine/final energy bits, post-filter fine
   pitch and gain, N = 1 and N = 2 stereo signs are raw bits.)
10. §4.3.2.2 final fine bits: the per-band `bits_left ≥ C` condition is only
    in the code (§3.2).
11. §4.3.7 IMDCT "scaling by 1/2": the normative scaling is given exactly in
    §10.4.
12. §4.4 says `celt_decode_lost()` is in `mdct.c`; it is in `celt.c`
    (informational).
13. `eMeans`' source comment says "Q6"; the values are in 1/16 log2 units
    (§10.1) (informational).

---

## 15. Tables

`tools/appendix_a_tables.py rfc6716.txt` extracts Appendix A in memory
exactly as RFC §A.1 describes (checking the RFC's SHA-1 of the archive) and
prints, as Rust `pub const` items with a comment naming file and symbol:

| Rust name | source | derivable? |
|---|---|---|
| `EBAND5MS` | `modes.c` `eband5ms` | = RFC Table 55 |
| `BAND_ALLOCATION` [11][21] | `modes.c` `band_allocation` | = RFC Table 57 (transposed) |
| `E_MEANS`, `E_MEANS_Q4` | `quant_bands.c` `eMeans` (float / fixed) | E_MEANS = E_MEANS_Q4/16 exactly |
| `PRED_COEF`, `PRED_COEF_Q15`, `BETA_COEF`, `BETA_COEF_Q15`, `BETA_INTRA`, `BETA_INTRA_Q15` | `quant_bands.c` | float = Q15/32768 exactly |
| `E_PROB_MODEL` [4][2][42] | `quant_bands.c` `e_prob_model` | no (data) |
| `SMALL_ENERGY_ICDF`, `TRIM_ICDF`, `SPREAD_ICDF`, `TAPSET_ICDF` | `quant_bands.c`, `celt.c` | = RFC Tables 56/58 PDFs |
| `TF_SELECT_TABLE` [4][8] | `celt.c` `tf_select_table` | = RFC Tables 60–63 |
| `COMB_FILTER_GAINS` [3][3] | `celt.c` `comb_filter` `gains` | = RFC §4.3.7.1 |
| `LOG2_FRAC_TABLE` [24] | `rate.c` | yes: ⌈8·log2(i+1)⌉ = log2_frac(i+1,3), exact |
| `ORDERY_TABLE` [30] | `bands.c` | yes (§8.2.2), exact |
| `BIT_INTERLEAVE_TABLE`, `BIT_DEINTERLEAVE_TABLE` [16] | `bands.c` | yes (§8.2.2), exact |
| `EXP2_TABLE8` [8] | `bands.c` `compute_qn` | yes: ⌊16384·2^(k/8)⌋, exact |
| `SPREAD_FACTOR` [3] | `vq.c` | = RFC Table 59 |
| `WINDOW120` [120] | `static_modes_float.h` `window120` | formula of RFC §4.3.7; 115/120 entries equal the float rounding, 5 differ by 1 ulp |
| `LOGN400` [21] | `static_modes_float.h` `logN400` | yes: log2_frac(width,3) = ⌈8·log2 width⌉, exact |
| `CACHE_INDEX50` [5][21], `CACHE_BITS50` [392] | `static_modes_float.h` | yes (§7.2, needs log2_frac exactly), exact |
| `CACHE_CAPS50` [8][21] | `static_modes_float.h` | yes (§7.3), exact |
| `PREEMPH` [4] | `static_modes_float.h` mode struct | coefficient of RFC §4.3.7.2 |

Not extracted (output-only, derivable): the FFT/MDCT twiddles
(`cos(2πi/N)`), FFT bit-reverse tables, and the fixed-point helper
coefficients (not used by the float build). `INV_TABLE` of `cwrs.c` is an
implementation device for V(N,K) and is not needed.

---

## Appendix: range decoder primitives used

All defined in RFC §4.1 (prose) and implemented in App. A `entdec.c`,
`entcode.c`: `ec_decode(ft)` / `ec_decode_bin(bits)` + `ec_dec_update(fl, fh,
ft)` (§4.1.2, §4.1.3.1), `ec_dec_bit_logp(logp)` (§4.1.3.2),
`ec_dec_icdf(icdf, ftb)` (§4.1.3.3), `ec_dec_bits(n)` raw bits from the end
(§4.1.4), `ec_dec_uint(ft)` (§4.1.5; values with ft − 1 ≥ 2^8 use 8
range-coded bits plus raw bits), `ec_tell()` = `nbits_total − ilog(rng)` and
`ec_tell_frac()` (§4.1.6). Raw-bit reads reduce the space available to the
range coder and are counted by `ec_tell`. Reads past the end of the data
return zero bits.
