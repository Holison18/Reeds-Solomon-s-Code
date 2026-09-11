# Reed-Solomon's Code Implementation

A from-scratch implementation of Reed-Solomon error correction over GF(256), the same field/parameters used by QR codes, plus a visual demo that runs a real image through it.

## Contents

| File | Purpose |
|---|---|
| [reed_solomon.py](reed_solomon.py) | The codec itself: GF(256) arithmetic, encoder, and error/erasure decoder. Running it directly executes a self-test against a worked QR example. |
| [image_demo.py](image_demo.py) | Encodes an image's bytes with the codec, corrupts them, and reconstructs the image with and without correction, saving before/after images and charts. |
| `demo_output/` | Created by `image_demo.py`: the generated PNGs and charts. |

## Requirements

- Python 3.9+
- `numpy`, `Pillow`, `matplotlib` (only needed for `image_demo.py` — `reed_solomon.py` has no dependencies beyond the standard library)

Install the demo dependencies:

```bash
pip install numpy Pillow matplotlib
```

## How to run

### 1. The codec's self-test

```bash
python reed_solomon.py
```

This encodes a known QR data block, checks the computed parity bytes against the article's worked example, injects 5 random byte errors and recovers them, then repeats the exercise with 10 known erasures. Every step ends in an `assert`, so a silent, clean exit means everything checked out.

### 2. The image corruption/recovery demo

```bash
python image_demo.py --image path/to/your/photo.jpg
```

`--image` is required — this runs your image through the codec, corrupts ~5% of the encoded bytes, and writes the results to `demo_output/`. Useful flags:

```bash
python image_demo.py --image photo.jpg --scatter-rate 0.02 --burst-fraction 0.01   # tune how much damage is introduced
python image_demo.py --image photo.jpg --nsym 32                                   # parity bytes per 255-byte block (default 32)
```

Output files (`demo_output/`):

- `1_original.png`, `2_corrupted_no_correction.png`, `3_recovered_with_rs.png` — the three images on their own
- `4_comparison.png` — the three side by side with corruption/correction stats in the title (the one worth putting in a slide)
- `5_success_rate_curve.png` — recovery success rate vs. corruption rate, with a line marking the theoretical correction limit

## How Reed-Solomon codes work

The explanation and terminology below follow James S. Plank's article ["A Tutorial on Reed-Solomon Coding for Fault-Tolerance in RAID-like Systems"](https://www.cs.cmu.edu/~guyb/realworld/reedsolomon/reed_solomon_codes.html), hosted on CMU's "Algorithms in the Real World" course page. Each step is mapped to the function in [reed_solomon.py](reed_solomon.py) that implements it.

### Specification

A Reed-Solomon code is denoted **RS(n, k)** with *s*-bit symbols (here, bytes, so *s = 8* and the field is GF(256)). The encoder takes *k* data symbols and appends *n − k* parity symbols to form an *n*-symbol codeword. A decoder can correct up to *t* errors, where **2t = n − k** — every 2 parity symbols buy you 1 correctable error. In this repo, `nsym = n - k` is the number of parity bytes, so the correction capacity is `nsym // 2`.

### Encoding

**Generator polynomial.** The codeword is constructed as *c(x) = g(x)·i(x)*, where *i(x)* is the information (data) block and *g(x)* is the generator polynomial, built from the roots of the field's primitive element *a*:

```
g(x) = (x - a^0)(x - a^1) ... (x - a^(2t-1))
```

→ `rs_generator_poly(nsym)`

**Systematic encoding.** In a systematic code "the data is left unchanged and the parity symbols are appended" after it, rather than the message being transformed into some other representation. The 2t parity symbols are computed via finite-field polynomial division of the data against the generator polynomial.

→ `rs_encode_msg(msg_in, nsym)` — performs the polynomial division in GF(256) and appends the remainder as parity, restoring the original data bytes afterward (the division process overwrites them as a side effect).

### Decoding

A received codeword is modeled as *r(x) = c(x) + e(x)*, i.e. the original codeword plus whatever error pattern corrupted it. Decoding proceeds in four stages:

**1. Syndrome calculation.** "A Reed-Solomon codeword has 2t syndromes that depend only on errors, not on the transmitted codeword." They're computed by evaluating *r(x)* at each root of the generator polynomial (`a^0, a^1, ..., a^(2t-1)`). If every syndrome is zero, the codeword is either clean or the errors are undetectable; if any is nonzero, errors are present and their signature is now isolated in the syndromes.

→ `rs_calc_syndromes(msg, nsym)`

**2. Finding error locations.** This is a two-step process:
- Solve for the **error locator polynomial** — this repo uses the Berlekamp-Massey algorithm (the article also mentions Euclid's algorithm as an alternative). If erasure positions are already known, their locator polynomial seeds the algorithm so it only has to search for the remaining, unlocated errors.
- Find the **roots** of that polynomial via the **Chien search** — a trial substitution of every codeword position, since a root corresponds to an actual error location.

→ `rs_find_error_locator(synd, nsym, erase_loc, erase_count)` (Berlekamp-Massey), `rs_find_errors(err_loc, nmess)` (Chien search)

**3. Finding error values.** Once the error *positions* are known, the **Forney algorithm** computes the error *magnitude* at each position — the byte value that must be XORed in to correct it — from the syndromes and the error locator polynomial.

→ `rs_correct_errata(msg, synd, pos)`

**4. Correction.** The computed magnitudes are applied at their positions, and the syndromes are recomputed on the corrected message as a final check that it's now a valid codeword.

→ `rs_correct_msg(msg_in, nsym, erase_pos)` — ties all four stages together and is the main decode entry point.

### Outcomes

Decoding one codeword ends in one of three ways, governed by **2s + r < 2t** (where *s* is the number of errors and *r* the number of erasures):
1. **Successful recovery** of the original codeword (this repo raises nothing and returns the corrected message).
2. **Detected but uncorrectable errors** — too many errors/erasures to fix (this repo raises `ValueError`, both up front when the error/erasure count is provably too high, and after correction if the result still isn't a valid codeword).
3. **Silent miscorrection** — an undetected error, i.e. the decoder "corrects" to the wrong codeword without noticing. This is a fundamental limitation of the math, not a bug: with enough errors, a corrupted codeword can look like a valid — but wrong — one. It's why `image_demo.py`'s success-rate chart shows a hard cliff at the theoretical limit rather than a graceful decline.
