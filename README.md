# Reed-Solomon Codec: Replication, Channel Simulation, and Sensitivity Analysis

**Course:** COE 592: Advanced Signal and Communication Theory  
**Topic:** Replication and Empirical Assessment of Reed-Solomon Error-Correcting Codes over $\text{GF}(256)$

---

## Overview

This repository provides an independent, from-scratch implementation of a **Reed-Solomon $(\text{RS})$ error-correcting codec** over Galois Field $\text{GF}(256)$. The project tests and replicates the central theoretical claims of Reed-Solomon codes—specifically the **Maximum Distance Separable (MDS)** property and the **Singleton bound** ($d_{\min} = n - k + 1$)—under both unknown errors ($2s \le n - k$) and known erasures ($r \le n - k$).

### What This Implementation Accomplishes

1. **Galois Field $\text{GF}(256)$ Engine:** Implements finite field arithmetic via log/exp lookup tables using primitive polynomial $p(x) = x^8 + x^4 + x^3 + x^2 + 1$ (`0x11D`).
2. **Algebraic Decoding Pipeline:** Complete 4-stage decoding architecture incorporating syndrome evaluation, the Berlekamp-Massey algorithm, Chien search, and Forney algorithm for error/erasure recovery.
3. **QR Standard Verification:** Validates the encoder and decoder against standard worked examples from ISO/IEC 18004.
4. **Information Loss Prevention (Text & Image):**
   - **Text Transmission:** Demonstrates zero-loss recovery of arbitrary string messages (e.g., `"He is a boy"`) under character noise and dropped packets.
   - **Image Channel Simulation:** Simulates noisy transmission of raw RGB images corrupted with both scattered byte errors and contiguous burst scratches.
5. **Beyond-the-Paper Sensitivity Study:** Generates empirical performance waterfalls and verifies the $2s + r \le n - k$ MDS trade-off frontier.

---

## Repository Structure

| File / Directory | Description |
|---|---|
| [`reed_solomon.py`](reed_solomon.py) | Pure-Python $\text{GF}(256)$ codec: arithmetic, generator polynomial, systematic encoder, Berlekamp-Massey, Chien search, and Forney decoder. Includes self-tests and text recovery demo. |
| [`image_demo.py`](image_demo.py) | Simulates image transmission across a degraded channel (random scatter + burst scratches) and reconstructs naive vs. RS-corrected output. |
| [`extended_analysis.py`](extended_analysis.py) | Extended sensitivity analysis: sweeps code rates ($nsym \in \{8, 16, 32, 64\}$) and maps the empirical MDS error/erasure trade-off frontier. |
| [`requirements.txt`](requirements.txt) | Environment dependencies (`numpy`, `Pillow`, `matplotlib`). |
| [`demo_output/`](demo_output/) | Generated visual artifacts, comparison collages, and empirical waterfall curves. |

---

## Setup & Requirements

- Python 3.9+
- Install dependencies:

```bash
pip install -r requirements.txt
```

---

## How to Run & Reproduce Results

### 1. Codec Self-Test & Text Recovery Demo

```bash
python reed_solomon.py
```
- Validates encoder parity against the published QR standard worked example.
- Verifies syndrome calculation on clean codewords.
- Tests recovery from 5 random byte errors and 10 known erasures.
- Runs an interactive text string recovery demonstration (`"He is a boy"`).

### 2. Image Channel Corruption & Recovery Demo

```bash
python image_demo.py --image cat-1045782_640.jpg
```
- Splits the image into $\text{RS}(255, 223)$ blocks ($nsym = 32$).
- Injects $3\%$ random scattered byte corruptions and a $2\%$ contiguous burst scratch.
- Saves comparison figures to `demo_output/` showing naive vs. Reed-Solomon reconstruction.

### 3. Extended Sensitivity & MDS Bound Analysis

```bash
python extended_analysis.py
```
- Generates multi-curve sensitivity plots comparing error tolerance across code rates.
- Evaluates the $2s + r \le n - k$ MDS capacity frontier with empirical simulation points.

---

## Key Generated Visualizations

All figures in `demo_output/` are generated directly from the implementation code:

- `4_comparison.png` — Three-panel side-by-side comparison: Original, Corrupted without correction, and Recovered with Reed-Solomon.
- `5_success_rate_curve.png` — Block recovery success rate vs. channel corruption rate for $\text{RS}(255, 223)$, highlighting the theoretical correction limit.
- `6_sensitivity_analysis.png` — Multi-rate waterfall curves comparing $nsym \in \{8, 16, 32, 64\}$.
- `7_mds_bound_verification.png` — Empirical verification of the Singleton MDS trade-off frontier ($2s + r \le 32$).

---

## References & Additional Resources

- **Primary Tutorial & Theoretical Guide:**  
  James S. Plank, [*A Tutorial on Reed-Solomon Coding for Fault-Tolerance in RAID-like Systems*](https://www.cs.cmu.edu/~guyb/realworld/reedsolomon/reed_solomon_codes.html), Technical Report UT-CS-96-332, University of Tennessee / CMU "Algorithms in the Real World".
- **Foundational Paper:**  
  Irving S. Reed and Gustave Solomon, *"Polynomial Codes Over Certain Finite Fields"*, Journal of the Society for Industrial and Applied Mathematics (SIAM), Vol. 8, No. 2, pp. 300–304, 1960.
- **Decoding Algorithms:**  
  - E. R. Berlekamp, *Algebraic Coding Theory*, McGraw-Hill, 1968.  
  - J. L. Massey, *"Shift-register synthesis and switch-circuit analysis"*, IEEE Transactions on Information Theory, 1969.  
  - G. D. Forney, *"On decoding BCH codes"*, IEEE Transactions on Information Theory, 1965.
- **Standard Application:**  
  ISO/IEC 18004: *Information technology — Automatic identification and data capture techniques — QR Code bar code symbology specification*.
