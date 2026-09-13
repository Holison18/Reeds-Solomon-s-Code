"""
Extended Analysis and Sensitivity Study for Reed-Solomon Codec
COE 592: Advanced Signal and Communication Theory

This script produces publication-quality figures that extend beyond
the baseline paper implementation:
1. Sensitivity Analysis: Code Rate vs Channel Error Rate (nsym = 8, 16, 32, 64)
2. Trade-off Curve: Theoretical Error vs Erasure Capacity (MDS Singleton Bound)
3. Burst vs Random Noise Tolerance: Symbol-level preservation
"""

import random
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from reed_solomon import rs_encode_msg, rs_correct_msg

OUTPUT_DIR = Path(__file__).parent / "demo_output"
OUTPUT_DIR.mkdir(exist_ok=True)


def run_sensitivity_analysis():
    print("Generating Figure 6: Sensitivity Analysis (Code Rate vs Error Rate)...")
    nsym_list = [8, 16, 32, 64]
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]
    rates = np.linspace(0.005, 0.20, 15)
    trials = 150
    rng = random.Random(42)

    plt.figure(figsize=(8, 5))

    for nsym, color in zip(nsym_list, colors):
        k = 255 - nsym
        block_size = 255
        theo_limit = (nsym // 2) / block_size
        success_rates = []

        for r in rates:
            ok = 0
            for _ in range(trials):
                data = [rng.randint(0, 255) for _ in range(k)]
                block = rs_encode_msg(data, nsym)
                n_errors = min(block_size, np.random.binomial(block_size, r))
                err_positions = rng.sample(range(block_size), n_errors)
                corrupted = list(block)
                for p in err_positions:
                    corrupted[p] ^= rng.randint(1, 255)
                try:
                    fixed = rs_correct_msg(corrupted, nsym)
                    if fixed == block:
                        ok += 1
                except ValueError:
                    pass
            success_rates.append(ok / trials * 100)

        plt.plot(rates * 100, success_rates, marker="o", label=f"RS(255,{k}) [nsym={nsym}]", color=color, lw=2)
        plt.axvline(theo_limit * 100, color=color, linestyle=":", alpha=0.7, label=f"Limit ({theo_limit*100:.1f}%)" if nsym in [8, 64] else "")

    plt.title("Sensitivity Analysis: Block Recovery vs. Error Rate Across Code Rates", fontsize=12, fontweight="bold")
    plt.xlabel("Channel Byte Error Rate (%)", fontsize=11)
    plt.ylabel("Block Recovery Success Rate (%)", fontsize=11)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(loc="lower left", fontsize=9)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "6_sensitivity_analysis.png", dpi=300)
    plt.close()
    print("Saved 6_sensitivity_analysis.png")


def run_mds_erasure_vs_error_tradeoff():
    print("Generating Figure 7: MDS Frontier (Errors vs Erasures Trade-off)...")
    nsym = 32
    block_size = 255
    k = block_size - nsym
    trials = 60
    rng = random.Random(99)

    # 2s + r <= nsym
    # We test combinations of s (errors) and r (erasures)
    s_vals = [0, 4, 8, 12, 16]
    r_max = nsym

    plt.figure(figsize=(8, 5))

    # Theoretical line
    s_line = np.linspace(0, nsym // 2, 100)
    r_line = nsym - 2 * s_line
    plt.plot(s_line, r_line, "r--", lw=2, label="Singleton MDS Bound: 2s + r = 32")
    plt.fill_between(s_line, 0, r_line, color="green", alpha=0.15, label="Correctable Region")

    # Empirical test points
    test_points = [
        (0, 32, "100% Corrected"),
        (4, 24, "100% Corrected"),
        (8, 16, "100% Corrected"),
        (12, 8, "100% Corrected"),
        (16, 0, "100% Corrected"),
        (8, 18, "Failed (exceeds bound)"),
        (12, 10, "Failed (exceeds bound)"),
    ]

    for s, r, note in test_points:
        ok = 0
        for _ in range(trials):
            data = [rng.randint(0, 255) for _ in range(k)]
            block = rs_encode_msg(data, nsym)

            # Pick r erasure positions and s error positions disjointly
            all_pos = rng.sample(range(block_size), s + r)
            erase_pos = all_pos[:r]
            err_pos = all_pos[r:]

            corrupted = list(block)
            # Erase: set to 0 and pass position
            for ep in erase_pos:
                corrupted[ep] = 0
            # Error: corrupt value, do not disclose position
            for ep in err_pos:
                corrupted[ep] ^= rng.randint(1, 255)

            try:
                fixed = rs_correct_msg(corrupted, nsym, erase_pos=erase_pos)
                if fixed == block:
                    ok += 1
            except (ValueError, Exception):
                pass

        rate = ok / trials
        marker = "o" if rate > 0.9 else "x"
        c = "blue" if rate > 0.9 else "darkred"
        plt.scatter(s, r, color=c, marker=marker, s=80, zorder=5)
        plt.annotate(f"s={s}, r={r}\n({rate*100:.0f}%)", (s + 0.3, r + 0.5), fontsize=8)

    plt.title("Empirical Verification of the MDS Trade-off: 2s + r ≤ n - k (nsym=32)", fontsize=12, fontweight="bold")
    plt.xlabel("Number of Unknown Errors (s)", fontsize=11)
    plt.ylabel("Number of Known Erasures (r)", fontsize=11)
    plt.xlim(-1, 18)
    plt.ylim(-1, 36)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "7_mds_bound_verification.png", dpi=300)
    plt.close()
    print("Saved 7_mds_bound_verification.png")


if __name__ == "__main__":
    run_sensitivity_analysis()
    run_mds_erasure_vs_error_tradeoff()
    print("All extended analysis figures generated in demo_output/ successfully.")

