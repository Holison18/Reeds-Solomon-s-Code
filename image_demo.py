"""
Visual demo of the Reed-Solomon codec in reed_solomon.py: encode an image's
raw bytes in RS(255, 255-nsym) blocks, corrupt the encoded stream (random
byte noise plus a contiguous "scratch"), then reconstruct the image two
ways -- naively (ignoring the parity, i.e. no error correction) and via
rs_correct_msg -- so the difference is visible side by side.

Usage:
    python image_demo.py --image photo.jpg
"""

import argparse
import random
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from reed_solomon import rs_encode_msg, rs_correct_msg

OUTPUT_DIR = Path(__file__).parent / "demo_output"


def load_image(path: str) -> Image.Image:
    return Image.open(path).convert("RGB")


# Stream <-> blocks helpers

def encode_stream(data: bytes, k: int, nsym: int):
    pad_len = (-len(data)) % k
    padded = list(data) + [0] * pad_len
    blocks = [padded[i:i + k] for i in range(0, len(padded), k)]
    return [rs_encode_msg(b, nsym) for b in blocks]


def corrupt_blocks(blocks, scatter_rate: float, burst_fraction: float, rng: random.Random):
    block_size = len(blocks[0])
    flat = bytearray(b for block in blocks for b in block)
    n = len(flat)

    n_scatter = int(n * scatter_rate)
    for p in rng.sample(range(n), n_scatter):
        orig = flat[p]
        new = rng.randint(0, 255)
        while new == orig:
            new = rng.randint(0, 255)
        flat[p] = new
    n_scattered = n_scatter

    n_burst = 0
    if burst_fraction > 0:
        burst_len = int(n * burst_fraction)
        start = rng.randint(0, max(0, n - burst_len))
        for i in range(start, start + burst_len):
            flat[i] = rng.randint(0, 255)
        n_burst = burst_len

    corrupted = [list(flat[i:i + block_size]) for i in range(0, n, block_size)]
    return corrupted, n_scattered, n_burst


def decode_naive(blocks, k: int) -> bytes:
    return bytes(b for block in blocks for b in block[:k])


def decode_rs(blocks, k: int, nsym: int):
    out = bytearray()
    stats = {"total": len(blocks), "corrected": 0, "failed": 0, "clean": 0}
    for block in blocks:
        try:
            fixed = rs_correct_msg(block, nsym)
            if fixed == block:
                stats["clean"] += 1
            else:
                stats["corrected"] += 1
            out.extend(fixed[:k])
        except ValueError:
            stats["failed"] += 1
            out.extend(block[:k])
    return bytes(out), stats


def bytes_to_image(data: bytes, size, length: int) -> Image.Image:
    trimmed = data[:length]
    arr = np.frombuffer(trimmed, dtype=np.uint8)
    return Image.frombytes("RGB", size, arr.tobytes())


# Sweep: success rate vs corruption rate (theoretical, block-level)


def sweep_success_rate(k: int, nsym: int, rates, trials: int, rng: random.Random):
    block_size = k + nsym
    max_correctable = nsym // 2
    successes = []
    for rate in rates:
        ok = 0
        for _ in range(trials):
            data = [rng.randint(0, 255) for _ in range(k)]
            block = rs_encode_msg(data, nsym)
            n_errors = min(block_size, np.random.binomial(block_size, rate))
            positions = rng.sample(range(block_size), n_errors)
            corrupted = list(block)
            for p in positions:
                corrupted[p] ^= rng.randint(1, 255)
            try:
                fixed = rs_correct_msg(corrupted, nsym)
                if fixed == block:
                    ok += 1
            except ValueError:
                pass
        successes.append(ok / trials)
    return successes, max_correctable / block_size


# Main demo


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", type=str, required=True, help="path to the image file to corrupt and recover")
    parser.add_argument("--nsym", type=int, default=32, help="parity bytes per 255-byte block")
    parser.add_argument("--scatter-rate", type=float, default=0.03, help="fraction of bytes randomly flipped")
    parser.add_argument("--burst-fraction", type=float, default=0.02, help="fraction of the stream corrupted as one contiguous scratch")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--sweep-trials", type=int, default=300)
    parser.add_argument("--sweep-points", type=int, default=12)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    np.random.seed(args.seed)
    OUTPUT_DIR.mkdir(exist_ok=True)

    img = load_image(args.image)
    size = img.size
    data = img.tobytes()
    k = 255 - args.nsym

    print(f"image: {size[0]}x{size[1]} RGB, {len(data)} bytes")
    print(f"RS block: k={k} data + nsym={args.nsym} parity = {k + args.nsym} bytes/block "
          f"(corrects up to {args.nsym // 2} bad bytes/block)")

    blocks = encode_stream(data, k, args.nsym)
    corrupted, n_scattered, n_burst = corrupt_blocks(blocks, args.scatter_rate, args.burst_fraction, rng)
    total_bytes = len(blocks) * (k + args.nsym)
    print(f"corrupted {n_scattered} scattered bytes + {n_burst} bytes in a contiguous scratch "
          f"out of {total_bytes} total ({(n_scattered + n_burst) / total_bytes:.1%})")

    naive_bytes = decode_naive(corrupted, k)
    rs_bytes, stats = decode_rs(corrupted, k, args.nsym)
    print(f"blocks: {stats['total']} total, {stats['clean']} clean, "
          f"{stats['corrected']} corrected, {stats['failed']} unrecoverable")

    original_img = img
    naive_img = bytes_to_image(naive_bytes, size, len(data))
    recovered_img = bytes_to_image(rs_bytes, size, len(data))

    if stats["failed"] == 0:
        assert rs_bytes[:len(data)] == data, "recovered bytes should exactly match the original"
        print("recovered image is byte-for-byte identical to the original")
    else:
        print(f"{stats['failed']} block(s) exceeded correction capacity -- recovered image will show residual damage there")

    original_img.save(OUTPUT_DIR / "1_original.png")
    naive_img.save(OUTPUT_DIR / "2_corrupted_no_correction.png")
    recovered_img.save(OUTPUT_DIR / "3_recovered_with_rs.png")

    # Side-by-side collage for slides
    fig, axes = plt.subplots(1, 3, figsize=(12, 5))
    for ax, im, title in zip(
        axes,
        [original_img, naive_img, recovered_img],
        ["Original", "Corrupted\n(no error correction)", "Recovered\n(Reed-Solomon corrected)"],
    ):
        ax.imshow(im)
        ax.set_title(title, fontsize=12)
        ax.axis("off")
    corruption_pct = (n_scattered + n_burst) / total_bytes
    fig.suptitle(
        f"RS({k + args.nsym},{k}) per block  |  {corruption_pct:.1%} of bytes corrupted  |  "
        f"{stats['corrected']}/{stats['total']} blocks corrected, {stats['failed']} failed",
        fontsize=11,
        y=0.99,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(OUTPUT_DIR / "4_comparison.png", dpi=150)
    plt.close(fig)

    # Sweep chart: success rate vs corruption rate
    max_theoretical = (args.nsym // 2) / (k + args.nsym)
    rates = np.linspace(0, min(0.5, max_theoretical * 3), args.sweep_points)
    successes, capacity = sweep_success_rate(k, args.nsym, rates, args.sweep_trials, rng)

    fig2, ax2 = plt.subplots(figsize=(7, 4.5))
    ax2.plot(rates * 100, np.array(successes) * 100, marker="o")
    ax2.axvline(capacity * 100, color="red", linestyle="--", label=f"theoretical limit ({capacity:.1%} bytes/block)")
    ax2.set_xlabel("Random byte corruption rate (%)")
    ax2.set_ylabel("Block recovery success rate (%)")
    ax2.set_title(f"RS({k + args.nsym},{k}) recovery success vs. corruption rate")
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(OUTPUT_DIR / "5_success_rate_curve.png", dpi=150)
    plt.close(fig2)

    print(f"\nsaved images and charts to {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
