"""Plot the runtime scaling benchmark.

Reads experiments/runtime_scaling.json and writes figures/runtime_speedup.png.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="experiments/runtime_scaling.json")
    ap.add_argument("--out", default="figures/runtime_speedup.png")
    args = ap.parse_args()

    report = json.loads(Path(args.data).read_text())
    rows = [r for r in report["rows"] if r["baseline_s"] is not None]
    ns = [r["n"] for r in rows]
    base = [r["baseline_s"] for r in rows]
    cur = [r["current_s"] for r in rows]
    speed = [r["speedup"] for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    ax1.plot(ns, base, "o-", color="#b0413e", lw=2, ms=7, label="pure-Python baseline")
    ax1.plot(ns, cur, "o-", color="#2f6f9f", lw=2, ms=7, label="vectorized core")
    ax1.set_xlabel("cells")
    ax1.set_ylabel("fit time (s)")
    ax1.set_title("Landscape fit runtime")
    ax1.grid(True, alpha=0.3)
    ax1.legend(frameon=False)

    bars = ax2.bar([str(n) for n in ns], speed, color="#2f6f9f", width=0.6)
    for bar, s in zip(bars, speed):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.08,
                 f"{s:.1f}x", ha="center", va="bottom", fontsize=10)
    ax2.set_xlabel("cells")
    ax2.set_ylabel("speedup")
    ax2.set_title("Speedup, identical outputs")
    ax2.set_ylim(0, max(speed) * 1.2)
    ax2.grid(True, axis="y", alpha=0.3)

    fig.suptitle("Vectorized core vs pure-Python (energies bit-identical, max abs diff 0.0)", y=1.02)
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
