"""Plot the pancreas entropy-production validation.

Reads experiments/pancreas_entropy_validation.json and writes
figures/pancreas_entropy_validation.png.
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
    ap.add_argument("--data", default="experiments/pancreas_entropy_validation.json")
    ap.add_argument("--out", default="figures/pancreas_entropy_validation.png")
    args = ap.parse_args()

    d = json.loads(Path(args.data).read_text())
    obs, ctl = d["observed"], d["shuffled_velocity_control"]
    ep = [obs["entropy_production_rate"], ctl["entropy_production_rate"]]
    pv = [obs["permutation_p_value"], ctl["permutation_p_value"]]
    names = ["observed", "shuffled velocity\n(negative control)"]
    colors = ["#2f6f9f", "#b0413e"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.4))

    ax1.bar(names, ep, color=colors, width=0.55)
    ax1.set_yscale("log")
    ax1.set_ylabel("entropy production rate (arb. units)")
    ax1.set_title("Entropy production")
    ax1.grid(True, axis="y", alpha=0.3)

    bars = ax2.bar(names, pv, color=colors, width=0.55)
    ax2.axhline(0.05, color="black", ls="--", lw=1)
    ax2.text(1.02, 0.05, "p = 0.05", va="center", fontsize=9, transform=ax2.get_yaxis_transform())
    for bar, p in zip(bars, pv):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02, f"p={p:.3f}", ha="center", fontsize=10)
    ax2.set_ylim(0, max(1.0, max(pv) * 1.15))
    ax2.set_ylabel("permutation p-value")
    ax2.set_title("Significance vs shuffled-velocity null")
    ax2.grid(True, axis="y", alpha=0.3)

    fig.suptitle(
        f"Pancreas endocrine branch, proxy velocity ({d['n_cells']} cells, {d['n_genes']} genes)",
        y=1.02,
    )
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
