"""Irreversibility detection across three independent datasets.

The entropy-production test (directed velocity vs a velocity-shuffled control) is
the pipeline's robust result. This summarizes it across three datasets from three
tissues, each measured with the identical setup (500 cells, 50 genes, 5 seeds,
100 permutations, proxy velocity):

  - pancreatic endocrinogenesis   (Bastidas-Ponce et al. 2019)
  - mouse gastrulation erythroid  (Pijuan-Sala et al. 2019)
  - human bone marrow             (Setty et al. 2019)

Each is reproduced by its own benchmark (pancreas_entropy_validation.py,
second_dataset_validation.py, and the bone-marrow run); this script records the
measured means and draws the figure.

    python benchmarks/generalization.py
"""
from __future__ import annotations

import json
from pathlib import Path

RESULTS = [
    {"dataset": "pancreas endocrine", "tissue": "pancreas", "directed_p": 0.010, "shuffled_p": 0.725},
    {"dataset": "gastrulation erythroid", "tissue": "gastrulation", "directed_p": 0.010, "shuffled_p": 0.347},
    {"dataset": "bone marrow", "tissue": "hematopoiesis", "directed_p": 0.010, "shuffled_p": 0.352},
]


def main():
    report = {
        "benchmark": "irreversibility_generalization",
        "setup": "500 cells, 50 genes, 5 seeds, 100 permutations, proxy velocity",
        "results": RESULTS,
        "conclusion": "Entropy production is significant (directed p=0.010, beating all permutations) on all three datasets, while the shuffled-velocity control is not. The irreversibility result generalizes across three tissues.",
    }
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/generalization.json").write_text(json.dumps(report, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    names = [r["dataset"] for r in RESULTS]
    directed = [r["directed_p"] for r in RESULTS]
    shuffled = [r["shuffled_p"] for r in RESULTS]
    x = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    ax.bar(x - 0.2, directed, 0.4, label="directed (real velocity)", color="#0f7d99")
    ax.bar(x + 0.2, shuffled, 0.4, label="shuffled control", color="#b0413e")
    ax.axhline(0.05, color="black", ls="--", lw=1)
    ax.text(len(names) - 0.5, 0.07, "p = 0.05", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['dataset']}\n({r['tissue']})" for r in RESULTS], fontsize=9)
    ax.set_ylabel("entropy-production permutation p-value")
    ax.set_ylim(0, 1)
    ax.set_title("Irreversibility detection generalizes across three tissues")
    ax.legend(frameon=False)
    for xi, d in zip(x, directed):
        ax.text(xi - 0.2, d + 0.02, f"{d:.3f}", ha="center", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    Path("figures").mkdir(exist_ok=True)
    fig.savefig("figures/generalization.png", dpi=150, bbox_inches="tight")
    print("wrote experiments/generalization.json and figures/generalization.png")
    for r in RESULTS:
        print(f"  {r['dataset']:24s} directed p={r['directed_p']:.3f}  shuffled p={r['shuffled_p']:.3f}")


if __name__ == "__main__":
    main()
