"""Commitment free-energy barrier along the committor, with uncertainty.

For each random subsample of the pancreas endocrine branch, fits the landscape,
computes the committor, and takes the potential of mean force F(q) = -kT ln P(q)
along it. Reports the commitment barrier (transition-state height in kT) as a
mean with a confidence interval, and plots the per-seed profiles with the mean.

    python benchmarks/commitment_barrier.py --seeds 8
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from thermodynamic_waddington import (
    FitConfig,
    committor_free_energy_profile,
    developmental_coordinate,
    fit_landscape,
)

from pancreas_entropy_validation import BRANCH, load_branch


def run(seeds, cells, genes, grid):
    curves = []
    barriers = []
    barrier_qs = []
    for s in range(seeds):
        expr, vel, labels = load_branch("data/real/endocrinogenesis_day15.h5ad", BRANCH, cells, genes, s)
        cfg = FitConfig(neighbors=20, dimensions=6, bootstrap_replicates=2, enable_cycle_decomposition=False,
                        entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2, seed=s)
        fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=labels, metadata={"velocity_status": "derived_proxy"})
        q = developmental_coordinate(fit, ["Ductal"], ["Beta"])
        prof = committor_free_energy_profile(q, temperature=1.0, grid=grid)
        curves.append(prof.free_energy_kt)
        barriers.append(prof.barrier_kt)
        barrier_qs.append(prof.barrier_q)
        print(f"  seed {s}: barrier {prof.barrier_kt:.2f} kT at q={prof.barrier_q:.2f}")
    b = np.asarray(barriers)
    n = b.size
    sem = float(b.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    grid_x = np.linspace(0.0, 1.0, grid).tolist()
    return {
        "benchmark": "commitment_barrier_pancreas",
        "seeds": seeds,
        "n_cells": cells,
        "barrier_kt_mean": round(float(b.mean()), 3),
        "barrier_kt_ci95": [round(float(b.mean() - 1.96 * sem), 3), round(float(b.mean() + 1.96 * sem), 3)],
        "barrier_q_mean": round(float(np.mean(barrier_qs)), 3),
        "grid": grid_x,
        "curves": curves,
        "note": "Free-energy barrier a cell must cross to commit, from the potential of mean force along the committor. kT is the density-based pseudopotential convention.",
    }


def plot(report, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    grid = np.asarray(report["grid"])
    curves = np.asarray(report["curves"])
    mean = curves.mean(axis=0)
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    for c in curves:
        ax.plot(grid, c, color="#0f7d99", alpha=0.18, lw=1)
    ax.plot(grid, mean, color="#0f7d99", lw=2.6, label="mean PMF")
    bq = report["barrier_q_mean"]
    bk = report["barrier_kt_mean"]
    ax.axvline(bq, color="#b0413e", ls="--", lw=1)
    ax.annotate(f"barrier {bk:.1f} kT", xy=(bq, mean.max()), xytext=(bq + 0.03, mean.max() * 0.9),
                color="#b0413e", fontsize=11)
    ax.set_xlabel("committor q  (0 = progenitor, 1 = terminal fate)")
    ax.set_ylabel("free energy along reaction coordinate (kT)")
    ci = report["barrier_kt_ci95"]
    ax.set_title(f"Commitment barrier: {bk:.1f} kT  [{ci[0]:.1f}, {ci[1]:.1f}]  ({report['seeds']} seeds)")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--cells", type=int, default=400)
    ap.add_argument("--genes", type=int, default=50)
    ap.add_argument("--grid", type=int, default=60)
    ap.add_argument("--out", default="experiments/commitment_barrier.json")
    ap.add_argument("--figure", default="figures/commitment_barrier.png")
    args = ap.parse_args()
    report = run(args.seeds, args.cells, args.genes, args.grid)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    ci = report["barrier_kt_ci95"]
    print(f"commitment barrier = {report['barrier_kt_mean']} kT  95% CI [{ci[0]}, {ci[1]}]  at q={report['barrier_q_mean']}")
    plot(report, args.figure)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
