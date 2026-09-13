"""Physical free-energy calibration on pancreatic endocrinogenesis.

Fits the landscape on the endocrine branch, expresses it in kT via the
density-anchored Boltzmann inversion, and reports:
  - the landscape depth in kT (the Waddington pseudopotential scale)
  - basin depth per cell type in kT
  - the agreement (R^2) between the path-work landscape and the density
    landscape, which is low precisely because the path-work landscape carries
    the non-equilibrium structure the density landscape cannot see.

Writes a summary JSON and a two-panel figure.

    python benchmarks/free_energy_calibration.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from thermodynamic_waddington.calibration import basin_barriers_kt, boltzmann_free_energy, calibrate
from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.model import fit_landscape

from pancreas_entropy_validation import BRANCH, load_branch


def run(path, cells, genes, seed, permutations):
    expr, vel, labels = load_branch(path, BRANCH, cells, genes, seed)
    cfg = FitConfig(
        neighbors=20, dimensions=6, bootstrap_replicates=6,
        entropy_production_bootstrap_replicates=8,
        entropy_production_permutation_replicates=permutations, seed=seed,
    )
    fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=labels,
                        metadata={"velocity_status": "derived_proxy"})
    rep = calibrate(fit.energies, fit.embedding)
    fb = np.asarray(rep.boltzmann_energies_kt)
    bb = basin_barriers_kt(rep.boltzmann_energies_kt, fit.labels)
    ep = fit.diagnostics.get("entropy_production", {})
    summary = {
        "benchmark": "free_energy_calibration_pancreas",
        "dataset": "scVelo endocrinogenesis_day15 (endocrine branch)",
        "n_cells": rep.n_cells,
        "landscape_range_kt": rep.boltzmann_energy_range_kt,
        "kt_per_work_unit": rep.kt_per_work_unit,
        "path_vs_density_r2": rep.r_squared,
        "entropy_production_p": ep.get("permutation_p_value"),
        "basin_depth_kt": {k: bb[k]["mean_kt"] for k in BRANCH if k in bb},
        "assumptions": rep.assumptions,
        "note": "Low path-vs-density R^2 reflects non-equilibrium structure the density landscape misses; entropy production is significant on the same data.",
    }
    return summary, fit, rep, fb, bb


def plot(summary, fit, rep, fb, bb, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    order = [k for k in BRANCH if k in bb]
    depth = [bb[k]["mean_kt"] for k in order]
    ax1.bar(range(len(order)), depth, color="#0f7d99", width=0.62)
    ax1.set_xticks(range(len(order)))
    ax1.set_xticklabels([k.replace(" ", "\n") for k in order], fontsize=9)
    ax1.set_ylabel("density landscape depth (kT)")
    ax1.set_title(f"Waddington pseudopotential, range {summary['landscape_range_kt']:.1f} kT")
    ax1.grid(True, axis="y", alpha=0.3)

    fp = np.asarray(fit.energies)
    fp = (fp - fp.min())
    ax2.scatter(fp, fb, s=12, alpha=0.5, color="#6a4c93", edgecolors="none")
    ax2.set_xlabel("path-work free energy (arb. units)")
    ax2.set_ylabel("density landscape (kT)")
    ax2.set_title(f"Path-work vs density,  R² = {rep.r_squared:.2f}")
    ax2.grid(True, alpha=0.3)

    pval = summary["entropy_production_p"]
    fig.suptitle(
        f"Free-energy calibration, pancreas endocrine branch ({summary['n_cells']} cells), EP p={pval:.4f}",
        y=1.02,
    )
    fig.tight_layout()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--path", default="data/real/endocrinogenesis_day15.h5ad")
    ap.add_argument("--cells", type=int, default=500)
    ap.add_argument("--genes", type=int, default=50)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--permutations", type=int, default=100)
    ap.add_argument("--out", default="experiments/free_energy_calibration.json")
    ap.add_argument("--figure", default="figures/free_energy_calibration.png")
    args = ap.parse_args()
    summary, fit, rep, fb, bb = run(args.path, args.cells, args.genes, args.seed, args.permutations)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(summary, indent=2))
    print(f"landscape range = {summary['landscape_range_kt']:.2f} kT")
    print(f"path-vs-density R2 = {summary['path_vs_density_r2']:.3f}   EP p = {summary['entropy_production_p']}")
    for k, v in summary["basin_depth_kt"].items():
        print(f"  {k:16s} {v:5.2f} kT")
    plot(summary, fit, rep, fb, bb, args.figure)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
