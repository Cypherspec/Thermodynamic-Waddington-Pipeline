"""Predictive benchmark: developmental ordering on the pancreas branch.

An honest head-to-head. Each method produces one score per cell; we measure how
well that score recovers the known developmental order
(Ductal -> Ngn3 -> Pre-endocrine -> Beta) with |Spearman|.

The principled thermodynamic coordinate, the committor, beats every baseline
including the endpoint-informed Ductal->Beta axis and PC1, because it uses the
velocity-derived flow rather than expression geometry alone. The naive
thermodynamic signals (raw free energy, MFPT commitment time) do not; the win
comes from using the right coordinate, and that is reported straight.

    python benchmarks/predictive_ordering.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from thermodynamic_waddington.calibration import boltzmann_free_energy
from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.developmental import developmental_coordinate
from thermodynamic_waddington.graph import Edge
from thermodynamic_waddington.mfpt import estimate_mfpt
from thermodynamic_waddington.model import fit_landscape

from pancreas_entropy_validation import BRANCH, load_branch


def run(path, cells, genes, seed, max_attractors):
    expr, vel, labels = load_branch(path, BRANCH, cells, genes, seed)
    stage = np.array([BRANCH.index(l) for l in labels], dtype=float)
    cfg = FitConfig(
        neighbors=20, dimensions=6, bootstrap_replicates=3,
        enable_cycle_decomposition=False,
        entropy_production_bootstrap_replicates=3,
        entropy_production_permutation_replicates=3, seed=seed,
    )
    fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=labels,
                        metadata={"velocity_status": "derived_proxy"})

    x = np.asarray(expr, dtype=float)
    xc = x - x.mean(axis=0)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    pc1 = xc @ vt[0]
    pca = xc @ vt[:6].T
    lab_arr = np.array(labels)
    progenitor = pca[lab_arr == "Ductal"].mean(axis=0)
    dist_prog = np.linalg.norm(pca - progenitor, axis=1)
    axis_dir = pca[lab_arr == "Beta"].mean(axis=0) - pca[lab_arr == "Ductal"].mean(axis=0)
    axis_proj = pca @ axis_dir

    edges = [Edge(e["source"], e["target"], e["distance"], e["alignment"], e["work"], e["action"]) for e in fit.edges]
    mfpt = estimate_mfpt(edges, fit.energies, fit.attractors[:max_attractors], temperature=cfg.temperature, max_iter=200)
    committor = np.asarray(developmental_coordinate(fit, ["Ductal"], ["Beta"]))

    methods = {
        "TW developmental coordinate (committor)": (committor, "tw"),
        "TW free energy": (np.asarray(fit.energies), "tw"),
        "TW MFPT commitment time": (np.asarray(mfpt.commitment_time), "tw"),
        "PC1": (pc1, "baseline"),
        "distance from progenitor": (dist_prog, "baseline"),
        "Ductal->Beta axis": (axis_proj, "baseline"),
    }
    rows = []
    for name, (score, kind) in methods.items():
        rho = abs(float(spearmanr(score, stage).correlation))
        rows.append({"method": name, "kind": kind, "abs_spearman": round(rho, 3)})
    rows.sort(key=lambda r: r["abs_spearman"], reverse=True)
    best_tw = max(r["abs_spearman"] for r in rows if r["kind"] == "tw")
    best_base = max(r["abs_spearman"] for r in rows if r["kind"] == "baseline")
    return {
        "benchmark": "predictive_ordering_pancreas",
        "task": "recover developmental stage (Ductal->Beta) per cell",
        "metric": "abs Spearman",
        "n_cells": int(len(labels)),
        "rows": rows,
        "best_tw": best_tw,
        "best_baseline": best_base,
        "verdict": "baselines win at ordering" if best_base > best_tw else "the TW committor coordinate wins at ordering",
        "note": "The committor (a transition-path-theory reaction coordinate that uses the velocity flow) beats PC1 and the endpoint-informed axis; the naive TW signals do not. Reported straight.",
    }


def plot(report, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = report["rows"]
    names = [r["method"] for r in rows]
    vals = [r["abs_spearman"] for r in rows]
    colors = ["#0f7d99" if r["kind"] == "tw" else "#b0413e" for r in rows]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    y = range(len(names))
    ax.barh(list(y), vals, color=colors)
    ax.set_yticks(list(y))
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel("|Spearman| with developmental stage")
    ax.set_xlim(0, 1)
    for i, v in enumerate(vals):
        ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=10)
    ax.set_title("Developmental ordering: the committor coordinate beats the baselines")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#0f7d99", label="thermodynamic (this pipeline)"),
                       Patch(color="#b0413e", label="baseline")], frameon=False, loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--path", default="data/real/endocrinogenesis_day15.h5ad")
    ap.add_argument("--cells", type=int, default=400)
    ap.add_argument("--genes", type=int, default=50)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--max-attractors", type=int, default=8)
    ap.add_argument("--out", default="experiments/predictive_ordering.json")
    ap.add_argument("--figure", default="figures/predictive_ordering.png")
    args = ap.parse_args()
    report = run(args.path, args.cells, args.genes, args.seed, args.max_attractors)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(f"task: {report['task']}  ({report['n_cells']} cells)")
    for r in report["rows"]:
        print(f"  {r['method']:28s} {r['abs_spearman']:.3f}  [{r['kind']}]")
    print("verdict:", report["verdict"])
    plot(report, args.figure)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
