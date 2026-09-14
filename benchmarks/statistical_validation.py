"""Statistical validation of the ordering result across seeds.

A single benchmark run is not evidence. This repeats the developmental-ordering
comparison over many random subsamples and reports, for each method, the mean
and 95% confidence interval of the Spearman correlation with developmental
stage, the fraction of seeds in which the committor beats it, and a paired
Wilcoxon test of the committor against the strongest baseline.

    python benchmarks/statistical_validation.py --seeds 10
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr, wilcoxon

from thermodynamic_waddington import FitConfig, developmental_coordinate, fit_landscape

from pancreas_entropy_validation import BRANCH, load_branch


def _pca(expr, dims=6):
    x = np.asarray(expr, dtype=float)
    xc = x - x.mean(axis=0)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    return xc, vt


def scores_for_seed(seed, cells, genes):
    expr, vel, labels = load_branch("data/real/endocrinogenesis_day15.h5ad", BRANCH, cells, genes, seed)
    stage = np.array([BRANCH.index(l) for l in labels], dtype=float)
    lab = np.array(labels)
    cfg = FitConfig(neighbors=20, dimensions=6, bootstrap_replicates=2, enable_cycle_decomposition=False,
                    entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2, seed=seed)
    fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=labels, metadata={"velocity_status": "derived_proxy"})
    q = np.array(developmental_coordinate(fit, ["Ductal"], ["Beta"]))
    xc, vt = _pca(expr)
    pca = xc @ vt[:6].T
    axis = pca[lab == "Beta"].mean(axis=0) - pca[lab == "Ductal"].mean(axis=0)
    return {
        "committor": abs(spearmanr(q, stage).correlation),
        "PC1": abs(spearmanr(xc @ vt[0], stage).correlation),
        "Ductal->Beta axis": abs(spearmanr(pca @ axis, stage).correlation),
        "TW free energy": abs(spearmanr(np.asarray(fit.energies), stage).correlation),
    }


def summarize(values):
    a = np.asarray(values, dtype=float)
    n = a.size
    mean = float(a.mean())
    sem = float(a.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    return {"mean": round(mean, 4), "std": round(float(a.std(ddof=1)) if n > 1 else 0.0, 4),
            "ci95_low": round(mean - 1.96 * sem, 4), "ci95_high": round(mean + 1.96 * sem, 4)}


def run(seeds, cells, genes):
    per = {m: [] for m in ["committor", "PC1", "Ductal->Beta axis", "TW free energy"]}
    for s in range(seeds):
        sc = scores_for_seed(s, cells, genes)
        for m, v in sc.items():
            per[m].append(v)
        print(f"  seed {s}: committor={sc['committor']:.3f}  PC1={sc['PC1']:.3f}")
    stats = {m: summarize(v) for m, v in per.items()}
    comm = np.asarray(per["committor"])
    baselines = ["PC1", "Ductal->Beta axis", "TW free energy"]
    best_base = max(baselines, key=lambda m: np.mean(per[m]))
    wins = int((comm > np.asarray(per[best_base])).sum())
    try:
        p = float(wilcoxon(comm, np.asarray(per[best_base])).pvalue)
    except ValueError:
        p = float("nan")
    return {
        "benchmark": "statistical_validation_ordering",
        "seeds": seeds,
        "n_cells": cells,
        "stats": stats,
        "committor_vs_best_baseline": {
            "best_baseline": best_base,
            "committor_wins": wins,
            "of_seeds": seeds,
            "paired_wilcoxon_p": p,
        },
        "raw": per,
    }


def plot(report, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = ["committor", "Ductal->Beta axis", "PC1", "TW free energy"]
    order = [m for m in order if m in report["stats"]]
    means = [report["stats"][m]["mean"] for m in order]
    errs = [[report["stats"][m]["mean"] - report["stats"][m]["ci95_low"] for m in order],
            [report["stats"][m]["ci95_high"] - report["stats"][m]["mean"] for m in order]]
    colors = ["#0f7d99" if m in ("committor", "TW free energy") else "#b0413e" for m in order]
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    y = list(range(len(order)))
    ax.barh(y, means, xerr=errs, color=colors, capsize=4)
    ax.set_yticks(y); ax.set_yticklabels(order); ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Spearman with developmental stage (mean, 95% CI)")
    for i, m in enumerate(order):
        ax.text(report["stats"][m]["mean"] + 0.02, i, f"{report['stats'][m]['mean']:.2f}", va="center", fontsize=10)
    cw = report["committor_vs_best_baseline"]
    ax.set_title(f"Ordering across {report['seeds']} seeds: committor wins {cw['committor_wins']}/{cw['of_seeds']} (p={cw['paired_wilcoxon_p']:.3f})")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--cells", type=int, default=400)
    ap.add_argument("--genes", type=int, default=50)
    ap.add_argument("--out", default="experiments/statistical_validation.json")
    ap.add_argument("--figure", default="figures/statistical_validation.png")
    args = ap.parse_args()
    report = run(args.seeds, args.cells, args.genes)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    for m, st in report["stats"].items():
        print(f"  {m:20s} {st['mean']:.3f}  [{st['ci95_low']:.3f}, {st['ci95_high']:.3f}]")
    cw = report["committor_vs_best_baseline"]
    print(f"committor beats {cw['best_baseline']} in {cw['committor_wins']}/{cw['of_seeds']} seeds, paired p={cw['paired_wilcoxon_p']:.4g}")
    plot(report, args.figure)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
