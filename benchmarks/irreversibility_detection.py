"""Predictive benchmark: detecting irreversibility (the task this pipeline is for).

The ordering benchmark asked "can you recover developmental stage", and simple
baselines win that. This asks the question the pipeline is actually built to
answer: "is the process irreversible, or could it run backwards?"

For each of N random subsamples of the pancreas endocrine branch we form a
directed dataset (real proxy velocity) and a matched equilibrium control (the
velocity vectors permuted across cells, which destroys the direction while
keeping every geometric feature). Each method scores each dataset, and we report
how well the scores separate directed from control with AUROC over the 2N
datasets.

The baselines are the same ones that won the ordering benchmark (first principal
component, distance from the progenitor). They never look at velocity direction,
so they score a directed dataset and its shuffled control identically and land at
chance. That is the point: expression-based trajectory methods cannot see
irreversibility. The entropy-production test can.

    python benchmarks/irreversibility_detection.py
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from thermodynamic_waddington.entropy_production import EntropyProductionConfig, estimate_entropy_production
from thermodynamic_waddington.graph import build_knn, estimate_local_diffusion, local_density
from thermodynamic_waddington.linalg import pca

from pancreas_entropy_validation import BRANCH, load_branch


def _core(expr, vel, seed, k=20, dims=6):
    pts, _, _ = pca(expr.tolist(), dims, seed)
    graph = build_knn(pts, k)
    dens = local_density(pts, graph, 0.65)
    diff = estimate_local_diffusion(graph, vel.tolist(), 1e-4)
    return pts, graph, dens, diff


def tw_score(expr, vel, seed, perms):
    pts, graph, dens, diff = _core(expr, vel, seed)
    cfg = EntropyProductionConfig(temperature=1.0, velocity_scale=1.0, bootstrap_replicates=0,
                                  permutation_replicates=perms, seed=seed + 11)
    rep = estimate_entropy_production(pts, vel.tolist(), dens, diff, graph, cfg).to_dict()
    p = rep.get("permutation_p_value", 1.0)
    return -math.log10(p + 1e-4)  # higher = more irreversible


def velocity_coherence(expr, vel, seed):
    pts, graph, _, _ = _core(expr, vel, seed)
    v = np.asarray(vel, dtype=float)
    norms = np.linalg.norm(v, axis=1) + 1e-12
    vn = v / norms[:, None]
    vals = []
    for i, nbrs in enumerate(graph.neighbors):
        if nbrs:
            vals.append(float(np.mean(vn[nbrs] @ vn[i])))
    return float(np.mean(vals)) if vals else 0.0


def blind_pc1(expr, vel, seed):
    x = np.asarray(expr, dtype=float)
    xc = x - x.mean(axis=0)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    return float(np.std(xc @ vt[0]))


def blind_dist(expr, vel, seed):
    x = np.asarray(expr, dtype=float)
    xc = x - x.mean(axis=0)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    pca6 = xc @ vt[:6].T
    return float(np.std(np.linalg.norm(pca6 - pca6.mean(axis=0), axis=1)))


METHODS = {
    "TW entropy production": (tw_score, "tw"),
    "velocity coherence": (velocity_coherence, "velocity"),
    "PC1 (expression only)": (blind_pc1, "blind"),
    "distance from progenitor": (blind_dist, "blind"),
}


def run(path, seeds, cells, genes, perms):
    labels_bin = []
    scores = {name: [] for name in METHODS}
    directed_p, control_p = [], []
    for s in range(seeds):
        expr, vel, labels = load_branch(path, BRANCH, cells, genes, s)
        rng = np.random.default_rng(1000 + s)
        vel_shuf = vel[rng.permutation(len(labels))]
        for tag, v in [("directed", vel), ("control", vel_shuf)]:
            labels_bin.append(1 if tag == "directed" else 0)
            for name, (fn, _) in METHODS.items():
                if name == "TW entropy production":
                    val = fn(expr, v, s, perms)
                else:
                    val = fn(expr, v, s)
                scores[name].append(val)
        print(f"  seed {s}: scored directed and control")
    rows = []
    for name, (_, kind) in METHODS.items():
        auc = float(roc_auc_score(labels_bin, scores[name]))
        rows.append({"method": name, "kind": kind, "auroc": round(auc, 3)})
    rows.sort(key=lambda r: r["auroc"], reverse=True)
    return {
        "benchmark": "irreversibility_detection_pancreas",
        "task": "separate directed differentiation from a velocity-shuffled equilibrium control",
        "metric": "AUROC over 2N datasets",
        "n_pairs": seeds,
        "n_cells": cells,
        "rows": rows,
        "note": "Expression-only baselines never use velocity direction, so they score a dataset and its shuffled control identically and sit at 0.5. This is the task the pipeline is built for.",
    }


def plot(report, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = report["rows"]
    names = [r["method"] for r in rows]
    vals = [r["auroc"] for r in rows]
    palette = {"tw": "#0f7d99", "velocity": "#6a4c93", "blind": "#b0413e"}
    colors = [palette[r["kind"]] for r in rows]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    y = list(range(len(names)))
    ax.barh(y, vals, color=colors)
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.axvline(0.5, color="black", ls="--", lw=1)
    ax.text(0.5, -0.6, "chance", ha="center", fontsize=9)
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("AUROC: directed vs equilibrium control")
    for i, v in enumerate(vals):
        ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=10)
    ax.set_title("Detecting irreversibility: the task this pipeline is for")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=palette["tw"], label="entropy production (this pipeline)"),
                       Patch(color=palette["velocity"], label="velocity heuristic"),
                       Patch(color=palette["blind"], label="expression-only baseline")],
              frameon=False, loc="lower right")
    ax.grid(True, axis="x", alpha=0.3)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--path", default="data/real/endocrinogenesis_day15.h5ad")
    ap.add_argument("--seeds", type=int, default=12)
    ap.add_argument("--cells", type=int, default=300)
    ap.add_argument("--genes", type=int, default=50)
    ap.add_argument("--permutations", type=int, default=60)
    ap.add_argument("--out", default="experiments/irreversibility_detection.json")
    ap.add_argument("--figure", default="figures/irreversibility_detection.png")
    args = ap.parse_args()
    report = run(args.path, args.seeds, args.cells, args.genes, args.permutations)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(f"task: {report['task']}")
    for r in report["rows"]:
        print(f"  {r['method']:28s} AUROC={r['auroc']:.3f}  [{r['kind']}]")
    plot(report, args.figure)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
