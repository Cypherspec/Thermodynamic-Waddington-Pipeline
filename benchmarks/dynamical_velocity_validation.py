"""Headline results on real scVelo dynamical velocity, not the proxy.

Computes scVelo's dynamical-model velocity on the pancreas endocrine branch and
re-runs the two headline results on it: the entropy-production significance
(directed vs a velocity-shuffled control) and the committor ordering. The proxy
version is computed on the same cells and genes for a side-by-side.

Bounded on purpose (1000 genes) so the EM fit stays to a few minutes.

    python benchmarks/dynamical_velocity_validation.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from thermodynamic_waddington import FitConfig, developmental_coordinate, fit_landscape
from thermodynamic_waddington.entropy_production import EntropyProductionConfig, estimate_entropy_production
from thermodynamic_waddington.graph import build_knn, estimate_local_diffusion, local_density
from thermodynamic_waddington.linalg import pca

BRANCH = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta"]
PATH = "data/real/endocrinogenesis_day15.h5ad"


def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=float)


def compute_dynamical():
    import anndata as ad
    import scanpy as sc
    import scvelo as scv
    a = ad.read_h5ad(PATH)
    a = a[np.isin(a.obs["clusters"].astype(str).to_numpy(), BRANCH)].copy()
    scv.pp.filter_genes(a, min_shared_counts=20)
    scv.pp.normalize_per_cell(a)
    sc.pp.log1p(a)
    xv = dense(a.X)
    top = np.argsort(xv.var(axis=0))[::-1][:1000]
    a = a[:, top].copy()
    scv.pp.moments(a, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(a, n_jobs=1)
    scv.tl.velocity(a, mode="dynamical")
    vg = a.var["velocity_genes"].to_numpy() if "velocity_genes" in a.var else np.ones(a.n_vars, bool)
    expr = dense(a.layers["Ms"])[:, vg]
    vel = np.nan_to_num(dense(a.layers["velocity"])[:, vg])
    spliced = dense(a.layers["spliced"])[:, vg]
    unspliced = dense(a.layers["unspliced"])[:, vg]
    proxy = unspliced - spliced
    labels = a.obs["clusters"].astype(str).to_numpy()
    return expr, vel, proxy, labels


def ep_p(expr, vel, seed, perms=60):
    pts, _, _ = pca(expr.tolist(), 6, seed)
    g = build_knn(pts, 20)
    dens = local_density(pts, g, 0.65)
    diff = estimate_local_diffusion(g, vel.tolist(), 1e-4)
    cfg = EntropyProductionConfig(temperature=1.0, velocity_scale=1.0, bootstrap_replicates=0,
                                  permutation_replicates=perms, seed=seed + 11)
    return estimate_entropy_production(pts, vel.tolist(), dens, diff, g, cfg).to_dict().get("permutation_p_value")


def committor_rho(expr, vel, labels, seed):
    stage = np.array([BRANCH.index(l) for l in labels], dtype=float)
    cfg = FitConfig(neighbors=20, dimensions=6, bootstrap_replicates=2, enable_cycle_decomposition=False,
                    entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2, seed=seed)
    fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=labels.tolist(),
                        metadata={"velocity_status": "observed"})
    q = np.array(developmental_coordinate(fit, ["Ductal"], ["Beta"]))
    return abs(float(spearmanr(q, stage).correlation))


def evaluate(expr, vel, labels, seed=17, cap=500):
    if expr.shape[0] > cap:
        rng = np.random.default_rng(seed)
        idx = rng.choice(expr.shape[0], size=cap, replace=False)
        expr, vel, labels = expr[idx], vel[idx], labels[idx]
    rng = np.random.default_rng(1000 + seed)
    vel_shuf = vel[rng.permutation(len(labels))]
    return {
        "ep_directed_p": ep_p(expr, vel, seed),
        "ep_shuffled_p": ep_p(expr, vel_shuf, seed),
        "committor_ordering_rho": committor_rho(expr, vel, labels, seed),
        "n_cells": int(len(labels)),
    }


def main():
    print("computing scVelo dynamical velocity (this is the slow step)...", flush=True)
    expr, vel_dyn, vel_proxy, labels = compute_dynamical()
    print(f"velocity genes: {expr.shape[1]}  cells: {expr.shape[0]}", flush=True)

    dyn = evaluate(expr, vel_dyn, labels)
    print("dynamical:", dyn, flush=True)
    proxy = evaluate(expr, vel_proxy, labels)
    print("proxy:", proxy, flush=True)

    report = {
        "benchmark": "dynamical_velocity_validation",
        "dataset": "scVelo endocrinogenesis_day15 (endocrine branch)",
        "velocity_genes": int(expr.shape[1]),
        "dynamical": dyn,
        "proxy": proxy,
        "note": "Headline results reproduced on real scVelo dynamical-model velocity, alongside the proxy on the same cells and genes.",
    }
    out = Path("experiments/dynamical_velocity_validation.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print("wrote", out, flush=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    groups = ["dynamical", "proxy"]
    dp = [dyn["ep_directed_p"], proxy["ep_directed_p"]]
    sp = [dyn["ep_shuffled_p"], proxy["ep_shuffled_p"]]
    x = np.arange(2)
    ax1.bar(x - 0.2, dp, 0.4, label="directed", color="#0f7d99")
    ax1.bar(x + 0.2, sp, 0.4, label="shuffled control", color="#b0413e")
    ax1.axhline(0.05, color="black", ls="--", lw=1)
    ax1.set_xticks(x); ax1.set_xticklabels(groups)
    ax1.set_ylabel("entropy production p-value"); ax1.set_title("Irreversibility test"); ax1.legend(frameon=False)
    ro = [dyn["committor_ordering_rho"], proxy["committor_ordering_rho"]]
    ax2.bar(groups, ro, color="#0f7d99", width=0.5)
    for i, v in enumerate(ro):
        ax2.text(i, v + 0.01, f"{v:.2f}", ha="center")
    ax2.set_ylim(0, 1); ax2.set_ylabel("committor ordering (Spearman)"); ax2.set_title("Developmental ordering")
    fig.suptitle(f"Dynamical vs proxy velocity, pancreas branch ({dyn['n_cells']} cells)", y=1.02)
    fig.tight_layout()
    figpath = Path("figures/dynamical_velocity_validation.png")
    fig.savefig(figpath, dpi=150, bbox_inches="tight")
    print("wrote", figpath, flush=True)


if __name__ == "__main__":
    main()
