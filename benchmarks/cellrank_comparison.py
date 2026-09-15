"""Head-to-head with CellRank on the pancreas endocrine branch.

A fair, apples-to-apples comparison on CellRank's own turf. Both methods estimate
a "probability of reaching the terminal fate": CellRank's fate probability to the
Beta terminal (VelocityKernel + GPCCA) and this pipeline's committor. We score
each by how well it recovers the known developmental stage, alongside PC1.

The point is not to beat CellRank at fate mapping (it is excellent at that); it is
to show the two agree on ordering while this pipeline additionally provides the
irreversibility test and kT barrier CellRank does not compute.

    python benchmarks/cellrank_comparison.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

BRANCH = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta"]
PATH = "data/real/endocrinogenesis_day15.h5ad"


def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=float)


def prepare(n_cells=600, n_genes=1000, seed=0):
    import anndata as ad
    import scanpy as sc
    import scvelo as scv
    a = ad.read_h5ad(PATH)
    a = a[np.isin(a.obs["clusters"].astype(str).to_numpy(), BRANCH)].copy()
    rng = np.random.default_rng(seed)
    if a.n_obs > n_cells:
        a = a[rng.choice(a.n_obs, size=n_cells, replace=False)].copy()
    scv.pp.filter_genes(a, min_shared_counts=20)
    scv.pp.normalize_per_cell(a)
    sc.pp.log1p(a)
    top = np.argsort(dense(a.X).var(axis=0))[::-1][:n_genes]
    a = a[:, top].copy()
    scv.pp.moments(a, n_pcs=30, n_neighbors=30)
    # unspliced-minus-spliced proxy on the smoothed moments, the same velocity the
    # pipeline uses, so CellRank and this pipeline compare on identical input.
    # (scVelo's stochastic/dynamical estimators are used elsewhere; the proxy keeps
    # this comparison fast and avoids a scVelo-numpy incompatibility.)
    a.layers["velocity"] = np.asarray(dense(a.layers["Mu"]) - dense(a.layers["Ms"]))
    scv.tl.velocity_graph(a)
    return a


def cellrank_fate_to_beta(a):
    import cellrank as cr
    vk = cr.kernels.VelocityKernel(a).compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(a).compute_transition_matrix()
    kernel = 0.8 * vk + 0.2 * ck
    g = cr.estimators.GPCCA(kernel)
    g.compute_schur(n_components=6)
    g.compute_macrostates(n_states=len(BRANCH), cluster_key="clusters")
    g.set_terminal_states(["Beta"])
    g.compute_fate_probabilities()
    fate = g.fate_probabilities
    names = list(fate.names)
    beta = names.index("Beta") if "Beta" in names else 0
    return np.asarray(fate.X[:, beta], dtype=float)


def tw_committor(a):
    from thermodynamic_waddington import FitConfig, developmental_coordinate, fit_landscape
    expr = dense(a.layers["Ms"]) if "Ms" in a.layers else dense(a.X)
    vel = np.nan_to_num(dense(a.layers["velocity"]))
    labels = a.obs["clusters"].astype(str).tolist()
    cfg = FitConfig(neighbors=20, dimensions=6, seed=0, enable_cycle_decomposition=False,
                    bootstrap_replicates=2, entropy_production_bootstrap_replicates=2,
                    entropy_production_permutation_replicates=2)
    fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=labels,
                        metadata={"velocity_status": "observed"})
    return np.asarray(developmental_coordinate(fit, ["Ductal"], ["Beta"]), dtype=float)


def main():
    a = prepare()
    stage = np.array([BRANCH.index(l) for l in a.obs["clusters"].astype(str)], dtype=float)
    xc = dense(a.X) - dense(a.X).mean(axis=0)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    pc1 = xc @ vt[0]

    scores = {"PC1": abs(float(spearmanr(pc1, stage).correlation))}
    cellrank_status = "ran"
    try:
        cr_fate = cellrank_fate_to_beta(a)
        scores["CellRank fate to Beta"] = abs(float(spearmanr(cr_fate, stage).correlation))
    except Exception as exc:
        scores["CellRank fate to Beta"] = None
        cellrank_status = f"could not run: {type(exc).__name__}: {exc}"
        print("CellRank step failed:", cellrank_status)
    scores["TW committor"] = abs(float(spearmanr(tw_committor(a), stage).correlation))

    report = {
        "benchmark": "cellrank_comparison_pancreas",
        "n_cells": int(a.n_obs),
        "metric": "abs Spearman with developmental stage",
        "scores": scores,
        "cellrank_status": cellrank_status,
        "note": "CellRank and the committor both estimate probability of reaching Beta; scored on ordering. This pipeline additionally provides the entropy-production test and kT barrier that CellRank does not. NOTE: CellRank 2.0 + numpy 2.x is currently incompatible (pygpcca calls a removed numpy signature); run this in a numpy<2 environment to get CellRank's number.",
    }
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/cellrank_comparison.json").write_text(json.dumps(report, indent=2))
    for k, v in scores.items():
        print(f"  {k:26s} {v}")
    print("wrote experiments/cellrank_comparison.json")


if __name__ == "__main__":
    main()
