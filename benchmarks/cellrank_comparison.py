"""Head-to-head with CellRank on the pancreas endocrine branch.

A fair, apples-to-apples comparison on CellRank's own turf. Both methods compute
the same quantity: the probability of reaching the Beta terminal before falling
back to the Ductal progenitor. For CellRank that is a two-boundary absorption
probability (VelocityKernel + ConnectivityKernel + GPCCA, with Ductal and Beta
both absorbing); for this pipeline it is the transition-path-theory committor. We
score each by how well it recovers the known developmental stage, alongside PC1.

The point is not to beat CellRank at fate mapping (it is excellent at that); it is
to show the two agree on ordering while this pipeline additionally provides the
irreversibility test and kT barrier CellRank does not compute. Run at a few
thousand cells (the regime GPCCA is built for), not the earlier 150.

    python benchmarks/cellrank_comparison.py --cells 2000 --genes 1000 --seeds 3
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

# numpy 2.x compatibility shim for pygpcca (CellRank's GPCCA backend), which calls
# numpy.testing.assert_array_equal with the x=/y= keyword names numpy 2.x removed.
# Installed before cellrank is imported so the real CellRank code runs unchanged.
import numpy.testing as _npt

_orig_aae = _npt.assert_array_equal


def _aae_compat(*args, **kwargs):
    if "x" in kwargs or "y" in kwargs:
        a = kwargs.pop("x", None)
        b = kwargs.pop("y", None)
        return _orig_aae(a, b, **kwargs)
    return _orig_aae(*args, **kwargs)


_npt.assert_array_equal = _aae_compat

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
    scv.tl.velocity_graph(a, n_jobs=1)  # n_jobs=1 avoids a Windows multiprocessing hang
    return a


def cellrank_fate_to_beta(a):
    # Two-boundary absorption: make the progenitor (Ductal) and the terminal (Beta)
    # both absorbing and read the probability of reaching Beta first. That is exactly
    # the committor P(reach Beta before Ductal), so the two methods compute the same
    # quantity on the same graph. A single terminal is degenerate on a linear lineage
    # (CellRank warns "only 1 terminal state, all cells prob 1"), which is why the
    # earlier single-terminal setup gave a noisy, near-constant fate. solver="direct"
    # avoids the gmres non-convergence seen at larger n.
    import cellrank as cr
    import pandas as pd
    vk = cr.kernels.VelocityKernel(a).compute_transition_matrix()
    ck = cr.kernels.ConnectivityKernel(a).compute_transition_matrix()
    kernel = 0.8 * vk + 0.2 * ck
    g = cr.estimators.GPCCA(kernel)
    cl = a.obs["clusters"].astype(str)
    ts = pd.Series(index=a.obs_names, dtype="object")
    ts[(cl == "Ductal").to_numpy()] = "Ductal"
    ts[(cl == "Beta").to_numpy()] = "Beta"
    g.set_terminal_states(ts.astype("category"))
    g.compute_fate_probabilities(solver="direct")
    fate = g.fate_probabilities
    names = list(fate.names)
    beta = names.index("Beta")
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


def score_seed(cells, genes, seed):
    a = prepare(n_cells=cells, n_genes=genes, seed=seed)
    stage = np.array([BRANCH.index(l) for l in a.obs["clusters"].astype(str)], dtype=float)
    xc = dense(a.X) - dense(a.X).mean(axis=0)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    out = {"PC1": abs(float(spearmanr(xc @ vt[0], stage).correlation))}
    try:
        out["CellRank fate to Beta"] = abs(float(spearmanr(cellrank_fate_to_beta(a), stage).correlation))
    except Exception as exc:
        out["CellRank fate to Beta"] = None
        print("CellRank step failed:", type(exc).__name__, exc, flush=True)
    out["TW committor"] = abs(float(spearmanr(tw_committor(a), stage).correlation))
    return out


def _summary(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return {"mean": None, "n": 0}
    a = np.asarray(vals)
    return {"mean": round(float(a.mean()), 3), "std": round(float(a.std()), 3), "n": int(a.size)}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, default=150)
    ap.add_argument("--genes", type=int, default=500)
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()

    methods = ["PC1", "CellRank fate to Beta", "TW committor"]
    per = {m: [] for m in methods}
    for s in range(args.seeds):
        sc = score_seed(args.cells, args.genes, s)
        for m in methods:
            per[m].append(sc[m])
        print(f"  seed {s}: " + "  ".join(f"{m}={sc[m]}" for m in methods), flush=True)

    stats = {m: _summary(per[m]) for m in methods}
    report = {
        "benchmark": "cellrank_comparison_pancreas",
        "n_cells": args.cells,
        "seeds": args.seeds,
        "metric": "abs Spearman with developmental stage (mean over seeds)",
        "stats": stats,
        "note": "CellRank two-boundary absorption to Beta (Ductal and Beta both absorbing) and the committor compute the same quantity on the same graph, scored on recovering developmental stage (mean over seeds). At 2000 cells, the regime GPCCA is built for, CellRank wins: CellRank ~0.99 > committor ~0.95 > PC1 ~0.94. That is the expected result on a linear lineage with enough cells; the committor beats PC1 and tracks CellRank closely, and ordering is not this pipeline's claim. A single terminal is degenerate on a linear lineage (all cells prob 1), which is why an earlier 150-cell single-terminal setup gave a noisy, non-meaningful fate; this uses two boundaries and supersedes it. This pipeline additionally provides the entropy-production test and kT barrier CellRank does not compute. Running CellRank under numpy 2 needs the pygpcca shim and n_jobs=1 in this script.",
    }
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/cellrank_comparison.json").write_text(json.dumps(report, indent=2))
    for m in methods:
        print(f"  {m:26s} {stats[m]}")
    print("wrote experiments/cellrank_comparison.json")


if __name__ == "__main__":
    main()
