"""Does scVelo dynamical velocity carry more circulation than the crude proxy?

The cyclic-fraction finding (linear lineages near-reversible) used the
unspliced-minus-spliced proxy. This checks whether a proper scVelo dynamical
velocity fit changes the picture, on the same pancreas cells: if real velocity
carries circulation the proxy misses, the cyclic fraction should rise above the
reversible floor. Reported straight either way.

    python benchmarks/dynamical_cyclic.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.irreversibility import cyclic_irreversibility

BRANCH = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta"]
PATH = "data/real/endocrinogenesis_day15.h5ad"


def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=float)


def cyclic(expr, vel, dims=6):
    Xc = expr - expr.mean(0)
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    comp = vt[:dims].T
    pts = Xc @ comp
    g = build_knn(pts.tolist(), 20)
    real = cyclic_irreversibility(pts, vel @ comp, g).cyclic_fraction
    floor = cyclic_irreversibility(pts, (pts.mean(0) - pts), g).cyclic_fraction
    return round(real, 4), round(floor, 4)


def main(n_cells=800, n_genes=100, seed=0):
    import anndata as ad
    import scanpy as sc
    import scvelo as scv

    a = ad.read_h5ad(PATH)
    a = a[np.isin(a.obs["clusters"].astype(str).to_numpy(), BRANCH)].copy()
    idx = np.random.default_rng(seed).choice(a.n_obs, size=n_cells, replace=False)
    a = a[idx].copy()
    scv.pp.filter_genes(a, min_shared_counts=20)
    scv.pp.normalize_per_cell(a)
    sc.pp.log1p(a)
    top = np.argsort(dense(a.X).var(0))[::-1][:n_genes]
    a = a[:, top].copy()
    scv.pp.moments(a, n_pcs=30, n_neighbors=30)

    expr = dense(a.layers["Ms"])
    proxy = dense(a.layers["Mu"]) - dense(a.layers["Ms"])

    scv.tl.recover_dynamics(a, n_jobs=1)
    scv.tl.velocity(a, mode="dynamical")
    dyn = dense(a.layers["velocity"])

    # keep genes with a finite dynamical velocity in every cell
    keep = np.isfinite(dyn).all(axis=0)
    expr_k, proxy_k, dyn_k = expr[:, keep], proxy[:, keep], dyn[:, keep]

    proxy_real, proxy_floor = cyclic(expr_k, proxy_k)
    dyn_real, dyn_floor = cyclic(expr_k, dyn_k)
    report = {
        "benchmark": "dynamical_vs_proxy_cyclic_fraction",
        "dataset": "pancreas endocrine",
        "n_cells": n_cells, "n_genes_fit": int(keep.sum()),
        "proxy": {"cyclic_fraction": proxy_real, "floor": proxy_floor, "excess": round(proxy_real - proxy_floor, 4)},
        "dynamical": {"cyclic_fraction": dyn_real, "floor": dyn_floor, "excess": round(dyn_real - dyn_floor, 4)},
        "note": "Cyclic fraction of proxy vs scVelo dynamical velocity on the same pancreas cells. If dynamical velocity carried circulation the proxy misses, its excess over the reversible floor would be larger.",
    }
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/dynamical_cyclic.json").write_text(json.dumps(report, indent=2))
    print(f"proxy      cyclic={proxy_real}  floor={proxy_floor}  excess={proxy_real - proxy_floor:+.4f}")
    print(f"dynamical  cyclic={dyn_real}  floor={dyn_floor}  excess={dyn_real - dyn_floor:+.4f}  (genes fit: {int(keep.sum())})")
    print("wrote experiments/dynamical_cyclic.json")


if __name__ == "__main__":
    main()
