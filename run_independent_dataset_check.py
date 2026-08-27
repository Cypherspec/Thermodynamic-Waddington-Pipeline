"""Independent-dataset check (item 1, the highest-priority item, from the
research note's next-steps list): does the entropy-production signal show
up in a completely different dataset -- different tissue, different lab,
different biology -- using the same method (real scVelo dynamical-model
velocity, same entropy-production estimator, same statistical design)?

Dataset: dentate gyrus neurogenesis (Hochgerner et al. 2018, via the scVelo
tutorial data mirror), NOT pancreatic endocrinogenesis. This is granule-cell
lineage development in the hippocampus: Radial Glia-like progenitors ->
nIPC -> Neuroblast -> Granule immature -> Granule mature. Nothing about this
biology, this lab, or this measurement is shared with the pancreas dataset
used in Addenda 1-3.

Honest limitation up front: the smallest stage (nIPC, n=19) forces
N_PER_STAGE=15, well below the 90/stage used for pancreas. This run will
necessarily be lower-powered. That's disclosed, not hidden -- a smaller,
honestly-reported effect here is more informative than pretending the
sample sizes match.
"""
from __future__ import annotations

import json
import time

import numpy as np
import scvelo as scv
import scanpy as sc
import anndata as ad

from thermodynamic_waddington.preprocessing import velocity_in_transformed_space
from thermodynamic_waddington.graph import build_knn, local_density, estimate_local_diffusion
from thermodynamic_waddington.entropy_production import (
    EntropyProductionConfig,
    estimate_entropy_production,
)

DATA_PATH = "data/DentateGyrus/10X43_1.h5ad"
STAGE_ORDER = ["Radial Glia-like", "nIPC", "Neuroblast", "Granule immature", "Granule mature"]
N_PER_STAGE = 15          # bottlenecked by nIPC (n=19 total)
N_TOP_GENES = 300
N_PCS = 30
N_NEIGHBORS_GRAPH = 30    # scanpy neighbor graph for moments, not the entropy-production kNN
PERMUTATION_REPLICATES = 400
BOOTSTRAP_REPLICATES = 200
EP_NEIGHBORS = 24          # matches the k that fit best in the Addendum 3 sensitivity sweep


def fit_dynamical(adata: ad.AnnData) -> ad.AnnData:
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)
    adata = adata[:, adata.var.highly_variable].copy()
    scv.pp.moments(adata, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS_GRAPH)
    scv.tl.recover_dynamics(adata, n_jobs=1)
    scv.tl.velocity(adata, mode="dynamical")
    return adata


def load_idx(adata: ad.AnnData, seed: int):
    rng = np.random.default_rng(seed)
    idx_by_stage = []
    for stage in STAGE_ORDER:
        pool = np.where(adata.obs["clusters"].values == stage)[0]
        if len(pool) < N_PER_STAGE:
            raise ValueError(f"stage {stage} only has {len(pool)} cells, need {N_PER_STAGE}")
        chosen = rng.choice(pool, size=N_PER_STAGE, replace=False)
        idx_by_stage.append(chosen)
    return np.concatenate(idx_by_stage)


def run(adata: ad.AnnData, seed: int) -> dict:
    converged_genes = adata.var_names[~adata.var["fit_alpha"].isna()]
    idx = load_idx(adata, seed=seed)
    sub = adata[idx][:, converged_genes]

    expression = np.asarray(sub.layers["Ms"], dtype=float)
    velocity_raw = np.asarray(sub.layers["velocity"], dtype=float)
    good_cols = ~np.isnan(velocity_raw).any(axis=0)
    expression = expression[:, good_cols]
    velocity_raw = velocity_raw[:, good_cols]

    points, velocity, _ = velocity_in_transformed_space(expression.tolist(), velocity_raw.tolist())

    graph = build_knn(points, EP_NEIGHBORS)
    densities = local_density(points, graph, 0.65)
    diffusions = estimate_local_diffusion(graph, velocity, 1e-4)

    cfg = EntropyProductionConfig(
        bootstrap_replicates=BOOTSTRAP_REPLICATES,
        permutation_replicates=PERMUTATION_REPLICATES,
        seed=seed + 9001,
    )
    report = estimate_entropy_production(points, velocity, densities, diffusions, graph, cfg)
    return {
        "dataset": "dentate_gyrus_hochgerner2018",
        "lineage": "granule cell (Radial Glia-like -> nIPC -> Neuroblast -> Granule immature -> Granule mature)",
        "n_cells": report.n_cells,
        "n_genes_used": int(good_cols.sum()),
        "n_genes_converged_total": int(len(converged_genes)),
        "entropy_production_rate": report.entropy_production_rate,
        "bootstrap_ci": list(report.bootstrap_ci) if report.bootstrap_ci else None,
        "permutation_p_value": report.permutation_p_value,
        "stationarity_residual": report.stationarity_residual,
    }


def main():
    adata = ad.read_h5ad(DATA_PATH)
    print(f"loaded {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    t0 = time.time()
    adata_fit = fit_dynamical(adata)
    print(f"dynamical fit done in {time.time()-t0:.1f}s, "
          f"{int((~adata_fit.var['fit_alpha'].isna()).sum())} genes converged", flush=True)
    adata_fit.write("data/DentateGyrus/dentategyrus_dynamical.h5ad")

    result = run(adata_fit, seed=2018)
    print(json.dumps(result, indent=2))
    with open("experiments/entropy_production_dentategyrus_independent.json", "w") as f:
        json.dump(result, f, indent=2)


if __name__ == "__main__":
    main()
