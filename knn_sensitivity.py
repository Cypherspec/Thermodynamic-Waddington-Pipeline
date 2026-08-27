"""Sensitivity check (item 3 from the research note's next-steps list):
does the entropy-production signal (real scVelo dynamical-model velocity)
survive changing the kNN graph's neighbor count, or is it sitting on a
knife-edge of k=24 specifically?

Design: same 450-cell Beta and Alpha branches used in v3 (same seeds, same
cell indices -- this is NOT a fresh resample, deliberately, so any change in
result is attributable to k alone, not also to a different subsample).
Sweep k in {12, 18, 24, 32, 40}. Permutation replicates reduced to 150 for
this sweep (10 total branch x k combinations at 400 replicates each would
be prohibitively slow single-core); this trades permutation resolution for
coverage across k -- appropriate for a sensitivity check where the question
is "does the qualitative conclusion hold," not "what is the precise p-value"
(that's what the dedicated v3 run at k=24, 400 replicates, is for).
"""
from __future__ import annotations

import json
import time

import numpy as np
import anndata as ad

from thermodynamic_waddington.preprocessing import velocity_in_transformed_space
from thermodynamic_waddington.graph import build_knn, local_density, estimate_local_diffusion
from thermodynamic_waddington.entropy_production import (
    EntropyProductionConfig,
    estimate_entropy_production,
)
from run_entropy_replication_v3_dynamical import load_branch_idx, DATA_PATH

K_VALUES = [12, 18, 24, 32, 40]
PERMUTATION_REPLICATES = 150
BOOTSTRAP_REPLICATES = 100


def run_one(adata, terminal_fate: str, seed: int, converged_genes, k: int) -> dict:
    idx = load_branch_idx(adata, terminal_fate, seed=seed)
    sub = adata[idx][:, converged_genes]

    expression = np.asarray(sub.layers["Ms"], dtype=float)
    velocity_raw = np.asarray(sub.layers["velocity"], dtype=float)
    good_cols = ~np.isnan(velocity_raw).any(axis=0)
    expression = expression[:, good_cols]
    velocity_raw = velocity_raw[:, good_cols]

    points, velocity, _ = velocity_in_transformed_space(expression.tolist(), velocity_raw.tolist())

    graph = build_knn(points, k)
    densities = local_density(points, graph, 0.65)
    diffusions = estimate_local_diffusion(graph, velocity, 1e-4)

    cfg = EntropyProductionConfig(
        bootstrap_replicates=BOOTSTRAP_REPLICATES,
        permutation_replicates=PERMUTATION_REPLICATES,
        seed=seed + 4242 + k,
    )
    report = estimate_entropy_production(points, velocity, densities, diffusions, graph, cfg)
    return {
        "k": k,
        "entropy_production_rate": report.entropy_production_rate,
        "permutation_p_value": report.permutation_p_value,
        "stationarity_residual": report.stationarity_residual,
    }


def main():
    adata = ad.read_h5ad(DATA_PATH)
    converged_genes = adata.var_names[~adata.var["fit_alpha"].isna()]

    results = {"Beta": [], "Alpha": []}
    for fate, seed in [("Beta", 7), ("Alpha", 1013)]:
        for k in K_VALUES:
            t0 = time.time()
            r = run_one(adata, fate, seed, converged_genes, k)
            r["elapsed_seconds"] = round(time.time() - t0, 1)
            results[fate].append(r)
            print(f"{fate} k={k}: EP={r['entropy_production_rate']:.3f} "
                  f"p={r['permutation_p_value']:.4f} "
                  f"residual={r['stationarity_residual']:.4f} "
                  f"({r['elapsed_seconds']}s)", flush=True)

    with open("experiments/knn_sensitivity_dynamical.json", "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
