# PRE-SPECIFIED CONFIRMATORY TEST
#
# This follows up the exploratory Radial Glia-like -> Granule mature result
# from the two-fate comparison (p=0.047, uncorrected, one of six tests run
# on that dataset -- see RESEARCH_NOTE Addendum 6). The design below was
# fixed BEFORE this script was run, specifically to avoid picking parameters
# after looking at results, which is the same mistake the multiple-comparisons
# caveat in Addendum 6 was written to guard against.
#
# Fixed design:
#   lineage:                 Radial Glia-like -> Granule mature
#   n_per_stage:              45 (ceiling both stages support)
#   k (kNN neighbors):        24 (the value that fit best in the sensitivity sweep)
#   permutation_replicates:   1000 (up from 300, for a non-floor-clipped p-value;
#                              floor at 1000 replicates is 1/1001 = 0.000999)
#   seed:                     424242 (NEW -- draws a different 90-cell subsample
#                              than the exploratory run, so this is a genuine
#                              held-out replication, not the same cells re-run
#                              with more permutations)
#   significance threshold:   0.01, not 0.05 -- stricter than usual because this
#                              is a confirmatory follow-up to an exploratory hit,
#                              which is exactly the situation where using the same
#                              threshold that produced the hit risks confirming
#                              noise (regression to the mean / winner's curse)
from __future__ import annotations

import json

import numpy as np
import anndata as ad

from thermodynamic_waddington.preprocessing import velocity_in_transformed_space
from thermodynamic_waddington.graph import build_knn, local_density, estimate_local_diffusion
from thermodynamic_waddington.entropy_production import (
    EntropyProductionConfig,
    estimate_entropy_production,
)

DATA_PATH = "data/DentateGyrus/dentategyrus_dynamical.h5ad"
N_PER_STAGE = 45
NEIGHBORS = 24
PERMUTATION_REPLICATES = 1000
BOOTSTRAP_REPLICATES = 300
SEED = 424242
SIGNIFICANCE_THRESHOLD = 0.01


def main():
    adata = ad.read_h5ad(DATA_PATH)
    converged_genes = adata.var_names[~adata.var["fit_alpha"].isna()]

    rng = np.random.default_rng(SEED)
    idx_by_stage = []
    for stage in ["Radial Glia-like", "Granule mature"]:
        pool = np.where(adata.obs["clusters"].values == stage)[0]
        idx_by_stage.append(rng.choice(pool, size=N_PER_STAGE, replace=False))
    idx = np.concatenate(idx_by_stage)

    sub = adata[idx][:, converged_genes]
    expression = np.asarray(sub.layers["Ms"], dtype=float)
    velocity_raw = np.asarray(sub.layers["velocity"], dtype=float)
    good_cols = ~np.isnan(velocity_raw).any(axis=0)
    expression = expression[:, good_cols]
    velocity_raw = velocity_raw[:, good_cols]

    points, velocity, _ = velocity_in_transformed_space(expression.tolist(), velocity_raw.tolist())
    graph = build_knn(points, NEIGHBORS)
    densities = local_density(points, graph, 0.65)
    diffusions = estimate_local_diffusion(graph, velocity, 1e-4)

    cfg = EntropyProductionConfig(
        bootstrap_replicates=BOOTSTRAP_REPLICATES,
        permutation_replicates=PERMUTATION_REPLICATES,
        seed=SEED + 1,
    )
    report = estimate_entropy_production(points, velocity, densities, diffusions, graph, cfg)

    result = {
        "design": "pre-specified confirmatory test, fixed before execution",
        "lineage": "Radial Glia-like -> Granule mature",
        "n_per_stage": N_PER_STAGE,
        "k": NEIGHBORS,
        "permutation_replicates": PERMUTATION_REPLICATES,
        "seed": SEED,
        "significance_threshold": SIGNIFICANCE_THRESHOLD,
        "n_cells": report.n_cells,
        "n_genes_used": int(good_cols.sum()),
        "entropy_production_rate": report.entropy_production_rate,
        "bootstrap_ci": list(report.bootstrap_ci) if report.bootstrap_ci else None,
        "permutation_p_value": report.permutation_p_value,
        "stationarity_residual": report.stationarity_residual,
        "passes_confirmatory_threshold": report.permutation_p_value < SIGNIFICANCE_THRESHOLD,
    }
    print(json.dumps(result, indent=2))
    with open("experiments/confirmatory_granule_test.json", "w") as f:
        json.dump(result, f, indent=2)


if __name__ == "__main__":
    main()
