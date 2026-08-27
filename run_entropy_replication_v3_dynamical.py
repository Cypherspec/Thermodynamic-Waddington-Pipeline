"""Entropy-production replication v3: same design as v2 (90 cells/stage,
400 permutation replicates, Fisher's combined test across branches) but
with the velocity input swapped from the spliced-minus-unspliced proxy to
real scVelo dynamical-model velocity (EM-fit two-state transcriptional
kinetics; see fit_dynamical_velocity.py and
experiments/dynamical_velocity_fit_summary.json for how it was produced).

This is the item flagged as the "strongest remaining threat" to the v2
result in the research note: if the entropy-production signal only shows up
with the crude proxy and disappears (or changes sign/loses significance)
with model-based velocity, that would mean the proxy was doing the work,
not real dynamics. This script checks that directly instead of assuming
either outcome.

Gene set: of the 300 highly-variable genes used to fit the dynamical model,
255 converged to a usable kinetic fit. This script uses those 255 genes
directly (not re-selected by variance) since gene selection here has
already happened once, at the HVG-selection step, and re-selecting again by
variance on top of that would double-filter in a way that's harder to
reason about.
"""
from __future__ import annotations

import json
import math

import numpy as np
import anndata as ad

from thermodynamic_waddington.preprocessing import velocity_in_transformed_space
from thermodynamic_waddington.graph import build_knn, local_density, estimate_local_diffusion
from thermodynamic_waddington.entropy_production import (
    EntropyProductionConfig,
    estimate_entropy_production,
)
from run_entropy_replication_v2 import fishers_method

DATA_PATH = "data/real/endocrinogenesis_day15_dynamical.h5ad"
STAGE_ORDER = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine"]
N_PER_STAGE = 90
PERMUTATION_REPLICATES = 400
BOOTSTRAP_REPLICATES = 200
NEIGHBORS = 24


def load_branch_idx(adata: ad.AnnData, terminal_fate: str, seed: int):
    rng = np.random.default_rng(seed)
    stages = STAGE_ORDER + [terminal_fate]
    idx_by_stage = []
    for stage in stages:
        pool = np.where(adata.obs["clusters"].values == stage)[0]
        if len(pool) < N_PER_STAGE:
            raise ValueError(f"stage {stage} only has {len(pool)} cells, need {N_PER_STAGE}")
        chosen = rng.choice(pool, size=N_PER_STAGE, replace=False)
        idx_by_stage.append(chosen)
    idx = np.concatenate(idx_by_stage)
    return idx


def run_branch(adata: ad.AnnData, terminal_fate: str, seed: int, converged_genes) -> dict:
    idx = load_branch_idx(adata, terminal_fate, seed=seed)
    sub = adata[idx][:, converged_genes]

    expression = np.asarray(sub.layers["Ms"], dtype=float)   # smoothed spliced (scVelo's own moment estimate)
    velocity_raw = np.asarray(sub.layers["velocity"], dtype=float)

    # a handful of genes can carry NaN velocity for individual cells even
    # among "converged" genes (steady-state-only fits); drop any gene column
    # with NaNs in this subsample rather than silently imputing zeros, which
    # would fabricate a value the model didn't actually produce.
    good_cols = ~np.isnan(velocity_raw).any(axis=0)
    expression = expression[:, good_cols]
    velocity_raw = velocity_raw[:, good_cols]

    points, velocity, norm_report = velocity_in_transformed_space(
        expression.tolist(), velocity_raw.tolist()
    )

    graph = build_knn(points, NEIGHBORS)
    densities = local_density(points, graph, 0.65)
    diffusions = estimate_local_diffusion(graph, velocity, 1e-4)

    cfg = EntropyProductionConfig(
        bootstrap_replicates=BOOTSTRAP_REPLICATES,
        permutation_replicates=PERMUTATION_REPLICATES,
        seed=seed + 8181,
    )
    report = estimate_entropy_production(points, velocity, densities, diffusions, graph, cfg)

    return {
        "terminal_fate": terminal_fate,
        "n_cells": report.n_cells,
        "n_genes_used": int(good_cols.sum()),
        "entropy_production_rate": report.entropy_production_rate,
        "bootstrap_ci": list(report.bootstrap_ci) if report.bootstrap_ci else None,
        "permutation_p_value": report.permutation_p_value,
        "stationarity_residual": report.stationarity_residual,
    }


def main():
    adata = ad.read_h5ad(DATA_PATH)
    converged_genes = adata.var_names[~adata.var["fit_alpha"].isna()]
    print(f"{len(converged_genes)} genes with converged dynamical fits available", flush=True)

    results = {}
    for fate, seed in [("Beta", 7), ("Alpha", 1013)]:
        print(f"running branch -> {fate} (dynamical-model velocity)...", flush=True)
        results[fate] = run_branch(adata, fate, seed=seed, converged_genes=converged_genes)
        r = results[fate]
        print(f"  n_genes_used={r['n_genes_used']} "
              f"entropy_production={r['entropy_production_rate']:.3f} "
              f"CI={r['bootstrap_ci']} p={r['permutation_p_value']:.4f} "
              f"residual={r['stationarity_residual']:.4f}", flush=True)

    p_values = [results["Beta"]["permutation_p_value"], results["Alpha"]["permutation_p_value"]]
    statistic, p_combined = fishers_method(p_values)

    output = {
        "description": "Same design as v2 (90 cells/stage, 400 permutation replicates, "
                        "Fisher combined test) but with real scVelo dynamical-model "
                        "velocity instead of the spliced-minus-unspliced proxy.",
        "velocity_source": "scvelo dynamical model (EM-fit two-state kinetics)",
        "branches": results,
        "combined_test": {
            "method": "fisher",
            "input_p_values": p_values,
            "chi2_statistic": statistic,
            "degrees_of_freedom": 4,
            "combined_p_value": p_combined,
        },
    }
    with open("experiments/entropy_production_replication_v3_dynamical.json", "w") as f:
        json.dump(output, f, indent=2)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
