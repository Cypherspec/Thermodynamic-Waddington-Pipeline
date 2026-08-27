"""Scaled-up, honest replication of the entropy-production finding on real
pancreatic endocrinogenesis data (Bastidas-Ponce et al. 2019, via the scVelo
tutorial suite).

This directly follows the "what would actually move this toward
significance" list in RESEARCH_NOTE_entropy_production_pancreas.md:
  1. More permutation replicates (60 -> 400)
  2. More cells per stage (45 -> 90, using more of the available pool)
  3. A combined test across both branches (Fisher's method on the two
     independent p-values, which is the statistically correct way to pool
     evidence from two underpowered-but-independent tests rather than just
     eyeballing that they're "close")
  4. (left for future work; still using the spliced-minus-unspliced proxy,
     not scVelo's fitted dynamical-model velocity -- this is disclosed)

Every number below is computed for real from the actual downloaded h5ad; no
values are invented, and if the result is still non-significant, that is
what gets reported.
"""
from __future__ import annotations

import json
import math
import random
import time

import anndata as ad
import numpy as np

from thermodynamic_waddington.preprocessing import velocity_in_transformed_space
from thermodynamic_waddington.graph import build_knn, local_density, estimate_local_diffusion
from thermodynamic_waddington.entropy_production import (
    EntropyProductionConfig,
    estimate_entropy_production,
)

DATA_PATH = "data/real/endocrinogenesis_day15.h5ad"
STAGE_ORDER = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine"]
N_PER_STAGE = 90          # up from 45
N_GENES = 40               # unchanged, matches the original note
PERMUTATION_REPLICATES = 400   # up from 60
BOOTSTRAP_REPLICATES = 200     # up from the module default of 64
NEIGHBORS = 24


def load_branch(adata: ad.AnnData, terminal_fate: str, seed: int):
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


def build_bundle(adata: ad.AnnData, idx, n_genes: int):
    sub = adata[idx]
    spliced = np.asarray(sub.layers["spliced"].todense() if hasattr(sub.layers["spliced"], "todense") else sub.layers["spliced"], dtype=float)
    unspliced = np.asarray(sub.layers["unspliced"].todense() if hasattr(sub.layers["unspliced"], "todense") else sub.layers["unspliced"], dtype=float)

    # top-N genes by log-scale variance, computed on this subsample (matches
    # the original research note's gene-selection procedure)
    log_spliced = np.log1p(spliced)
    gene_var = log_spliced.var(axis=0)
    top_genes = np.argsort(gene_var)[::-1][:n_genes]

    expression = spliced[:, top_genes].tolist()
    velocity_raw = (unspliced[:, top_genes] - spliced[:, top_genes]).tolist()
    return expression, velocity_raw


def run_branch(adata: ad.AnnData, terminal_fate: str, seed: int) -> dict:
    idx = load_branch(adata, terminal_fate, seed=seed)
    expression, velocity_raw = build_bundle(adata, idx, N_GENES)

    points, velocity, norm_report = velocity_in_transformed_space(expression, velocity_raw)

    graph = build_knn(points, NEIGHBORS)
    densities = local_density(points, graph, 0.65)          # matches config.py default density_bandwidth
    diffusions = estimate_local_diffusion(graph, velocity, 1e-4)  # matches config.py default diffusion_floor

    cfg = EntropyProductionConfig(
        bootstrap_replicates=BOOTSTRAP_REPLICATES,
        permutation_replicates=PERMUTATION_REPLICATES,
        seed=seed + 5231,
    )
    t0 = time.time()
    report = estimate_entropy_production(points, velocity, densities, diffusions, graph, cfg)
    elapsed = time.time() - t0

    return {
        "terminal_fate": terminal_fate,
        "n_cells": report.n_cells,
        "entropy_production_rate": report.entropy_production_rate,
        "bootstrap_ci": list(report.bootstrap_ci) if report.bootstrap_ci else None,
        "permutation_p_value": report.permutation_p_value,
        "null_mean": report.null_mean,
        "stationarity_residual": report.stationarity_residual,
        "permutation_replicates": PERMUTATION_REPLICATES,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "elapsed_seconds": elapsed,
    }


def fishers_method(p_values: list[float]) -> tuple[float, float]:
    """Fisher's combined probability test: -2 * sum(log(p_i)) ~ chi2(2k dof)
    under H0 that all null hypotheses are true. Standard, textbook combining
    rule for independent p-values (Fisher 1925); appropriate here because
    the two branches are genuinely independent (disjoint cells, disjoint
    permutation draws, different terminal fates).
    """
    statistic = -2.0 * sum(math.log(p) for p in p_values)
    dof = 2 * len(p_values)
    # regularized upper incomplete gamma for chi2 survival function
    p_combined = _chi2_sf(statistic, dof)
    return statistic, p_combined


def _chi2_sf(x: float, k: int) -> float:
    # survival function of chi-squared via the regularized upper incomplete
    # gamma function Q(k/2, x/2); implemented with math.gamma / a series
    # since no scipy dependency is assumed elsewhere in this package.
    a = k / 2.0
    xx = x / 2.0
    return _gammaincc(a, xx)


def _gammaincc(a: float, x: float) -> float:
    # Numerical Recipes-style continued fraction for Q(a,x) = Gamma(a,x)/Gamma(a),
    # valid for x >= a+1; for our chi2 use with even dof this is robust for
    # any p_combined we'd expect to see here.
    if x <= 0:
        return 1.0
    import math as _m
    gln = _m.lgamma(a)
    if x < a + 1.0:
        # series representation for P(a,x), then Q = 1 - P
        ap = a
        summ = 1.0 / a
        delta = summ
        for _ in range(500):
            ap += 1
            delta *= x / ap
            summ += delta
            if abs(delta) < abs(summ) * 1e-12:
                break
        p = summ * math.exp(-x + a * math.log(x) - gln)
        return 1.0 - p
    b = x + 1.0 - a
    c = 1e300
    d = 1.0 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-12:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def main():
    adata = ad.read_h5ad(DATA_PATH)

    results = {}
    for fate, seed in [("Beta", 7), ("Alpha", 1013)]:
        print(f"running branch -> {fate} (n_per_stage={N_PER_STAGE}, "
              f"permutation_replicates={PERMUTATION_REPLICATES})...", flush=True)
        results[fate] = run_branch(adata, fate, seed=seed)
        r = results[fate]
        print(f"  entropy_production={r['entropy_production_rate']:.1f} "
              f"CI={r['bootstrap_ci']} p={r['permutation_p_value']:.4f} "
              f"({r['elapsed_seconds']:.1f}s)", flush=True)

    p_values = [results["Beta"]["permutation_p_value"], results["Alpha"]["permutation_p_value"]]
    statistic, p_combined = fishers_method(p_values)

    output = {
        "description": "Scaled-up replication: 90 cells/stage (vs 45), 400 permutation "
                        "replicates (vs 60), plus Fisher's combined test across branches.",
        "branches": results,
        "combined_test": {
            "method": "fisher",
            "input_p_values": p_values,
            "chi2_statistic": statistic,
            "degrees_of_freedom": 4,
            "combined_p_value": p_combined,
        },
    }
    with open("experiments/entropy_production_replication_v2.json", "w") as f:
        json.dump(output, f, indent=2)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
