import sys, json, numpy as np
sys.path.insert(0, '.')
from run_independent_dataset_check import *

adata_fit = ad.read_h5ad("data/DentateGyrus/dentategyrus_dynamical.h5ad")
converged_genes = adata_fit.var_names[~adata_fit.var["fit_alpha"].isna()]

STAGES = ["Radial Glia-like", "Neuroblast", "Granule immature", "Granule mature"]

def run_single_age(adata, age, n_per_stage, seed):
    age_mask = adata.obs["age(days)"].astype(str) == str(age)
    rng = np.random.default_rng(seed)
    idx_by_stage = []
    for stage in STAGES:
        pool = np.where(age_mask.values & (adata.obs["clusters"].values == stage))[0]
        chosen = rng.choice(pool, size=n_per_stage, replace=False)
        idx_by_stage.append(chosen)
    idx = np.concatenate(idx_by_stage)
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
    cfg = EntropyProductionConfig(bootstrap_replicates=150, permutation_replicates=300, seed=seed+7777)
    report = estimate_entropy_production(points, velocity, densities, diffusions, graph, cfg)
    return {
        "age_days": age,
        "n_cells": report.n_cells,
        "n_genes_used": int(good_cols.sum()),
        "entropy_production_rate": report.entropy_production_rate,
        "bootstrap_ci": list(report.bootstrap_ci) if report.bootstrap_ci else None,
        "permutation_p_value": report.permutation_p_value,
        "stationarity_residual": report.stationarity_residual,
    }

r12 = run_single_age(adata_fit, 12, n_per_stage=30, seed=1200)
print(json.dumps(r12, indent=2), flush=True)
json.dump(r12, open("experiments/dentategyrus_single_timepoint_p12.json","w"), indent=2)
