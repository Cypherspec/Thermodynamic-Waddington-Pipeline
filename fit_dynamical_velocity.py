"""Fit real scVelo dynamical-model velocity on the pancreatic endocrinogenesis
data, replacing the spliced-minus-unspliced proxy used in the original
research note and its v2 replication.

What "dynamical model" actually means here (this is the real math, not a
marketing term): for each gene, scVelo fits a two-state transcriptional
kinetics model

    du/dt = alpha(t) - beta * u
    ds/dt = beta * u - gamma * s

(u = unspliced counts, s = spliced counts, alpha = transcription rate,
beta = splicing rate, gamma = degradation rate) under a piecewise-constant
"on/off" transcriptional-state assumption, and estimates each cell's latent
time t and each gene's (alpha, beta, gamma, switching time) by maximum
likelihood via an EM algorithm (Bergen et al. 2020, Nature Biotechnology).
Velocity for a cell/gene is then the fitted ds/dt at that cell's inferred
latent time -- a model-based estimate of transcriptional dynamics, not a
finite-difference proxy.

This is fit once on a broader cell pool (for numerical stability of the EM
fit) and then the same branch-specific cells used in the v2 replication are
pulled back out, so the entropy-production comparison stays apples-to-apples
with the earlier note.
"""
from __future__ import annotations

import json
import time

import numpy as np
import scvelo as scv
import scanpy as sc
import anndata as ad

DATA_PATH = "data/real/endocrinogenesis_day15.h5ad"
N_TOP_GENES = 300   # kept modest: single-core EM fit, one gene at a time
N_PCS = 30
N_NEIGHBORS = 30


def main():
    scv.settings.verbosity = 2
    adata = ad.read_h5ad(DATA_PATH)
    print(f"loaded {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

    t0 = time.time()
    # scv.pp.filter_and_normalize's convenience n_top_genes kwarg is broken in
    # scvelo 0.3.4 against modern scanpy (it forwards n_top_genes into
    # normalize_per_cell, which doesn't accept it) -- so the steps are done
    # explicitly here instead: filter by shared counts, normalize, log-transform,
    # then select highly-variable genes via scanpy directly.
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)
    adata = adata[:, adata.var.highly_variable].copy()
    scv.pp.moments(adata, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS)
    print(f"preprocessing done in {time.time()-t0:.1f}s, "
          f"{adata.n_vars} genes retained after filtering", flush=True)

    t0 = time.time()
    scv.tl.recover_dynamics(adata, n_jobs=1)
    print(f"recover_dynamics (EM fit per gene) done in {time.time()-t0:.1f}s", flush=True)

    t0 = time.time()
    scv.tl.velocity(adata, mode="dynamical")
    print(f"dynamical-model velocity computed in {time.time()-t0:.1f}s", flush=True)

    n_genes_fit = int((~adata.var["fit_alpha"].isna()).sum()) if "fit_alpha" in adata.var else None
    n_velocity_genes = int(np.sum(~np.isnan(adata.layers["velocity"]).all(axis=0)))

    print(f"genes with a converged kinetic fit: {n_genes_fit}", flush=True)
    print(f"genes with usable velocity values: {n_velocity_genes}", flush=True)

    adata.write("data/real/endocrinogenesis_day15_dynamical.h5ad")

    summary = {
        "n_cells": int(adata.n_obs),
        "n_genes_after_filter": int(adata.n_vars),
        "n_genes_with_converged_fit": n_genes_fit,
        "n_genes_with_usable_velocity": n_velocity_genes,
        "mode": "dynamical",
        "model": "two-state transcriptional kinetics (Bergen et al. 2020), "
                 "EM-fit alpha/beta/gamma/switching-time per gene, "
                 "velocity = fitted ds/dt at each cell's inferred latent time",
    }
    with open("experiments/dynamical_velocity_fit_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
