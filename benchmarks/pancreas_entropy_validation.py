"""Real-data entropy-production validation on pancreatic endocrinogenesis.

Loads the public scVelo endocrinogenesis dataset, restricts to the endocrine
differentiation branch (Ductal -> Ngn3 -> Pre-endocrine -> Beta), builds a
transparent unspliced-minus-spliced velocity proxy on the most variable genes,
and runs the landscape fit. The pipeline estimates entropy production with a
label-permutation null; a system with directed differentiation flow should show
entropy production that the shuffled-velocity null cannot reproduce.

This uses the proxy velocity, not a scVelo dynamical fit, so it validates the
estimator and its significance test on real expression dynamics. The dynamical
-velocity result on the Beta and Alpha branches is in the research note.

    python benchmarks/pancreas_entropy_validation.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.model import fit_landscape

BRANCH = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta"]


def _dense(x):
    return x.toarray() if hasattr(x, "toarray") else np.asarray(x)


def load_branch(path, branch, n_cells, n_genes, seed):
    import anndata as ad
    a = ad.read_h5ad(path)
    labels = a.obs["clusters"].astype(str).to_numpy()
    keep = np.isin(labels, branch)
    a = a[keep].copy()
    labels = labels[keep]
    spliced = _dense(a.layers["spliced"]).astype(float)
    unspliced = _dense(a.layers["unspliced"]).astype(float)
    var = spliced.var(axis=0)
    top = np.argsort(var)[::-1][:n_genes]
    spliced = spliced[:, top]
    unspliced = unspliced[:, top]
    rng = np.random.default_rng(seed)
    if spliced.shape[0] > n_cells:
        idx = rng.choice(spliced.shape[0], size=n_cells, replace=False)
        spliced, unspliced, labels = spliced[idx], unspliced[idx], labels[idx]
    velocity = unspliced - spliced
    return spliced, velocity, labels.tolist()


def run(path, n_cells, n_genes, seed, permutations):
    expr, vel, labels = load_branch(path, BRANCH, n_cells, n_genes, seed)
    cfg = FitConfig(
        neighbors=20,
        dimensions=6,
        bootstrap_replicates=8,
        entropy_production_bootstrap_replicates=24,
        entropy_production_permutation_replicates=permutations,
        seed=seed,
    )
    meta = {"velocity_status": "derived_proxy", "gene_names": [f"g{i}" for i in range(n_genes)]}
    fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=labels, metadata=meta)
    ep = fit.diagnostics.get("entropy_production", {})

    # explicit negative control: destroy the velocity-displacement coupling by
    # permuting which cell each velocity vector belongs to, then re-estimate
    rng = np.random.default_rng(seed + 1)
    perm = rng.permutation(len(labels))
    shuffled = fit_landscape(expr.tolist(), vel[perm].tolist(), config=cfg, labels=labels, metadata=meta)
    ep_shuf = shuffled.diagnostics.get("entropy_production", {})

    report = {
        "benchmark": "pancreas_entropy_production_proxy",
        "dataset": "scVelo endocrinogenesis_day15",
        "branch": BRANCH,
        "n_cells": len(labels),
        "n_genes": n_genes,
        "observed": {
            "entropy_production_rate": ep.get("entropy_production_rate"),
            "permutation_p_value": ep.get("permutation_p_value"),
            "null_mean": ep.get("null_mean"),
            "bootstrap_ci": ep.get("bootstrap_ci"),
            "n_pairs": ep.get("n_pairs"),
        },
        "shuffled_velocity_control": {
            "entropy_production_rate": ep_shuf.get("entropy_production_rate"),
            "permutation_p_value": ep_shuf.get("permutation_p_value"),
        },
        "velocity": "unspliced minus spliced proxy (not scVelo dynamical)",
        "note": "Proxy-velocity execution and significance check. Dynamical-velocity Beta/Alpha result is in RESEARCH_NOTE_entropy_production_pancreas.md.",
    }
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--path", default="data/real/endocrinogenesis_day15.h5ad")
    ap.add_argument("--cells", type=int, default=500)
    ap.add_argument("--genes", type=int, default=50)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--permutations", type=int, default=200)
    ap.add_argument("--out", default="experiments/pancreas_entropy_validation.json")
    args = ap.parse_args()
    report = run(args.path, args.cells, args.genes, args.seed, args.permutations)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    obs = report["observed"]
    ctl = report["shuffled_velocity_control"]
    print(f"cells={report['n_cells']} genes={report['n_genes']}")
    print(f"observed  EP={obs['entropy_production_rate']}  p={obs['permutation_p_value']}")
    print(f"shuffled  EP={ctl['entropy_production_rate']}  p={ctl['permutation_p_value']}")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
