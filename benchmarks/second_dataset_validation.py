"""Generalization to a second, independent dataset: gastrulation erythroid.

Every other result is pancreas. This runs the two headline results on the mouse
gastrulation erythroid lineage (Blood progenitors -> Erythroid): does the
committor still beat PC1 at ordering, and does entropy production still separate
directed from a shuffled control? Reported across seeds.

    python benchmarks/second_dataset_validation.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from thermodynamic_waddington import FitConfig, developmental_coordinate, fit_landscape
from thermodynamic_waddington.entropy_production import EntropyProductionConfig, estimate_entropy_production
from thermodynamic_waddington.graph import build_knn, estimate_local_diffusion, local_density
from thermodynamic_waddington.linalg import pca

PATH = "data/real/gastrulation_erythroid.h5ad"
ORDER = ["Blood progenitors 1", "Blood progenitors 2", "Erythroid1", "Erythroid2", "Erythroid3"]


def dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=float)


def load(n_cells, n_genes, seed):
    import anndata as ad
    a = ad.read_h5ad(PATH)
    lab = a.obs["celltype"].astype(str).to_numpy()
    keep = np.where(np.isin(lab, ORDER))[0]
    rng = np.random.default_rng(seed)
    if keep.size > n_cells:
        keep = rng.choice(keep, size=n_cells, replace=False)
    a = a[keep].copy()               # subsample cells before densifying
    lab = lab[keep]
    sp = dense(a.layers["spliced"])
    un = dense(a.layers["unspliced"])
    top = np.argsort(sp.var(axis=0))[::-1][:n_genes]
    return sp[:, top], (un - sp)[:, top], lab


def ep_p(expr, vel, seed, perms):
    pts, _, _ = pca(expr.tolist(), 6, seed)
    g = build_knn(pts, 20)
    dens = local_density(pts, g, 0.65)
    diff = estimate_local_diffusion(g, vel.tolist(), 1e-4)
    cfg = EntropyProductionConfig(temperature=1.0, velocity_scale=1.0, bootstrap_replicates=0,
                                  permutation_replicates=perms, seed=seed + 11)
    return estimate_entropy_production(pts, vel.tolist(), dens, diff, g, cfg).to_dict().get("permutation_p_value")


def run(seeds, cells, genes, perms):
    comm, pc, ep_dir, ep_shuf = [], [], [], []
    for s in range(seeds):
        expr, vel, lab = load(cells, genes, s)
        stage = np.array([ORDER.index(l) for l in lab], dtype=float)
        cfg = FitConfig(neighbors=20, dimensions=6, bootstrap_replicates=2, enable_cycle_decomposition=False,
                        entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2, seed=s)
        fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=lab.tolist(), metadata={"velocity_status": "derived_proxy"})
        q = np.array(developmental_coordinate(fit, [ORDER[0]], [ORDER[-1]]))
        xc = np.asarray(expr) - np.asarray(expr).mean(0)
        _, _, vt = np.linalg.svd(xc, full_matrices=False)
        comm.append(abs(float(spearmanr(q, stage).correlation)))
        pc.append(abs(float(spearmanr(xc @ vt[0], stage).correlation)))
        rng = np.random.default_rng(1000 + s)
        ep_dir.append(ep_p(expr, vel, s, perms))
        ep_shuf.append(ep_p(expr, vel[rng.permutation(len(lab))], s, perms))
        print(f"  seed {s}: committor={comm[-1]:.3f} PC1={pc[-1]:.3f} ep_dir={ep_dir[-1]:.3f} ep_shuf={ep_shuf[-1]:.3f}")
    return {
        "benchmark": "second_dataset_gastrulation_erythroid",
        "lineage": ORDER,
        "seeds": seeds,
        "committor_ordering_mean": round(float(np.mean(comm)), 3),
        "pc1_ordering_mean": round(float(np.mean(pc)), 3),
        "ep_directed_mean_p": round(float(np.mean(ep_dir)), 3),
        "ep_shuffled_mean_p": round(float(np.mean(ep_shuf)), 3),
        "committor_beats_pc1": int(np.sum(np.array(comm) > np.array(pc))),
        "note": "Independent dataset. Same pipeline, no tuning.",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--cells", type=int, default=500)
    ap.add_argument("--genes", type=int, default=50)
    ap.add_argument("--permutations", type=int, default=100)
    ap.add_argument("--out", default="experiments/second_dataset_validation.json")
    args = ap.parse_args()
    report = run(args.seeds, args.cells, args.genes, args.permutations)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "lineage"}, indent=2))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
