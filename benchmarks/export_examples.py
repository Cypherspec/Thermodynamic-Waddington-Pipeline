"""Write example per-cell exports (examples/exports/) for two real lineages.

    python benchmarks/export_examples.py

Each export is a cell table plus a summary (see thermodynamic_waddington.export),
ready to load into a tabular or map-based platform.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from thermodynamic_waddington.export import export_cells

OUT = Path("examples/exports")


def fly(cells=3000, seed=0):
    import anndata as ad
    import fly_testis as ft

    if not Path(ft.PATH).exists():
        print(f"skip fly: {ft.PATH} not found (run benchmarks/fly_testis_prepare.py)")
        return
    d = ft.load()
    idx = np.sort(np.random.default_rng(seed).choice(len(d["labels"]), size=min(cells, len(d["labels"])), replace=False))
    names = ad.read_h5ad(ft.PATH, backed="r").obs_names.to_numpy()[d["rows"][idx]]
    order = [k for k, _ in sorted(ft.STAGES.items(), key=lambda kv: kv[1])]
    res = export_cells(d["spliced"][idx], d["velocity"][idx], d["labels"][idx], ft.SOURCE, ft.TARGET, OUT,
                       prefix="fly_testis", coords=d["umap"][idx], cell_ids=[f"FCA59_{n}" for n in names],
                       stage_order=order, seed=seed,
                       dataset="Fly Cell Atlas testis germline (library FCA59), velocity from raw reads")
    s = res.summary
    print(f"fly: {s['n_cells']} cells, commits at {s['commitment']['commitment_stage']}, "
          f"barrier {s['commitment']['barrier_kt']} kT, noise-matched excess {s['irreversibility']['noise_matched_excess']}")


def pancreas(cells=1500, genes=100, seed=0):
    import anndata as ad
    from generalization import sparse_top_genes

    path = Path("data/real/endocrinogenesis_day15.h5ad")
    if not path.exists():
        print(f"skip pancreas: {path} not found")
        return
    order = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta"]
    a = ad.read_h5ad(path)
    lab = a.obs["clusters"].astype(str).to_numpy()
    keep = np.where(np.isin(lab, order))[0]
    keep = np.sort(np.random.default_rng(seed).choice(keep, size=min(cells, len(keep)), replace=False))
    a = a[keep].copy()
    top = sparse_top_genes(a.layers["spliced"], genes)
    sp = np.asarray(a.layers["spliced"][:, top].toarray(), dtype=float)
    un = np.asarray(a.layers["unspliced"][:, top].toarray(), dtype=float)
    res = export_cells(sp, un - sp, a.obs["clusters"].astype(str).to_numpy(), ["Ductal"], ["Beta"], OUT,
                       prefix="pancreas", coords=np.asarray(a.obsm["X_umap"]), cell_ids=list(a.obs_names),
                       stage_order=order, seed=seed, dataset="Pancreatic endocrinogenesis E15.5 (scVelo)")
    s = res.summary
    print(f"pancreas: {s['n_cells']} cells, commits at {s['commitment']['commitment_stage']}, "
          f"barrier {s['commitment']['barrier_kt']} kT, noise-matched excess {s['irreversibility']['noise_matched_excess']}")


if __name__ == "__main__":
    fly()
    pancreas()
