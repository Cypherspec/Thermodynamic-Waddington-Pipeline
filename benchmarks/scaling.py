"""Scaling benchmark for the graph core and the fit.

Two things are measured:

1. build_knn: graph construction. This used to be an O(n^2) pure-Python
   pairwise-distance pass, the hard wall on dataset size. It is now a KD-tree
   (n log n) and scales to tens of thousands of cells in seconds, with output
   bit-identical to the old brute force.

2. Core fit (cycle decomposition off): the full landscape plus entropy
   production. It is faster and no longer crashes on large graphs (the cycle
   SCC pass is now iterative), but it is still super-linear because the
   free-energy propagation and bootstrap are. Practical to a few thousand
   cells; the propagation is the next scaling target. This is reported honestly
   rather than claimed as atlas-scale.

    python benchmarks/scaling.py
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.model import fit_landscape
from thermodynamic_waddington.synthetic import make_synthetic_dataset


def time_knn(sizes, k=20, dims=6, seed=0):
    rows = []
    rng = np.random.default_rng(seed)
    build_knn(rng.normal(size=(64, dims)).tolist(), k)  # warm up the KD-tree backend
    for n in sizes:
        x = rng.normal(size=(n, dims)).tolist()
        t = time.perf_counter()
        g = build_knn(x, k)
        rows.append({"n": n, "seconds": round(time.perf_counter() - t, 3), "edges": len(g.edges)})
        print(f"  knn  n={n:6d}  {rows[-1]['seconds']:7.2f}s")
    return rows


def time_fit(sizes, seed=7):
    rows = []
    for n in sizes:
        ds = make_synthetic_dataset(cells=n, genes=16, seed=seed)
        cfg = FitConfig(
            neighbors=20, dimensions=6, bootstrap_replicates=2,
            enable_cycle_decomposition=False,
            entropy_production_bootstrap_replicates=2,
            entropy_production_permutation_replicates=2,
        )
        t = time.perf_counter()
        fit_landscape(ds.expression, ds.velocity, config=cfg, labels=ds.labels)
        rows.append({"n": n, "seconds": round(time.perf_counter() - t, 2)})
        print(f"  fit  n={n:6d}  {rows[-1]['seconds']:7.2f}s")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--knn-sizes", type=int, nargs="+", default=[1000, 4000, 16000, 50000])
    ap.add_argument("--fit-sizes", type=int, nargs="+", default=[500, 1000, 2000, 4000])
    ap.add_argument("--out", default="experiments/scaling.json")
    ap.add_argument("--figure", default="figures/scaling.png")
    args = ap.parse_args()

    print("=== build_knn (KD-tree) ===")
    knn = time_knn(args.knn_sizes)
    print("=== core fit (cycle decomposition off) ===")
    fit = time_fit(args.fit_sizes)
    report = {
        "benchmark": "scaling",
        "knn": knn,
        "core_fit": fit,
        "notes": "build_knn is n log n and bit-identical to the old brute force; core fit is still super-linear (propagation + bootstrap).",
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    ax1.plot([r["n"] for r in knn], [r["seconds"] for r in knn], "o-", color="#0f7d99", lw=2, ms=7)
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("cells"); ax1.set_ylabel("seconds")
    ax1.set_title("Graph construction (KD-tree), to 50k cells")
    ax1.grid(True, which="both", alpha=0.3)
    ax2.plot([r["n"] for r in fit], [r["seconds"] for r in fit], "o-", color="#6a4c93", lw=2, ms=7)
    ax2.set_xlabel("cells"); ax2.set_ylabel("seconds")
    ax2.set_title("Core fit (landscape + EP)")
    ax2.grid(True, alpha=0.3)
    fig.suptitle("Scaling: graph core is n log n; full fit practical to a few thousand cells", y=1.02)
    fig.tight_layout()
    Path(args.figure).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure, dpi=150, bbox_inches="tight")
    print("wrote", args.out, "and", args.figure)


if __name__ == "__main__":
    main()
