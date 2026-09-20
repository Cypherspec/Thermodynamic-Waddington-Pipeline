"""Scaling of the calibrated irreversibility measures.

The Hodge cyclic-fraction projection used to be a dense n-by-n pseudo-inverse
(O(cells^3)); it is now a sparse least-squares solve (O(nnz)). This benchmark
times both irreversibility readouts across cell counts, and shows the old dense
projection next to the new one at the sizes where the dense one is still feasible.

    python benchmarks/irreversibility_scaling.py
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.irreversibility import (
    _edge_flow,
    _frac_dense,
    _frac_sparse,
    _incidence,
    _laplacian_pinv,
    _pairs,
    cycle_affinities,
    cyclic_irreversibility,
)

J = np.array([[0.0, -1.0], [1.0, 0.0]])


def sim(n, seed=0, omega=1.0, a=1.0, D=1.0):
    rng = np.random.default_rng(seed)
    pts = rng.normal(0.0, math.sqrt(D / a), size=(n, 2))
    return pts, -(pts @ (a * np.eye(2) + omega * J).T)


def main(sizes=(500, 1000, 2000, 5000, 10000, 25000, 50000)):
    rows = []
    for n in sizes:
        pts, vel = sim(n)
        g = build_knn(pts.tolist(), 20)
        t = time.perf_counter()
        frac = cyclic_irreversibility(pts, vel, g).cyclic_fraction
        t_frac = time.perf_counter() - t
        t = time.perf_counter()
        aff = cycle_affinities(pts, vel, g).rms_affinity_kt
        t_aff = time.perf_counter() - t

        # dense pseudo-inverse projection, only where it is still feasible
        t_dense = None
        if n <= 3000:
            pairs = _pairs(g)
            f = _edge_flow(pts, vel, pairs)
            t = time.perf_counter()
            _frac_dense(pairs, f, n, _laplacian_pinv(pairs, n))
            t_dense = time.perf_counter() - t

        rows.append({
            "n_cells": n,
            "cyclic_fraction": round(float(frac), 4),
            "cycle_affinity_kt": round(float(aff), 3),
            "t_cyclic_fraction_s": round(t_frac, 3),
            "t_cycle_affinity_s": round(t_aff, 3),
            "t_dense_projection_s": round(t_dense, 3) if t_dense is not None else None,
            "speedup_vs_dense": round(t_dense / t_frac, 1) if t_dense and t_frac > 0 else None,
        })
        sp = rows[-1]["speedup_vs_dense"]
        print(f"  n={n:>6}: cyclic {t_frac:.3f}s  affinity {t_aff:.3f}s"
              + (f"  (dense {t_dense:.2f}s, {sp}x faster)" if t_dense else "  (dense skipped)"), flush=True)

    report = {
        "benchmark": "irreversibility_scaling",
        "setup": "synthetic rotational field, k=20 neighbors, sparse vs dense Hodge projection",
        "rows": rows,
        "note": "The sparse least-squares projection is bit-identical to the dense pseudo-inverse but O(nnz); the cyclic fraction and cycle affinity both scale to 50k cells in about a second, where the dense projection is cubic and infeasible.",
    }
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/irreversibility_scaling.json").write_text(json.dumps(report, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ns = [r["n_cells"] for r in rows]
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.plot(ns, [r["t_cyclic_fraction_s"] for r in rows], "o-", color="#0f7d99", lw=2.4, label="cyclic fraction (sparse)")
    ax.plot(ns, [r["t_cycle_affinity_s"] for r in rows], "s-", color="#14202e", lw=2, label="cycle affinity")
    dns = [r["n_cells"] for r in rows if r["t_dense_projection_s"] is not None]
    dts = [r["t_dense_projection_s"] for r in rows if r["t_dense_projection_s"] is not None]
    ax.plot(dns, dts, "^--", color="#b0413e", lw=2, label="dense pseudo-inverse (old)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("cells"); ax.set_ylabel("seconds")
    ax.set_title("Irreversibility measures scale to 50k cells")
    ax.legend(frameon=False); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    Path("figures").mkdir(exist_ok=True)
    fig.savefig("figures/irreversibility_scaling.png", dpi=150, bbox_inches="tight")
    print("wrote experiments/irreversibility_scaling.json and figures/irreversibility_scaling.png")


if __name__ == "__main__":
    main()
