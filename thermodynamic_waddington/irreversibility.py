"""Calibrated irreversibility measure: the cyclic fraction of the velocity flow.

The Seifert estimator in `entropy_production.py` scores EP against the kNN density
proxy, not the stationary distribution of the estimated rates, so the net current
is not divergence-free and the statistic also counts the relaxation (gradient)
part of the flow. That part is large for any coherent field, including a
conservative, reversible one, so on a synthetic equilibrium field it
false-positives (see benchmarks/synthetic_ground_truth.py).

This measure is calibrated by construction via a discrete Hodge decomposition.
The flow along each graph edge, f_ij = 0.5 (v_i + v_j) . e_ij, is split into a
gradient (curl-free, reversible) part and a cyclic (divergence-free, irreversible)
part by projecting onto the graph Laplacian's range:

    f = B phi  +  f_cyclic,    phi = L^+ B^T f,   L = B^T B

The irreversibility is the cyclic energy fraction ||f_cyclic||^2 / ||f||^2 in
[0, 1]. A gradient field v = -grad U has f in the range of B, so the fraction is
~0 (up to a finite-sample discretization floor); a rotational field carries
circulation the gradient part cannot absorb, so the fraction is large. The
statistic is scale-invariant and needs no rate model or density proxy.

It is validated on ground truth by discrimination (equilibrium vs non-equilibrium
AUROC) and monotonic recovery of a known drive, not by a single-dataset p-value:
robustly testing broken detailed balance in finite, noisy, geometric data is an
open problem, so the equilibrium floor is reported as a calibration reference.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class IrreversibilityReport:
    cyclic_fraction: float
    bootstrap_ci: tuple[float, float] | None
    n_pairs: int
    n_cells: int

    def to_dict(self):
        return {
            "cyclic_fraction": self.cyclic_fraction,
            "bootstrap_ci": list(self.bootstrap_ci) if self.bootstrap_ci else None,
            "n_pairs": self.n_pairs,
            "n_cells": self.n_cells,
            "method": "hodge_cyclic_flow_fraction",
        }


def _pairs(graph):
    seen, out = set(), []
    for i, nbrs in enumerate(graph.neighbors):
        for j in nbrs:
            key = (i, j) if i < j else (j, i)
            if key not in seen:
                seen.add(key)
                out.append(key)
    return np.asarray(out, dtype=int) if out else np.zeros((0, 2), dtype=int)


def _edge_flow(pts, vels, pairs):
    i, j = pairs[:, 0], pairs[:, 1]
    e = pts[j] - pts[i]
    norm = np.linalg.norm(e, axis=1)
    norm[norm == 0] = 1.0
    e = e / norm[:, None]
    return 0.5 * (np.einsum("kd,kd->k", vels[i], e) + np.einsum("kd,kd->k", vels[j], e))


def _laplacian_pinv(pairs, n):
    L = np.zeros((n, n))
    i, j = pairs[:, 0], pairs[:, 1]
    np.add.at(L, (i, i), 1.0)
    np.add.at(L, (j, j), 1.0)
    np.add.at(L, (i, j), -1.0)
    np.add.at(L, (j, i), -1.0)
    return np.linalg.pinv(L)


def _cyclic_fraction(pairs, f, n, Lp):
    i, j = pairs[:, 0], pairs[:, 1]
    Btf = np.zeros(n)
    np.add.at(Btf, i, f)
    np.add.at(Btf, j, -f)
    phi = Lp @ Btf
    grad = phi[i] - phi[j]          # B phi
    cyclic = f - grad
    denom = float(np.sum(f * f)) + 1e-12
    return float(np.sum(cyclic * cyclic) / denom)


def cyclic_irreversibility(pts, vels, graph, bootstrap=0, seed=0):
    """Cyclic (irreversible) fraction of the velocity flow on the kNN graph.

    Returns a value in [0, 1]: ~0 for a conservative (reversible) field, larger
    for a field with real circulation. bootstrap>0 adds a resampled CI over edges.
    """
    pts = np.asarray(pts, dtype=float)
    vels = np.asarray(vels, dtype=float)
    n = len(pts)
    pairs = _pairs(graph)
    if len(pairs) == 0 or n < 3:
        return IrreversibilityReport(0.0, None, len(pairs), n)
    Lp = _laplacian_pinv(pairs, n)
    f = _edge_flow(pts, vels, pairs)
    obs = _cyclic_fraction(pairs, f, n, Lp)

    ci = None
    if bootstrap > 1:
        rng = np.random.default_rng(seed)
        m = len(pairs)
        vals = np.empty(bootstrap)
        for b in range(bootstrap):
            sel = rng.integers(0, m, size=m)
            vals[b] = _cyclic_fraction(pairs[sel], f[sel], n, Lp)
        ci = (float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975)))
    return IrreversibilityReport(obs, ci, len(pairs), n)
