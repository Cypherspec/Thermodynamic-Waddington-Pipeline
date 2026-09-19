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

from collections import deque
from dataclasses import dataclass

import numpy as np


@dataclass
class CycleAffinityReport:
    """Schnakenberg cycle affinities of the velocity flow, in kT units.

    Each fundamental cycle of the kNN graph carries a thermodynamic affinity
    A = sum around the cycle of log(k_ij / k_ji), the entropy produced per turn of
    that cycle. With rates from the edge-projected velocity this equals twice the
    circulation of the flow, so it is exactly zero for a gradient (reversible)
    field (Kolmogorov's criterion) and non-zero when the flow has curl.
    """

    n_cycles: int
    rms_affinity_kt: float
    max_affinity_kt: float
    mean_abs_affinity_kt: float

    def to_dict(self):
        return {
            "n_cycles": self.n_cycles,
            "rms_affinity_kt": self.rms_affinity_kt,
            "max_affinity_kt": self.max_affinity_kt,
            "mean_abs_affinity_kt": self.mean_abs_affinity_kt,
            "method": "schnakenberg_fundamental_cycle_affinities",
        }


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


def _incidence(pairs, n):
    # sparse oriented incidence B (edges x nodes): +1 at i, -1 at j
    from scipy.sparse import csr_matrix

    m = len(pairs)
    rows = np.repeat(np.arange(m), 2)
    cols = pairs.reshape(-1)
    data = np.tile(np.array([1.0, -1.0]), m)
    return csr_matrix((data, (rows, cols)), shape=(m, n))


def _frac_sparse(B, f):
    # gradient (curl-free) projection by sparse least squares: phi = argmin ||B phi - f||^2.
    # scales as O(nnz) per iteration instead of the O(n^3) dense pseudo-inverse.
    from scipy.sparse.linalg import lsqr

    phi = lsqr(B, f, atol=1e-9, btol=1e-9, iter_lim=5000)[0]
    cyclic = f - B.dot(phi)
    return float(cyclic @ cyclic / (float(f @ f) + 1e-12))


def _laplacian_pinv(pairs, n):
    L = np.zeros((n, n))
    i, j = pairs[:, 0], pairs[:, 1]
    np.add.at(L, (i, i), 1.0)
    np.add.at(L, (j, j), 1.0)
    np.add.at(L, (i, j), -1.0)
    np.add.at(L, (j, i), -1.0)
    return np.linalg.pinv(L)


def _frac_dense(pairs, f, n, Lp):
    i, j = pairs[:, 0], pairs[:, 1]
    Btf = np.zeros(n)
    np.add.at(Btf, i, f)
    np.add.at(Btf, j, -f)
    phi = Lp @ Btf
    cyclic = f - (phi[i] - phi[j])
    return float(np.sum(cyclic * cyclic) / (float(np.sum(f * f)) + 1e-12))


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
    f = _edge_flow(pts, vels, pairs)

    try:  # sparse least-squares scales to large graphs; dense pinv is the fallback
        B = _incidence(pairs, n)
        estimate = lambda pr, fv, sub: _frac_sparse(sub, fv)  # noqa: E731
        obs = estimate(pairs, f, B)
    except Exception:
        Lp = _laplacian_pinv(pairs, n)
        B = None
        estimate = lambda pr, fv, sub: _frac_dense(pr, fv, n, Lp)  # noqa: E731
        obs = estimate(pairs, f, None)

    ci = None
    if bootstrap > 1:
        rng = np.random.default_rng(seed)
        m = len(pairs)
        vals = np.empty(bootstrap)
        for b in range(bootstrap):
            sel = rng.integers(0, m, size=m)
            vals[b] = estimate(pairs[sel], f[sel], B[sel] if B is not None else None)
        ci = (float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975)))
    return IrreversibilityReport(obs, ci, len(pairs), n)


def _spanning_tree_potential(pairs, f, n):
    # BFS spanning forest; psi accumulates the flow along tree edges so that
    # psi[a]-psi[b] is the flow along the unique tree path b->a.
    adj = [[] for _ in range(n)]
    for e, (a, b) in enumerate(pairs):
        adj[a].append((b, f[e], e))
        adj[b].append((a, -f[e], e))
    psi = np.full(n, np.nan)
    is_tree = np.zeros(len(pairs), dtype=bool)
    for s in range(n):
        if not np.isnan(psi[s]):
            continue
        psi[s] = 0.0
        dq = deque([s])
        while dq:
            u = dq.popleft()
            for v, fl, e in adj[u]:
                if np.isnan(psi[v]):
                    psi[v] = psi[u] + fl
                    is_tree[e] = True
                    dq.append(v)
    return psi, is_tree


def cycle_affinities(pts, vels, graph, temperature=1.0, velocity_scale=1.0):
    """Thermodynamic affinity of every fundamental cycle, in kT.

    Zero for a conservative (reversible) field, non-zero with real circulation.
    """
    pts = np.asarray(pts, dtype=float)
    vels = np.asarray(vels, dtype=float)
    n = len(pts)
    pairs = _pairs(graph)
    if len(pairs) < 1 or n < 3:
        return CycleAffinityReport(0, 0.0, 0.0, 0.0)
    f = _edge_flow(pts, vels, pairs)
    psi, is_tree = _spanning_tree_potential(pairs, f, n)
    nt = ~is_tree
    if not nt.any():
        return CycleAffinityReport(0, 0.0, 0.0, 0.0)
    a, b = pairs[nt, 0], pairs[nt, 1]
    circ = f[nt] + psi[a] - psi[b]                # circulation around each fundamental cycle
    beta = velocity_scale / max(temperature, 1e-9)
    aff = np.abs(2.0 * beta * circ)               # affinity = 2 * beta * circulation
    return CycleAffinityReport(
        n_cycles=int(nt.sum()),
        rms_affinity_kt=float(np.sqrt(np.mean(aff * aff))),
        max_affinity_kt=float(aff.max()),
        mean_abs_affinity_kt=float(aff.mean()),
    )
