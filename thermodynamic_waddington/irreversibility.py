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


def _cyclic_edge_flow(pts, vels, graph):
    """Edge pairs, edge flow, and its cyclic (non-gradient) component."""
    pts = np.asarray(pts, dtype=float)
    vels = np.asarray(vels, dtype=float)
    n = len(pts)
    pairs = _pairs(graph)
    if len(pairs) == 0 or n < 3:
        return pairs, None, None
    f = _edge_flow(pts, vels, pairs)
    try:
        from scipy.sparse.linalg import lsqr
        B = _incidence(pairs, n)
        cyclic = f - B.dot(lsqr(B, f, atol=1e-9, btol=1e-9, iter_lim=5000)[0])
    except Exception:
        Lp = _laplacian_pinv(pairs, n)
        i, j = pairs[:, 0], pairs[:, 1]
        Btf = np.zeros(n)
        np.add.at(Btf, i, f)
        np.add.at(Btf, j, -f)
        phi = Lp @ Btf
        cyclic = f - (phi[i] - phi[j])
    return pairs, f, cyclic


def _incident(pairs, values, n):
    per = np.zeros(n)
    np.add.at(per, pairs[:, 0], values)
    np.add.at(per, pairs[:, 1], values)
    return per


def cyclic_flow_per_cell(pts, vels, graph):
    """Per-cell irreversibility density: the cyclic (non-gradient) flow energy
    incident on each cell, normalized to sum to 1. Zero for a gradient field;
    concentrates where the velocity field carries local circulation. Useful as an
    `obs` map showing *where* along a lineage irreversibility lives.

    It is an energy, so it also grows with local velocity magnitude; use
    `cyclic_fraction_per_cell` to compare regions with different velocity scales.
    """
    n = len(pts)
    pairs, f, cyclic = _cyclic_edge_flow(pts, vels, graph)
    if f is None:
        return np.zeros(n)
    per = _incident(pairs, cyclic * cyclic, n)
    total = per.sum()
    return per / total if total > 0 else per


def cyclic_fraction_per_cell(pts, vels, graph):
    """Per-cell cyclic fraction in [0, 1]: on the edges incident to each cell, the
    cyclic flow energy over cyclic plus gradient energy. The local analogue of
    `cyclic_irreversibility`; unlike `cyclic_flow_per_cell` it does not grow with
    velocity magnitude, so regions with different sequencing depth or speed
    compare on one scale. (Gradient and cyclic parts are orthogonal only over the
    whole graph, so the local denominator is their energy sum, not the raw flow.)
    """
    n = len(pts)
    pairs, f, cyclic = _cyclic_edge_flow(pts, vels, graph)
    if f is None:
        return np.zeros(n)
    gradient = f - cyclic
    cyc = _incident(pairs, cyclic * cyclic, n)
    total = cyc + _incident(pairs, gradient * gradient, n)
    return np.divide(cyc, total, out=np.zeros(n), where=total > 0)


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


@dataclass
class NoiseMatchedReport:
    """Cyclic fraction read against a null built from the data's own noise.

    `excess` is the observed cyclic fraction minus the noise-matched null mean,
    bias-corrected by running the same procedure on reference datasets that are
    reversible by construction. ~0 means the circulation is what this dataset's
    noise produces on a reversible flow; > 0 means circulation beyond it.
    """

    cyclic_fraction: float
    null_mean: float
    bias: float
    excess: float
    z: float
    draws: int
    references: int
    n_cells: int
    per_cell_excess: np.ndarray | None = None

    def to_dict(self):
        return {
            "cyclic_fraction": self.cyclic_fraction, "null_mean": self.null_mean, "bias": self.bias,
            "excess": self.excess, "z": self.z, "draws": self.draws, "references": self.references,
            "n_cells": self.n_cells, "method": "noise_matched_hodge_cyclic_fraction",
        }


def _neighbor_matrix(graph, n):
    # each cell plus its kNN, padded with self so rows are equal length
    width = 1 + max((len(nb) for nb in graph.neighbors), default=0)
    out = np.empty((n, width), dtype=int)
    for i, nb in enumerate(graph.neighbors):
        row = [i, *nb]
        out[i] = row + [i] * (width - len(row))
    return out


def _reversible_velocity(pts, vels, pairs, nbm, ridge=1e-8):
    """Per-cell velocity of the best potential flow: fit a node potential U to the
    work along every edge, then take its local least-squares gradient."""
    from scipy.sparse.linalg import lsqr

    i, j = pairs[:, 0], pairs[:, 1]
    work = 0.5 * np.einsum("kd,kd->k", vels[i] + vels[j], pts[j] - pts[i])
    U = -lsqr(_incidence(pairs, len(pts)), work, atol=1e-10, btol=1e-10, iter_lim=10000)[0]
    D = pts[nbm[:, 1:]] - pts[:, None, :]                  # n x k x d displacements
    y = U[nbm[:, 1:]] - U[:, None]                          # n x k potential differences
    G = np.einsum("nkd,nke->nde", D, D) + ridge * np.eye(pts.shape[1])
    return np.linalg.solve(G, np.einsum("nkd,nk->nd", D, y)[..., None])[..., 0]


def noise_matched_irreversibility(pts, vels, graph, draws=10, references=5, seed=0, per_cell=False,
                                  systematic=0.01):
    """Cyclic fraction against a null that carries the data's own noise.

    The matched gradient floor (a clean radial field) has no noise, so noisy or
    shallow-sequenced regions read as circulating. Here the null is built from
    the data itself: a reversible velocity (the local gradient of the potential
    that best explains the work along each edge) plus the observed residual
    noise (each cell's velocity minus its neighbourhood mean), re-drawn from
    neighbouring cells with random signs. Running the identical procedure on
    `references` reversible-by-construction datasets estimates the method's own
    bias, which is subtracted.

    Returns a NoiseMatchedReport; z uses the null and reference spread plus a
    `systematic` term (the largest residual bias seen on ground truth) so that
    nearly noise-free data cannot produce an inflated z.
    """
    pts = np.asarray(pts, dtype=float)
    vels = np.asarray(vels, dtype=float)
    n = len(pts)
    pairs = _pairs(graph)
    if len(pairs) == 0 or n < 4:
        return NoiseMatchedReport(0.0, 0.0, 0.0, 0.0, 0.0, draws, references, n)
    B = _incidence(pairs, n)
    nbm = _neighbor_matrix(graph, n)
    rng = np.random.default_rng(seed)

    def frac(v):
        return _frac_sparse(B, _edge_flow(pts, v, pairs))

    def redraw(vrev, resid):
        pick = nbm[np.arange(n), rng.integers(0, nbm.shape[1], size=n)]
        return vrev + rng.choice([-1.0, 1.0], size=n)[:, None] * resid[pick]

    def null_of(v, k, local):
        vrev = _reversible_velocity(pts, v, pairs, nbm)
        resid = v - v[nbm].mean(axis=1)
        fields = [redraw(vrev, resid) for _ in range(k)]
        glob = np.array([frac(f) for f in fields])
        loc = np.mean([cyclic_fraction_per_cell(pts, f, graph) for f in fields], axis=0) if local else None
        return glob, loc, vrev, resid

    obs = frac(vels)
    nulls, null_loc, vrev, resid = null_of(vels, draws, per_cell)
    bias, bias_loc = [], []
    for _ in range(references):
        ref = redraw(vrev, resid)
        rn, rl, _, _ = null_of(ref, draws, per_cell)
        bias.append(frac(ref) - rn.mean())
        if per_cell:
            bias_loc.append(cyclic_fraction_per_cell(pts, ref, graph) - rl)
    bias = np.asarray(bias)
    excess = obs - nulls.mean() - bias.mean()
    spread = float(np.sqrt(nulls.var() + bias.var() + systematic ** 2))
    local = None
    if per_cell:
        local = cyclic_fraction_per_cell(pts, vels, graph) - null_loc - np.mean(bias_loc, axis=0)
    return NoiseMatchedReport(
        cyclic_fraction=float(obs), null_mean=float(nulls.mean()), bias=float(bias.mean()),
        excess=float(excess), z=float(excess / spread), draws=draws, references=references, n_cells=n,
        per_cell_excess=local,
    )
