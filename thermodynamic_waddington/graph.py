from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .arrays import dot, norm, safe_log, sub, variance


@dataclass(frozen=True)
class Edge:
    source: int
    target: int
    distance: float
    alignment: float
    work: float = 0.0
    action: str = "forward"


@dataclass
class NeighborGraph:
    neighbors: list[list[int]]
    distances: list[list[float]]
    edges: list[Edge]


def _knn(x: np.ndarray, k: int) -> tuple[list[list[int]], list[list[float]]]:
    # exact k nearest neighbors, ordered by (distance, index) to match the old
    # brute-force sort exactly. Uses a KD-tree when scipy is present so the cost
    # is n log n and no n-by-n matrix is built; falls back to numpy otherwise.
    n = x.shape[0]
    k = max(1, min(k, n - 1))
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(x)
        m = min(n, k + 1)
        dd, ii = tree.query(x, k=m)
        if m == 1:
            ii = ii.reshape(-1, 1)
            dd = dd.reshape(-1, 1)
        nbrs: list[list[int]] = []
        dists: list[list[float]] = []
        for i in range(n):
            idx = ii[i]
            dst = dd[i]
            keep = idx != i
            idx = idx[keep][:k]
            dst = dst[keep][:k]
            order = np.lexsort((idx, dst))
            nbrs.append(idx[order].tolist())
            dists.append(dst[order].tolist())
        return nbrs, dists
    except Exception:
        sq_norm = np.einsum("ij,ij->i", x, x)
        d2 = sq_norm[:, None] + sq_norm[None, :] - 2.0 * (x @ x.T)
        np.fill_diagonal(d2, np.inf)
        np.clip(d2, 0.0, None, out=d2)
        part = np.argpartition(d2, k - 1, axis=1)[:, :k]
        nbrs = []
        dists = []
        for i in range(n):
            cand = part[i]
            order = np.lexsort((cand, d2[i, cand]))
            ci = cand[order]
            nbrs.append(ci.tolist())
            dists.append(np.sqrt(d2[i, ci]).tolist())
        return nbrs, dists


def build_knn(pts: Sequence[Sequence[float]], k: int) -> NeighborGraph:
    x = np.asarray(pts, dtype=float)
    n = x.shape[0]
    nbrs, dists = _knn(x, k)
    edges = [Edge(i, nbrs[i][t], dists[i][t], 0.0) for i in range(n) for t in range(len(nbrs[i]))]
    return NeighborGraph(nbrs, dists, edges)


def velocity_alignment(pt: Sequence[float], vel: Sequence[float], tgt: Sequence[float]) -> float:
    # cosine similarity between velocity and displacement to neighbor
    disp = sub(tgt, pt)
    d = norm(disp)
    s = norm(vel)
    if d == 0 or s == 0:
        return 0.0
    return dot(vel, disp) / (s * d)


def estimate_local_diffusion(graph: NeighborGraph, vels: Sequence[Sequence[float]], floor: float) -> list[float]:
    out: list[float] = []
    for i, nbr in enumerate(graph.neighbors):
        if not nbr:
            out.append(floor)
            continue
        # variance across neighbor velocities per feature, averaged
        comp_vars = [variance(vels[j][f] for j in nbr) for f in range(len(vels[i]))]
        out.append(max(floor, sum(comp_vars) / max(1, len(comp_vars))))
    return out


def local_density(pts: Sequence[Sequence[float]], graph: NeighborGraph, bw: float) -> list[float]:
    # Gaussian kernel density using precomputed neighbor distances
    denom = 2.0 * bw * bw
    dens: list[float] = []
    for i, nbr in enumerate(graph.neighbors):
        ksum = sum(math.exp(-(d * d) / denom) for d in graph.distances[i])
        dens.append((1.0 + ksum) / (1.0 + len(nbr)))
    return dens


def directed_edges(pts: Sequence[Sequence[float]], vels: Sequence[Sequence[float]],
                   graph: NeighborGraph, thresh: float) -> list[Edge]:
    out: list[Edge] = []
    for i, nbr in enumerate(graph.neighbors):
        for j, d in zip(nbr, graph.distances[i]):
            a = velocity_alignment(pts[i], vels[i], pts[j])
            if a >= thresh:
                out.append(Edge(i, j, d, a))
    return out


def edge_work(e: Edge, pts: Sequence[Sequence[float]], vels: Sequence[Sequence[float]],
              dens: Sequence[float], diff: Sequence[float],
              temp: float, vscale: float) -> Edge:
    src, tgt = e.source, e.target
    disp = sub(pts[tgt], pts[src])
    d = max(diff[src], 1e-12)
    # three terms: density ratio, velocity drift, noise
    density_term = -temp * math.log(max(dens[tgt], 1e-12) / max(dens[src], 1e-12))
    drift_term = -vscale * dot(vels[src], disp)
    noise_term = 0.5 * dot(disp, disp) / d
    w = (density_term + drift_term + noise_term) / max(temp, 1e-12)
    return Edge(src, tgt, e.distance, e.alignment, w, "forward")


def annotate_edges(pts, vels, graph, dens, diff, temp, vscale, thresh) -> list[Edge]:
    candidates = directed_edges(pts, vels, graph, thresh)
    return [edge_work(e, pts, vels, dens, diff, temp, vscale) for e in candidates]


def adjacency(edges: Sequence[Edge], n: int) -> list[list[Edge]]:
    adj: list[list[Edge]] = [[] for _ in range(n)]
    for e in edges:
        adj[e.source].append(e)
    return adj
