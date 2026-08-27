from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import dot, norm, pairwise_squared_distances, safe_log, sub, variance


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


def build_knn(pts: Sequence[Sequence[float]], k: int) -> NeighborGraph:
    sq = pairwise_squared_distances(pts)
    nbrs: list[list[int]] = []
    dists: list[list[float]] = []
    for i, row in enumerate(sq):
        order = sorted((j for j in range(len(row)) if j != i), key=lambda j: row[j])[:k]
        nbrs.append(order)
        dists.append([math.sqrt(row[j]) for j in order])
    edges = [Edge(i, j, math.sqrt(sq[i][j]), 0.0) for i, js in enumerate(nbrs) for j in js]
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
