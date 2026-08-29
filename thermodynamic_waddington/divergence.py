"""Local velocity-field divergence as an independent EP estimator.

div(v)(x) = sum_i d/dx_i v_i(x)

For a Fokker-Planck system at steady state, the EP rate density is
proportional to |v - D*grad(log p)|^2 / D (Tome & de Oliveira 2012).
The divergence gives a simpler proxy: div(v) < 0 means local convergence
(attractor-like), div(v) > 0 means local divergence (source-like). Net
positive mean divergence across a trajectory is inconsistent with a
detailed-balance equilibrium -- a necessary (not sufficient) condition for
irreversibility.

This is completely independent of the kNN flux formula in entropy_production.py.
Agreement between the two is stronger evidence than either alone.

References:
    Tome & de Oliveira (2012) Phys Rev E 86:021136
    Fang et al. (2019) Phys Rev Lett 122:170602
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .graph import NeighborGraph
from .arrays import mean


@dataclass
class DivergenceResult:
    per_cell: list[float]          # div(v) estimated at each cell
    mean_divergence: float         # mean across all cells
    std_divergence: float          # std across all cells
    fraction_converging: float     # fraction with div(v) < 0
    fraction_diverging: float      # fraction with div(v) > 0
    weighted_mean: float           # density-weighted mean
    sign_consistency: float        # |mean| / std, higher = more consistent direction
    audit: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        return {
            "mean_divergence": self.mean_divergence,
            "std_divergence": self.std_divergence,
            "fraction_converging": self.fraction_converging,
            "fraction_diverging": self.fraction_diverging,
            "weighted_mean": self.weighted_mean,
            "sign_consistency": self.sign_consistency,
            "audit": self.audit,
        }


def _local_divergence(
    i: int,
    pts: Sequence[Sequence[float]],
    vels: Sequence[Sequence[float]],
    nbrs: list[int],
    dists: list[float],
    bw: float,
) -> float:
    """Kernel-weighted finite-difference divergence at cell i.

    For each neighbor j, estimate d/dx_k v_k ~ (v_j[k] - v_i[k]) / dx[k]
    for each dimension k, weighted by a Gaussian kernel.
    Sum over k gives the local divergence.
    """
    if not nbrs:
        return 0.0

    n_dim = len(pts[i])
    denom = 2.0 * bw * bw
    total_weight = 0.0
    div_sum = 0.0

    for j, d in zip(nbrs, dists):
        w = math.exp(-(d * d) / denom)
        # displacement vector from i to j
        dx = [pts[j][k] - pts[i][k] for k in range(n_dim)]
        # velocity difference
        dv = [vels[j][k] - vels[i][k] for k in range(n_dim)]
        # dot(dv, dx) / dot(dx, dx) approximates div(v) along this direction
        dx2 = sum(x * x for x in dx)
        if dx2 > 1e-12:
            div_sum += w * sum(dv[k] * dx[k] for k in range(n_dim)) / dx2
        total_weight += w

    if total_weight < 1e-12:
        return 0.0
    return div_sum / total_weight


def estimate_divergence(
    pts: Sequence[Sequence[float]],
    vels: Sequence[Sequence[float]],
    graph: NeighborGraph,
    dens: Sequence[float],
    bandwidth: float = 1.0,
) -> DivergenceResult:
    """Estimate local velocity-field divergence at every cell.

    Uses kernel-weighted finite differences over each cell's kNN neighborhood.
    bandwidth should be on the order of the median neighbor distance.
    """
    n = len(pts)
    per_cell = [
        _local_divergence(i, pts, vels, graph.neighbors[i], graph.distances[i], bandwidth)
        for i in range(n)
    ]

    mu = mean(per_cell) if per_cell else 0.0
    variance = mean((v - mu) ** 2 for v in per_cell) if per_cell else 0.0
    sigma = math.sqrt(variance)

    n_conv = sum(1 for v in per_cell if v < 0)
    n_div  = sum(1 for v in per_cell if v > 0)

    # density-weighted mean -- cells in dense regions count more
    total_d = sum(dens) or 1.0
    w_mean = sum(per_cell[i] * dens[i] / total_d for i in range(n))

    consistency = abs(mu) / max(sigma, 1e-12)

    return DivergenceResult(
        per_cell=per_cell,
        mean_divergence=mu,
        std_divergence=sigma,
        fraction_converging=n_conv / max(n, 1),
        fraction_diverging=n_div / max(n, 1),
        weighted_mean=w_mean,
        sign_consistency=consistency,
        audit={
            "n_cells": n,
            "bandwidth": bandwidth,
            "method": "kernel_weighted_finite_difference",
            "interpretation": (
                "mean_divergence > 0 means net source behavior (inconsistent with "
                "detailed balance). mean_divergence < 0 means net convergence "
                "(attractor-dominated). sign_consistency > 1 means the mean direction "
                "is larger than cell-to-cell variation."
            ),
        },
    )


def cross_validate_divergence_vs_ep(
    div_result: DivergenceResult,
    ep_rate: float,
) -> dict[str, object]:
    """Check whether div and EP estimators agree on sign and rough magnitude.

    They measure related but distinct things:
    - div(v) is local and doesn't need a graph structure
    - EP rate integrates pairwise fluxes globally
    Both should be positive if the system is nonequilibrium, zero at equilibrium.
    """
    # sign agreement: both methods should say positive EP if nonequilibrium
    div_sign = 1 if div_result.weighted_mean > 0 else (-1 if div_result.weighted_mean < 0 else 0)
    ep_sign  = 1 if ep_rate > 0 else (-1 if ep_rate < 0 else 0)
    signs_agree = div_sign == ep_sign

    return {
        "div_weighted_mean": div_result.weighted_mean,
        "div_sign_consistency": div_result.sign_consistency,
        "ep_rate": ep_rate,
        "signs_agree": signs_agree,
        "fraction_converging": div_result.fraction_converging,
        "fraction_diverging": div_result.fraction_diverging,
        "note": (
            "Sign agreement is the meaningful check. Magnitude comparison is not "
            "valid because div(v) and EP have different units and normalizations."
        ),
    }
