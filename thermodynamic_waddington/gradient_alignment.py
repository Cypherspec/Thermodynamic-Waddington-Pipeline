"""Velocity-landscape gradient alignment.

Tests whether the RNA velocity field points "downhill" on the fitted
free-energy landscape -- i.e., whether the landscape is actually predictive
of the direction cells move, not just a static density description.

The gradient of the free energy at cell i is estimated from the landscape
values at neighboring cells: grad F(i) ~ sum_j w_ij * (F_j - F_i) * disp_ij
where disp_ij is the unit vector from i to j and w_ij is a kernel weight.

Alignment score: cos(angle between velocity and -grad F). Score = 1 means
velocity points exactly downhill (landscape fully predictive). Score = 0
means orthogonal (landscape says nothing about direction). Score = -1 means
velocity points uphill (landscape is anti-predictive, sign error somewhere).

If the landscape is a true Lyapunov function for the velocity dynamics,
we expect mean alignment > 0 across cells. This is a necessary condition
for the landscape to have predictive validity -- a weaker claim than the
thermodynamic one, and directly testable from the same data.

Reference:
    Ao (2004) J Phys A 37:L25  (decomposition of dynamics into gradient + curl)
    Fang et al. (2019) Phys Rev Lett 122:060601
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import dot, norm, mean
from .graph import NeighborGraph


@dataclass
class GradientAlignmentResult:
    per_cell: list[float]          # cos(angle) between v_i and -grad F_i
    mean_alignment: float          # mean across cells (> 0 = downhill on average)
    std_alignment: float
    fraction_downhill: float       # fraction with cos > 0
    fraction_uphill: float         # fraction with cos < 0
    weighted_mean: float           # energy-weighted: high-energy cells count more
    gradient_magnitudes: list[float]  # ||grad F|| per cell
    velocity_magnitudes: list[float]
    audit: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        return {
            "mean_alignment": self.mean_alignment,
            "std_alignment": self.std_alignment,
            "fraction_downhill": self.fraction_downhill,
            "fraction_uphill": self.fraction_uphill,
            "weighted_mean_alignment": self.weighted_mean,
            "mean_gradient_magnitude": mean(self.gradient_magnitudes),
            "mean_velocity_magnitude": mean(self.velocity_magnitudes),
            "audit": self.audit,
        }


def _landscape_gradient(
    i: int,
    pts: Sequence[Sequence[float]],
    energies: Sequence[float],
    nbrs: list[int],
    dists: list[float],
    bw: float,
) -> list[float]:
    """Kernel-weighted gradient of F at cell i using neighbor energy differences.

    For each neighbor j: contribution = w_ij * (F_j - F_i) * unit_vec(i->j)
    Sum and normalize by total weight to get the gradient vector.
    """
    n_dim = len(pts[i])
    if not nbrs:
        return [0.0] * n_dim

    denom = 2.0 * bw * bw
    grad = [0.0] * n_dim
    total_w = 0.0

    for j, d in zip(nbrs, dists):
        if d < 1e-12:
            continue
        w = math.exp(-(d * d) / denom)
        dE = energies[j] - energies[i]
        # unit displacement from i to j
        disp = [(pts[j][k] - pts[i][k]) / d for k in range(n_dim)]
        for k in range(n_dim):
            grad[k] += w * dE * disp[k]
        total_w += w

    if total_w > 1e-12:
        grad = [g / total_w for g in grad]
    return grad


def estimate_gradient_alignment(
    pts: Sequence[Sequence[float]],
    vels: Sequence[Sequence[float]],
    energies: Sequence[float],
    graph: NeighborGraph,
    bandwidth: float = 1.0,
) -> GradientAlignmentResult:
    """Compute per-cell alignment between velocity and negative landscape gradient.

    alignment_i = cos(v_i, -grad F_i)
                = -dot(v_i, grad F_i) / (||v_i|| * ||grad F_i||)

    Positive alignment means velocity points downhill -- the landscape
    predicts the direction of transcriptional change for that cell.
    """
    n = len(pts)
    alignments: list[float] = []
    grad_mags: list[float] = []
    vel_mags: list[float] = []

    for i in range(n):
        grad = _landscape_gradient(i, pts, energies, graph.neighbors[i],
                                   graph.distances[i], bandwidth)
        v = vels[i]
        gn = norm(grad)
        vn = norm(v)
        grad_mags.append(gn)
        vel_mags.append(vn)

        if gn < 1e-12 or vn < 1e-12:
            alignments.append(0.0)
            continue
        # negative gradient = downhill direction; alignment with velocity
        neg_grad = [-g for g in grad]
        alignments.append(dot(v, neg_grad) / (vn * gn))

    mu = mean(alignments)
    variance = mean((a - mu) ** 2 for a in alignments)
    sigma = math.sqrt(variance)

    n_down = sum(1 for a in alignments if a > 0)
    n_up   = sum(1 for a in alignments if a < 0)

    # energy-weighted mean: weight cells by their free energy
    # (high-energy progenitor cells matter more for testing predictivity)
    finite_e = [max(0.0, e) for e in energies]
    total_e  = sum(finite_e) or 1.0
    w_mean   = sum(alignments[i] * finite_e[i] / total_e for i in range(n))

    return GradientAlignmentResult(
        per_cell=alignments,
        mean_alignment=mu,
        std_alignment=sigma,
        fraction_downhill=n_down / max(n, 1),
        fraction_uphill=n_up / max(n, 1),
        weighted_mean=w_mean,
        gradient_magnitudes=grad_mags,
        velocity_magnitudes=vel_mags,
        audit={
            "n_cells": n,
            "bandwidth": bandwidth,
            "n_zero_gradient": sum(1 for g in grad_mags if g < 1e-12),
            "n_zero_velocity": sum(1 for v in vel_mags if v < 1e-12),
            "interpretation": (
                "mean_alignment > 0: velocity points downhill on average "
                "(landscape has predictive validity). "
                "mean_alignment ~ 0: landscape and velocity are orthogonal "
                "(landscape is static density, not predictive of dynamics). "
                "mean_alignment < 0: velocity points uphill (sign inconsistency "
                "between landscape construction and velocity direction)."
            ),
        },
    )


def permutation_test_alignment(
    pts: Sequence[Sequence[float]],
    vels: Sequence[Sequence[float]],
    energies: Sequence[float],
    graph: NeighborGraph,
    bandwidth: float = 1.0,
    n_permutations: int = 200,
    seed: int = 17,
) -> dict[str, object]:
    """Permutation test for gradient alignment.

    Null: shuffle cell positions while keeping velocity-energy pairs intact.
    If the observed mean alignment exceeds the null distribution, the
    landscape is genuinely predictive of velocity direction, not a spurious
    geometric correlation.
    """
    import random
    rng = random.Random(seed)

    obs = estimate_gradient_alignment(pts, vels, energies, graph, bandwidth)
    observed_stat = obs.mean_alignment

    idx = list(range(len(pts)))
    null_stats: list[float] = []
    for _ in range(n_permutations):
        shuffled = idx[:]
        rng.shuffle(shuffled)
        perm_energies = [energies[shuffled[i]] for i in range(len(pts))]
        null_r = estimate_gradient_alignment(pts, vels, perm_energies, graph, bandwidth)
        null_stats.append(null_r.mean_alignment)

    exceed = sum(1 for v in null_stats if v >= observed_stat)
    p_val  = (exceed + 1) / (len(null_stats) + 1)

    return {
        "observed_alignment": observed_stat,
        "null_mean": sum(null_stats) / max(len(null_stats), 1),
        "null_std": math.sqrt(sum((v - sum(null_stats)/max(len(null_stats),1))**2
                                  for v in null_stats) / max(len(null_stats), 1)),
        "p_value": p_val,
        "n_permutations": n_permutations,
        "significant": p_val < 0.05,
        "interpretation": (
            "p < 0.05: landscape gradient significantly predicts velocity direction. "
            "This is evidence the free-energy landscape captures real dynamics, "
            "not just a smooth summary of static cell density."
        ),
    }
