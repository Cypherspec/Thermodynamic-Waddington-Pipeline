"""Normalizes raw count data before feeding it into the pipeline.

Standard two-step (Luecken & Theis 2019): per-cell size normalization,
then log1p. Without this, raw UMI counts overflow the Arrhenius exp() in
edge_work and produce garbage results -- found empirically on real pancreas
data, not caught by any synthetic test.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, shape


@dataclass(frozen=True)
class NormalizationReport:
    """What was done to the data, so transforms are never hidden from the caller."""
    method: str
    target_sum: float
    n_cells: int
    n_genes: int
    size_factors: list[float]
    pre_normalization_total_range: tuple[float, float]
    post_normalization_total_range: tuple[float, float]
    log_transformed: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "method": self.method,
            "target_sum": self.target_sum,
            "n_cells": self.n_cells,
            "n_genes": self.n_genes,
            "size_factor_range": [min(self.size_factors), max(self.size_factors)] if self.size_factors else [0.0, 0.0],
            "pre_total_range": list(self.pre_normalization_total_range),
            "post_total_range": list(self.post_normalization_total_range),
            "log_transformed": self.log_transformed,
            "warning": (
                "Raw UMI counts should be normalized before fit_landscape. "
                "Input is assumed to be roughly unit-scale. "
                "Unnormalized input can overflow edge_work's exp() and produce "
                "meaningless entropy production estimates."
            ),
        }


def _row_totals(expr: Sequence[Sequence[float]]) -> list[float]:
    return [sum(row) for row in expr]


def _median(vals: Sequence[float]) -> float:
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return 0.0
    mid = n // 2
    return s[mid] if n % 2 else 0.5 * (s[mid - 1] + s[mid])


def normalize_total(
    expr: Sequence[Sequence[float]],
    target_sum: float | None = None,
    floor: float = 1e-9,
) -> tuple[list[list[float]], list[float]]:
    """Divide each cell by its total, rescale to target_sum.

    Returns (normalized_expression, size_factors). Size factors are reused
    by velocity transforms to keep the two on compatible scales.
    """
    totals = _row_totals(expr)
    if target_sum is None:
        pos = [t for t in totals if t > 0]
        target_sum = _median(pos) if pos else 1.0
    sf = [max(t, floor) / target_sum for t in totals]
    norm = [[v / max(f, floor) for v in row] for row, f in zip(expr, sf)]
    return norm, sf


def log1p(expr: Sequence[Sequence[float]]) -> list[list[float]]:
    """log(1+x) elementwise. Only valid for non-negative input (counts)."""
    return [[math.log1p(max(0.0, v)) for v in row] for row in expr]


def normalize_expression(
    expr: Sequence[Sequence[float]],
    target_sum: float | None = None,
    apply_log1p: bool = True,
) -> tuple[list[list[float]], NormalizationReport]:
    """Standard pipeline: size normalization then log1p."""
    n_cells, n_genes = shape(expr)
    pre = _row_totals(expr)
    norm, sf = normalize_total(expr, target_sum=target_sum)
    post = _row_totals(norm)
    if apply_log1p:
        norm = log1p(norm)
    ts = target_sum if target_sum is not None else (_median([t for t in pre if t > 0]) if any(t > 0 for t in pre) else 1.0)
    report = NormalizationReport(
        method="total_count" + ("_log1p" if apply_log1p else ""),
        target_sum=ts,
        n_cells=n_cells,
        n_genes=n_genes,
        size_factors=sf,
        pre_normalization_total_range=(min(pre) if pre else 0.0, max(pre) if pre else 0.0),
        post_normalization_total_range=(min(post) if post else 0.0, max(post) if post else 0.0),
        log_transformed=apply_log1p,
    )
    return norm, report


def scale_velocity_to_expression(vel: Sequence[Sequence[float]], sf: Sequence[float]) -> list[list[float]]:
    """Apply expression size factors to velocity without log1p (velocity is signed).
    Quick approximation -- prefer velocity_in_transformed_space for real data."""
    return [[v / max(f, 1e-9) for v in row] for row, f in zip(vel, sf)]


def velocity_in_transformed_space(
    expr: Sequence[Sequence[float]],
    vel: Sequence[Sequence[float]],
    target_sum: float | None = None,
) -> tuple[list[list[float]], list[list[float]], NormalizationReport]:
    """Puts velocity in log1p space via finite difference: log1p(x+v) - log1p(x).

    This is the scVelo approach (La Manno 2018, Bergen 2020). It avoids the
    scale mismatch between raw velocity and log-normalized expression that
    causes entropy production estimates to blow up on real count data.
    """
    norm_expr, report = normalize_expression(expr, target_sum=target_sum, apply_log1p=True)
    # extrapolate one step using the same size factors -- not recomputed from
    # the extrapolated state, to avoid velocity magnitude affecting normalization
    extrap_raw = [[e + v for e, v in zip(er, vr)] for er, vr in zip(expr, vel)]
    extrap_scaled = [[v / max(f, 1e-9) for v in row] for row, f in zip(extrap_raw, report.size_factors)]
    extrap_log = log1p(extrap_scaled)
    vel_t = [[e - b for e, b in zip(er, br)] for er, br in zip(extrap_log, norm_expr)]
    return norm_expr, vel_t, report
