from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, percentile
from .graph import Edge


@dataclass(frozen=True)
class SpectralSummary:
    stationary_distribution: list[float]
    mixing_proxy: float
    spectral_gap_proxy: float
    directed_asymmetry: float

    def to_dict(self) -> dict[str, object]:
        return {"stationary_distribution": self.stationary_distribution, "mixing_proxy": self.mixing_proxy, "spectral_gap_proxy": self.spectral_gap_proxy, "directed_asymmetry": self.directed_asymmetry}


def transition_matrix(edges: Sequence[Edge], n: int, temperature: float = 1.0) -> list[list[float]]:
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for edge in edges:
        matrix[edge.source][edge.target] += math.exp(-max(0.0, edge.distance) / max(temperature, 1e-9))
    for i, row in enumerate(matrix):
        total = sum(row)
        if total:
            matrix[i] = [value / total for value in row]
        else:
            matrix[i][i] = 1.0
    return matrix


def stationary_distribution(matrix: Sequence[Sequence[float]], iterations: int = 256) -> list[float]:
    n = len(matrix)
    if not n:
        return []
    distribution = [1.0 / n for _ in range(n)]
    for _ in range(iterations):
        next_distribution = [sum(distribution[i] * matrix[i][j] for i in range(n)) for j in range(n)]
        total = sum(next_distribution)
        distribution = [value / total for value in next_distribution] if total else distribution
    return distribution


def summarize_spectrum(edges: Sequence[Edge], n: int, temperature: float = 1.0) -> SpectralSummary:
    matrix = transition_matrix(edges, n, temperature)
    stationary = stationary_distribution(matrix)
    self_transition = mean([matrix[i][i] for i in range(n)]) if n else 0.0
    asymmetry = mean([abs(edge.alignment) for edge in edges]) if edges else 0.0
    gap = max(0.0, 1.0 - self_transition)
    return SpectralSummary(stationary, 1.0 / max(gap, 1e-9), gap, asymmetry)


def spectral_summary(edges: Sequence[Edge], n: int | Sequence[float], temperature: float = 1.0) -> dict[str, object]:
    if not isinstance(n, int):
        n = len(n)
    return summarize_spectrum(edges, n, temperature).to_dict()
