from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .arrays import matrix_multiply, matvec, transpose


@dataclass(frozen=True)
class KineticSummary:
    stationary_distribution: list[float]
    mixing_time: int
    entropy_rate: float
    spectral_gap: float


def row_normalize(weights: Sequence[Sequence[float]], self_loop: float = 1e-6) -> list[list[float]]:
    matrix = [list(row) for row in weights]
    for i, row in enumerate(matrix):
        if not row:
            continue
        row[i] += self_loop
        total = sum(max(0.0, value) for value in row)
        if total:
            matrix[i] = [max(0.0, value) / total for value in row]
    return matrix


def propagate(matrix: Sequence[Sequence[float]], initial: Sequence[float], steps: int) -> list[float]:
    vector = list(initial)
    for _ in range(max(0, steps)):
        vector = matvec(transpose(matrix), vector)
        total = sum(vector)
        if total:
            vector = [value / total for value in vector]
    return vector


def entropy_rate(matrix: Sequence[Sequence[float]], stationary: Sequence[float]) -> float:
    total = 0.0
    for weight, row in zip(stationary, matrix):
        for probability in row:
            if probability > 0:
                import math
                total -= weight * probability * math.log(probability)
    return total


def summarize(matrix: Sequence[Sequence[float]], steps: int = 250) -> KineticSummary:
    n = len(matrix)
    if not n:
        return KineticSummary([], 0, 0.0, 0.0)
    transition = row_normalize(matrix)
    stationary = [1.0 / n] * n
    previous = stationary
    mixing = steps
    for step in range(1, steps + 1):
        current = propagate(transition, previous, 1)
        distance = sum(abs(a - b) for a, b in zip(current, previous))
        if distance < 1e-8:
            mixing = step
            stationary = current
            break
        previous = current
        stationary = current
    diagonal = sum(transition[i][i] for i in range(n)) / n
    gap = max(0.0, 1.0 - diagonal)
    return KineticSummary(stationary, mixing, entropy_rate(transition, stationary), gap)
