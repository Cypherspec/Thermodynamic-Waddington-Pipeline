from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, variance


@dataclass(frozen=True)
class VelocityDiagnostics:
    velocity_norm_mean: float
    velocity_norm_p95: float
    expression_velocity_correlation: float
    phase_coherence: float

    def as_dict(self) -> dict:
        return {"velocity_norm_mean": self.velocity_norm_mean, "velocity_norm_p95": self.velocity_norm_p95, "expression_velocity_correlation": self.expression_velocity_correlation, "phase_coherence": self.phase_coherence}


def velocity_norms(velocity: Sequence[Sequence[float]]) -> list[float]:
    return [math.sqrt(sum(value * value for value in row)) for row in velocity]


def correlation(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right) or not left:
        return 0.0
    left_mean = mean(left)
    right_mean = mean(right)
    denominator = math.sqrt(variance(left) * variance(right))
    return sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right)) / max(1e-12, len(left) * denominator)


def diagnose(expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]]) -> VelocityDiagnostics:
    norms = velocity_norms(velocity)
    expression_signal = [mean(row) for row in expression]
    velocity_signal = [mean(row) for row in velocity]
    ordered = sorted(norms)
    p95 = ordered[min(len(ordered) - 1, int(0.95 * len(ordered)))] if ordered else 0.0
    coherence = abs(correlation(expression_signal, velocity_signal))
    return VelocityDiagnostics(mean(norms), p95, correlation(expression_signal, velocity_signal), coherence)
