from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence


@dataclass(frozen=True)
class SensitivityPoint:
    parameter: str
    value: float
    mean_energy: float
    energy_range: float
    attractor_count: int


def sweep(parameter: str, values: Sequence[float], fit_factory: Callable[[float], object]) -> list[SensitivityPoint]:
    points: list[SensitivityPoint] = []
    for value in values:
        result = fit_factory(float(value))
        energies = list(getattr(result, "energies"))
        points.append(SensitivityPoint(parameter, float(value), sum(energies) / len(energies) if energies else 0.0, max(energies) - min(energies) if energies else 0.0, len(getattr(result, "attractors", []))))
    return points


def normalized_instability(points: Sequence[SensitivityPoint]) -> float:
    if not points:
        return 0.0
    baseline = max(1e-12, abs(points[0].energy_range))
    return sum(abs(point.energy_range - points[0].energy_range) for point in points[1:]) / baseline
