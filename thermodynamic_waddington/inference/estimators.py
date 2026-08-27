from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Sequence

from ..arrays import mean, safe_exp
from ..graph import Edge


@dataclass(frozen=True)
class EstimatorComparison:
    method: str
    estimates: tuple[float, ...]
    bias_proxy: float
    effective_sample_size: float
    variance: float

    def to_dict(self) -> dict[str, object]:
        return {"method": self.method, "estimates": list(self.estimates), "bias_proxy": self.bias_proxy, "effective_sample_size": self.effective_sample_size, "variance": self.variance}


def naive_drift_integral(edges: Sequence[Edge], n: int, reference: int) -> list[float]:
    energy = [float("inf")] * n
    energy[reference] = 0.0
    for _ in range(n):
        changed = False
        for edge in edges:
            candidate = energy[edge.source] + edge.work
            if candidate < energy[edge.target]:
                energy[edge.target] = candidate
                changed = True
        if not changed:
            break
    finite = [value for value in energy if math.isfinite(value)]
    fill = max(finite) if finite else 0.0
    return [value if math.isfinite(value) else fill for value in energy]


def density_potential(densities: Sequence[float], temperature: float) -> list[float]:
    scale = max(temperature, 1e-12)
    return [-scale * math.log(max(density, 1e-300)) for density in densities]


def weighted_work_profile(edges: Sequence[Edge], temperature: float) -> EstimatorComparison:
    if not edges:
        return EstimatorComparison("weighted_work", (), 0.0, 0.0, 0.0)
    scale = max(temperature, 1e-12)
    values = [edge.work for edge in edges]
    weights = [safe_exp(-value / scale) for value in values]
    total = sum(weights) or 1.0
    estimate = tuple(-scale * math.log(max(weight / total, 1e-300)) for weight in weights)
    effective = total * total / max(1e-300, sum(weight * weight for weight in weights))
    center = mean(estimate)
    variance = mean([(value - center) ** 2 for value in estimate])
    return EstimatorComparison("weighted_work", estimate, 0.0, effective, variance)


def bootstrap_paths(edges: Sequence[Edge], temperature: float, replicates: int = 128, seed: int = 11) -> list[EstimatorComparison]:
    rng = random.Random(seed)
    result = []
    for _ in range(max(1, replicates)):
        sampled = [edges[rng.randrange(len(edges))] for _ in edges] if edges else []
        result.append(weighted_work_profile(sampled, temperature))
    return result


def estimator_disagreement(comparisons: Sequence[EstimatorComparison]) -> float:
    centers = [mean(comparison.estimates) for comparison in comparisons if comparison.estimates]
    return max(centers) - min(centers) if centers else 0.0
