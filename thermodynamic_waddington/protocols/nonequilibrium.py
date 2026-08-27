from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from ..arrays import dot, mean, safe_exp


@dataclass(frozen=True)
class ProtocolSegment:
    source: int
    target: int
    displacement: tuple[float, ...]
    velocity: tuple[float, ...]
    duration: float
    diffusion: float

    @property
    def mechanical_work(self) -> float:
        return -sum(a * b for a, b in zip(self.velocity, self.displacement)) * self.duration

    @property
    def dissipated_work(self) -> float:
        return max(0.0, self.mechanical_work - self.diffusive_scale)

    @property
    def diffusive_scale(self) -> float:
        return max(self.diffusion * self.duration, 1e-12)

    def to_dict(self) -> dict[str, object]:
        return {"source": self.source, "target": self.target, "displacement": list(self.displacement), "velocity": list(self.velocity), "duration": self.duration, "diffusion": self.diffusion, "mechanical_work": self.mechanical_work, "dissipated_work": self.dissipated_work}


def crooks_log_ratio(forward_work: float, reverse_work: float, temperature: float) -> float:
    return (forward_work - reverse_work) / max(temperature, 1e-12)


def work_histogram(work: Sequence[float], bins: int = 24) -> dict[str, list[float]]:
    if not work:
        return {"edges": [], "counts": []}
    lo, hi = min(work), max(work)
    if hi <= lo:
        return {"edges": [lo, hi], "counts": [len(work)]}
    width = (hi - lo) / bins
    counts = [0] * bins
    for value in work:
        counts[min(bins - 1, int((value - lo) / width))] += 1
    return {"edges": [lo + index * width for index in range(bins + 1)], "counts": counts}


def jarzynski_cumulants(work: Sequence[float], temperature: float) -> dict[str, float]:
    if not work:
        return {"free_energy": 0.0, "mean_work": 0.0, "variance_work": 0.0, "dissipation": 0.0}
    scale = max(temperature, 1e-12)
    shifted = [-value / scale for value in work]
    pivot = max(shifted)
    log_average = pivot + math.log(mean([safe_exp(value - pivot) for value in shifted]))
    average = mean(work)
    variance = mean([(value - average) ** 2 for value in work])
    free_energy = -scale * log_average
    return {"free_energy": free_energy, "mean_work": average, "variance_work": variance, "dissipation": average - free_energy}
