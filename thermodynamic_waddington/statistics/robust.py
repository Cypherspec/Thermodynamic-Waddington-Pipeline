from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Iterable, Sequence

from ..arrays import mean, median, percentile, variance


@dataclass(frozen=True)
class RobustSummary:
    count: int
    mean: float
    median: float
    mad: float
    lower: float
    upper: float
    finite_fraction: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "count": self.count,
            "mean": self.mean,
            "median": self.median,
            "mad": self.mad,
            "lower": self.lower,
            "upper": self.upper,
            "finite_fraction": self.finite_fraction,
        }


def finite_values(values: Iterable[float]) -> list[float]:
    return [float(value) for value in values if math.isfinite(float(value))]


def median_absolute_deviation(values: Sequence[float]) -> float:
    center = median(values)
    return median(abs(value - center) for value in values)


def summarize(values: Iterable[float], lower_probability: float = 0.025, upper_probability: float = 0.975) -> RobustSummary:
    raw = list(values)
    finite = finite_values(raw)
    if not finite:
        return RobustSummary(len(raw), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return RobustSummary(
        count=len(raw),
        mean=mean(finite),
        median=median(finite),
        mad=median_absolute_deviation(finite),
        lower=percentile(finite, lower_probability),
        upper=percentile(finite, upper_probability),
        finite_fraction=len(finite) / max(1, len(raw)),
    )


def bootstrap_mean(values: Sequence[float], replicates: int, seed: int) -> list[float]:
    finite = finite_values(values)
    if not finite:
        return []
    rng = random.Random(seed)
    result = []
    for _ in range(max(1, replicates)):
        sample = [finite[rng.randrange(len(finite))] for _ in finite]
        result.append(mean(sample))
    return result


def winsorize(values: Sequence[float], lower_probability: float = 0.01, upper_probability: float = 0.99) -> list[float]:
    finite = finite_values(values)
    if not finite:
        return []
    lower = percentile(finite, lower_probability)
    upper = percentile(finite, upper_probability)
    return [min(upper, max(lower, value)) for value in values]


def weighted_mean(values: Sequence[float], weights: Sequence[float]) -> float:
    denominator = sum(max(0.0, weight) for weight in weights)
    if denominator == 0.0:
        return mean(values)
    return sum(value * max(0.0, weight) for value, weight in zip(values, weights)) / denominator


def weighted_variance(values: Sequence[float], weights: Sequence[float]) -> float:
    center = weighted_mean(values, weights)
    denominator = sum(max(0.0, weight) for weight in weights)
    if denominator == 0.0:
        return variance(values)
    return sum(max(0.0, weight) * (value - center) ** 2 for value, weight in zip(values, weights)) / denominator
