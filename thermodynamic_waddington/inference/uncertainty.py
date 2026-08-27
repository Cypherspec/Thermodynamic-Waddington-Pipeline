from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from ..arrays import mean, percentile


@dataclass(frozen=True)
class ConfidenceBand:
    lower: tuple[float, ...]
    median: tuple[float, ...]
    upper: tuple[float, ...]
    coverage: float

    def to_dict(self) -> dict[str, object]:
        return {"lower": list(self.lower), "median": list(self.median), "upper": list(self.upper), "coverage": self.coverage}


def percentile_band(samples: Sequence[Sequence[float]], lower: float = 0.05, upper: float = 0.95) -> ConfidenceBand:
    if not samples:
        return ConfidenceBand((), (), (), upper - lower)
    width = len(samples[0])
    columns = [[sample[index] for sample in samples] for index in range(width)]
    return ConfidenceBand(tuple(percentile(column, lower) for column in columns), tuple(percentile(column, 0.5) for column in columns), tuple(percentile(column, upper) for column in columns), upper - lower)


def calibration_error(observed: Sequence[float], predicted: Sequence[float], bins: int = 10) -> float:
    if len(observed) != len(predicted) or not observed:
        return float("nan")
    errors = []
    for bucket in range(bins):
        left = bucket / bins
        right = (bucket + 1) / bins
        indices = [index for index, value in enumerate(predicted) if left <= value < right or bucket == bins - 1 and value == right]
        if indices:
            errors.append(abs(mean([observed[index] for index in indices]) - mean([predicted[index] for index in indices])))
    return mean(errors) if errors else 0.0


def confidence_width(band: ConfidenceBand) -> float:
    if not band.lower:
        return 0.0
    return mean([upper - lower for lower, upper in zip(band.lower, band.upper)])
