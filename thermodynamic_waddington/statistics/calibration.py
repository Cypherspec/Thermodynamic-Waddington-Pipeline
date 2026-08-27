from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from ..arrays import mean


@dataclass(frozen=True)
class CalibrationBin:
    lower: float
    upper: float
    count: int
    predicted: float
    observed: float
    absolute_error: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "lower": self.lower,
            "upper": self.upper,
            "count": self.count,
            "predicted": self.predicted,
            "observed": self.observed,
            "absolute_error": self.absolute_error,
        }


def brier_score(predicted: Sequence[float], observed: Sequence[int]) -> float:
    if len(predicted) != len(observed):
        raise ValueError("predicted and observed must have equal lengths")
    if not predicted:
        return 0.0
    return mean((p - o) ** 2 for p, o in zip(predicted, observed))


def log_loss(predicted: Sequence[float], observed: Sequence[int], floor: float = 1e-9) -> float:
    if len(predicted) != len(observed):
        raise ValueError("predicted and observed must have equal lengths")
    if not predicted:
        return 0.0
    return mean(-o * math.log(max(floor, p)) - (1 - o) * math.log(max(floor, 1 - p)) for p, o in zip(predicted, observed))


def reliability_bins(predicted: Sequence[float], observed: Sequence[int], bins: int = 10) -> list[CalibrationBin]:
    if len(predicted) != len(observed):
        raise ValueError("predicted and observed must have equal lengths")
    result = []
    for index in range(max(1, bins)):
        lower = index / max(1, bins)
        upper = (index + 1) / max(1, bins)
        selected = [i for i, probability in enumerate(predicted) if lower <= probability <= upper or (index == bins - 1 and probability <= upper)]
        if not selected:
            continue
        prediction = mean(predicted[i] for i in selected)
        observation = mean(float(observed[i]) for i in selected)
        result.append(CalibrationBin(lower, upper, len(selected), prediction, observation, abs(prediction - observation)))
    return result
