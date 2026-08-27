from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, percentile, variance


@dataclass(frozen=True)
class Interval:
    median: float
    lower: float
    upper: float
    samples: int

    def to_dict(self) -> dict[str, float | int]:
        return {"median": self.median, "lower": self.lower, "upper": self.upper, "samples": self.samples}


def interval(values: Sequence[float], lower: float = 0.025, upper: float = 0.975) -> Interval:
    if not values:
        return Interval(0.0, 0.0, 0.0, 0)
    return Interval(percentile(values, 0.5), percentile(values, lower), percentile(values, upper), len(values))


def effective_sample_size(weights: Sequence[float]) -> float:
    total = sum(weights)
    if total <= 0:
        return 0.0
    normalized = [weight / total for weight in weights]
    return 1.0 / sum(weight * weight for weight in normalized)


def coefficient_of_variation(values: Sequence[float]) -> float:
    average = mean(values)
    return math.sqrt(variance(values)) / max(1e-12, abs(average))


def calibrate_scores(scores: Sequence[float], labels: Sequence[int], bins: int = 10) -> dict[str, object]:
    if len(scores) != len(labels):
        raise ValueError("scores and labels must have equal length")
    groups: list[dict[str, float | int]] = []
    for bucket in range(bins):
        lo = bucket / bins
        hi = (bucket + 1) / bins
        selected = [label for score, label in zip(scores, labels) if lo <= score < hi or bucket == bins - 1 and score <= hi]
        if selected:
            groups.append({"lower": lo, "upper": hi, "count": len(selected), "predicted": (lo + hi) / 2, "observed": mean(selected)})
    error = mean([abs(float(group["predicted"]) - float(group["observed"])) for group in groups]) if groups else 0.0
    return {"bins": groups, "expected_calibration_error": error}


def bootstrap(values: Sequence[float], draws: int = 128, seed: int = 7) -> Interval:
    import random
    if not values:
        return interval([])
    rng = random.Random(seed)
    samples = [mean(rng.choice(list(values)) for _ in values) for _ in range(max(1, draws))]
    return interval(samples)


def bootstrap(values: Sequence[float], draws: int = 128, seed: int = 7):
    from dataclasses import dataclass
    import random
    @dataclass(frozen=True)
    class BootstrapSummary:
        mean: float
        lower: float
        upper: float
        standard_error: float
        draws: int
        def as_dict(self) -> dict[str, float | int]:
            return {"mean": self.mean, "lower": self.lower, "upper": self.upper, "standard_error": self.standard_error, "draws": self.draws}
    if not values:
        return BootstrapSummary(0.0, 0.0, 0.0, 0.0, 0)
    rng = random.Random(seed)
    samples = [sum(values[rng.randrange(len(values))] for _ in values) / len(values) for _ in range(max(1, draws))]
    center = sum(samples) / len(samples)
    spread = (sum((value - center) ** 2 for value in samples) / len(samples)) ** 0.5
    return BootstrapSummary(center, sorted(samples)[max(0, int(len(samples) * 0.025) - 1)], sorted(samples)[min(len(samples) - 1, int(len(samples) * 0.975))], spread, len(samples))
