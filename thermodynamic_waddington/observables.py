from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, percentile, safe_log
from .graph import Edge


@dataclass(frozen=True)
class ObservableReport:
    entropy: float
    effective_temperature: float
    basin_occupancy: dict[str, float]
    current_entropy_production: float
    fluctuation_asymmetry: float

    def to_dict(self) -> dict[str, object]:
        return {"entropy": self.entropy, "effective_temperature": self.effective_temperature, "basin_occupancy": self.basin_occupancy, "current_entropy_production": self.current_entropy_production, "fluctuation_asymmetry": self.fluctuation_asymmetry}


def _entropy(probabilities: Sequence[float]) -> float:
    return -sum(value * safe_log(value) for value in probabilities if value > 0)


def report(energies: Sequence[float], edges: Sequence[Edge], temperature: float, labels: Sequence[str] | None = None) -> ObservableReport:
    weights = [math.exp(max(-60.0, min(60.0, -energy / max(temperature, 1e-12)))) for energy in energies]
    total = sum(weights) or 1.0
    probabilities = [weight / total for weight in weights]
    occupancy: dict[str, float] = {}
    for index, probability in enumerate(probabilities):
        key = labels[index] if labels and index < len(labels) else "unlabeled"
        occupancy[key] = occupancy.get(key, 0.0) + probability
    entropy_production = sum(abs(edge.alignment) / max(edge.distance, 1e-9) for edge in edges) / max(1, len(edges))
    positive = [edge.alignment for edge in edges if edge.alignment > 0]
    negative = [-edge.alignment for edge in edges if edge.alignment < 0]
    asymmetry = abs(mean(positive) - mean(negative)) if positive and negative else 0.0
    effective_temperature = percentile(energies, 0.75) - percentile(energies, 0.25)
    return ObservableReport(_entropy(probabilities), effective_temperature, occupancy, entropy_production, asymmetry)
