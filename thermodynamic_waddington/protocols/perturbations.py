from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..arrays import mean


@dataclass(frozen=True)
class PerturbationResult:
    name: str
    baseline_energy: float
    perturbed_energy: float
    delta_energy: float
    response: tuple[float, ...]

    def to_dict(self) -> dict[str, object]:
        return {"name": self.name, "baseline_energy": self.baseline_energy, "perturbed_energy": self.perturbed_energy, "delta_energy": self.delta_energy, "response": list(self.response)}


def apply_transcriptional_shift(energy: Sequence[float], weights: Sequence[float], strength: float, name: str) -> PerturbationResult:
    baseline = mean(energy)
    response = tuple(value - strength * weight for value, weight in zip(energy, weights))
    perturbed = mean(response)
    return PerturbationResult(name, baseline, perturbed, perturbed - baseline, response)


def rank_responses(results: Sequence[PerturbationResult]) -> list[PerturbationResult]:
    return sorted(results, key=lambda result: abs(result.delta_energy), reverse=True)
