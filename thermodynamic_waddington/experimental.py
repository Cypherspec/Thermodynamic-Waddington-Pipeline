from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, percentile
from .graph import Edge


@dataclass(frozen=True)
class TransitionRisk:
    source: int
    target: int
    activation: float
    probability: float
    entropy: float
    metastability: float

    def to_dict(self) -> dict[str, float | int]:
        return self.__dict__.copy()


def transition_risk(edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0) -> list[TransitionRisk]:
    if not edges:
        return []
    scale = max(1e-9, temperature)
    result: list[TransitionRisk] = []
    outgoing: dict[int, list[Edge]] = {}
    for edge in edges:
        outgoing.setdefault(edge.source, []).append(edge)
    for source, source_edges in outgoing.items():
        logits = [-(energies[edge.target] - energies[source]) / scale for edge in source_edges]
        maximum = max(logits)
        weights = [math.exp(logit - maximum) for logit in logits]
        normalizer = sum(weights)
        for edge, weight in zip(source_edges, weights):
            probability = weight / max(1e-12, normalizer)
            activation = max(0.0, energies[edge.target] - energies[source])
            entropy = -probability * math.log(max(1e-12, probability))
            metastability = 1.0 / (1.0 + activation + edge.local_diffusion)
            result.append(TransitionRisk(source, edge.target, activation, probability, entropy, metastability))
    return result


def rare_transition_candidates(risks: Sequence[TransitionRisk], quantile: float = 0.9) -> list[TransitionRisk]:
    if not risks:
        return []
    cutoff = percentile([risk.activation for risk in risks], quantile)
    return [risk for risk in risks if risk.activation >= cutoff]


def path_entropy(risks: Sequence[TransitionRisk]) -> float:
    return sum(risk.entropy for risk in risks)


def basin_free_energy(energies: Sequence[float], members: Sequence[int], temperature: float = 1.0) -> float:
    if not members:
        return 0.0
    scale = max(1e-9, temperature)
    weights = [math.exp(-energies[index] / scale) for index in members]
    return -scale * math.log(max(1e-12, sum(weights)))


def estimate_escape_time(risks: Sequence[TransitionRisk], source: int) -> float:
    candidates = [risk for risk in risks if risk.source == source]
    if not candidates:
        return float("inf")
    rate = sum(risk.probability * math.exp(-risk.activation) for risk in candidates)
    return 1.0 / max(1e-12, rate)


def summarize_risk(risks: Sequence[TransitionRisk]) -> dict[str, float]:
    if not risks:
        return {"mean_activation": 0.0, "mean_probability": 0.0, "mean_metastability": 0.0, "path_entropy": 0.0}
    return {"mean_activation": mean([risk.activation for risk in risks]), "mean_probability": mean([risk.probability for risk in risks]), "mean_metastability": mean([risk.metastability for risk in risks]), "path_entropy": path_entropy(risks)}
