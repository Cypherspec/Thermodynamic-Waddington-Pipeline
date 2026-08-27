from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean
from .graph import Edge


@dataclass(frozen=True)
class InterventionEffect:
    source: int
    target: int
    baseline_work: float
    perturbed_work: float
    delta_work: float
    confidence: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "source": self.source,
            "target": self.target,
            "baseline_work": self.baseline_work,
            "perturbed_work": self.perturbed_work,
            "delta_work": self.delta_work,
            "confidence": self.confidence,
        }


def edge_work(edge: Edge, energies: Sequence[float]) -> float:
    return energies[edge.target] - energies[edge.source]


def counterfactual_edges(edges: Sequence[Edge], target_node: int, damping: float = 0.5) -> list[Edge]:
    result: list[Edge] = []
    for edge in edges:
        if edge.target == target_node or edge.source == target_node:
            result.append(Edge(edge.source, edge.target, edge.distance, edge.alignment * damping))
        else:
            result.append(edge)
    return result


def intervention_scan(edges: Sequence[Edge], energies: Sequence[float], targets: Sequence[int] | None = None, damping: float = 0.5) -> list[InterventionEffect]:
    if targets is None:
        targets = sorted({edge.target for edge in edges})
    baseline = {(edge.source, edge.target): edge_work(edge, energies) for edge in edges}
    effects: list[InterventionEffect] = []
    for target in targets:
        perturbed = counterfactual_edges(edges, target, damping)
        for edge in perturbed:
            if (edge.source, edge.target) not in baseline:
                continue
            new_work = edge_work(edge, energies) * max(0.01, edge.alignment + 1.0)
            old_work = baseline[(edge.source, edge.target)]
            delta = new_work - old_work
            confidence = min(1.0, abs(delta) / (1.0 + abs(old_work)))
            effects.append(InterventionEffect(edge.source, edge.target, old_work, new_work, delta, confidence))
    return effects


def rank_targets(effects: Sequence[InterventionEffect], limit: int = 10) -> list[int]:
    scores: dict[int, list[float]] = {}
    for effect in effects:
        scores.setdefault(effect.target, []).append(abs(effect.delta_work) * effect.confidence)
    ranked = sorted(scores, key=lambda target: mean(scores[target]), reverse=True)
    return ranked[:limit]


def intervention_matrix(edges: Sequence[Edge], energies: Sequence[float], targets: Sequence[int], damping: float = 0.5) -> list[list[float]]:
    effects = intervention_scan(edges, energies, targets, damping)
    by_target: dict[int, list[InterventionEffect]] = {target: [] for target in targets}
    for effect in effects:
        for target in targets:
            if effect.target == target:
                by_target[target].append(effect)
    return [[mean([abs(effect.delta_work) for effect in by_target[target]]) if by_target[target] else 0.0 for target in targets]]


def bootstrap_target_scores(edges: Sequence[Edge], energies: Sequence[float], targets: Sequence[int], replicates: int = 32, seed: int = 7) -> dict[int, dict[str, float]]:
    rng = random.Random(seed)
    samples: dict[int, list[float]] = {target: [] for target in targets}
    for _ in range(max(1, replicates)):
        resampled = [edges[rng.randrange(len(edges))] for _ in range(len(edges))] if edges else []
        effects = intervention_scan(resampled, energies, targets)
        for target in targets:
            values = [abs(effect.delta_work) for effect in effects if effect.target == target]
            samples[target].append(mean(values) if values else 0.0)
    summary: dict[int, dict[str, float]] = {}
    for target, values in samples.items():
        ordered = sorted(values)
        low = ordered[int(0.025 * (len(ordered) - 1))]
        high = ordered[int(0.975 * (len(ordered) - 1))]
        summary[target] = {"mean": mean(values), "lower": low, "upper": high, "stability": 1.0 - (high - low) / (1.0 + abs(mean(values)))}
    return summary


def softmax_policy(scores: Sequence[float], temperature: float = 1.0) -> list[float]:
    if not scores:
        return []
    scale = max(1e-9, temperature)
    shifted = [score / scale for score in scores]
    maximum = max(shifted)
    weights = [math.exp(value - maximum) for value in shifted]
    total = sum(weights)
    return [weight / total for weight in weights]
