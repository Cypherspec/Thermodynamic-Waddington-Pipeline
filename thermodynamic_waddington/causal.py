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


def counterfactual_ranking(edges: Sequence[Edge], energies: Sequence[float], damping: float = 0.5, limit: int = 10) -> list[int]:
    # O(E) equivalent of rank_targets(intervention_scan(...)). The scan damps
    # edges incident to each candidate node and scores nodes by mean
    # |delta_work| * confidence over the edges that land on them. Every scan
    # target contributes the same term to an edge unless it is incident, so the
    # full T x E sweep collapses to a single scatter over edges.
    import numpy as np
    if not edges:
        return []
    src = np.fromiter((e.source for e in edges), np.int64, len(edges))
    tgt = np.fromiter((e.target for e in edges), np.int64, len(edges))
    algn = np.fromiter((e.alignment for e in edges), float, len(edges))
    en = np.asarray(energies, float)
    base = en[tgt] - en[src]
    ab = np.abs(base)
    f0 = np.maximum(0.01, algn + 1.0)
    fd = np.maximum(0.01, algn * damping + 1.0)
    d0 = base * (f0 - 1.0)
    dd = base * (fd - 1.0)
    term0 = np.abs(d0) * np.minimum(1.0, np.abs(d0) / (1.0 + ab))
    termd = np.abs(dd) * np.minimum(1.0, np.abs(dd) / (1.0 + ab))
    tset = set(int(x) for x in tgt.tolist())
    total_targets = len(tset)
    src_in = np.fromiter((1.0 if int(x) in tset else 0.0 for x in src.tolist()), float, len(edges))
    dcount = np.where(src_in > 0.0, 2.0, 1.0)
    contrib = dcount * termd + (total_targets - dcount) * term0
    size = int(tgt.max()) + 1
    gtot = np.zeros(size)
    gcnt = np.zeros(size)
    np.add.at(gtot, tgt, contrib)
    np.add.at(gcnt, tgt, float(total_targets))
    order: list[int] = []
    seen: set[int] = set()
    for x in tgt.tolist():
        xi = int(x)
        if xi not in seen:
            seen.add(xi)
            order.append(xi)
    score = {m: gtot[m] / gcnt[m] for m in order}
    return sorted(order, key=lambda m: score[m], reverse=True)[:limit]


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
