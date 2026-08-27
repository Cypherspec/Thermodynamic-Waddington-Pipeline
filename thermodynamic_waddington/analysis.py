from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, quantile, variance
from .graph import Edge, adjacency
from .trajectory import basin_membership, trajectory_summary, stochastic_paths
from .uncertainty import bootstrap


@dataclass(frozen=True)
class LandscapeReport:
    score: float
    basin_sizes: dict[int, int]
    transition_flux: list[dict]
    critical_edges: list[dict]
    uncertainty: dict

    def as_dict(self) -> dict:
        return {"score": self.score, "basin_sizes": self.basin_sizes, "transition_flux": self.transition_flux, "critical_edges": self.critical_edges, "uncertainty": self.uncertainty}


def entropy(probabilities: Sequence[float]) -> float:
    return -sum(p * math.log(p) for p in probabilities if p > 0.0)


def basin_statistics(labels: Sequence[int], energies: Sequence[float], velocities: Sequence[Sequence[float]]) -> dict[int, dict]:
    groups: dict[int, list[int]] = {}
    for index, label in enumerate(labels):
        groups.setdefault(label, []).append(index)
    return {label: {"size": len(nodes), "energy_mean": mean([energies[i] for i in nodes]), "energy_spread": math.sqrt(max(0.0, variance([energies[i] for i in nodes]))), "velocity_mean": mean([math.sqrt(sum(value * value for value in velocities[i])) for i in nodes])} for label, nodes in groups.items()}


def flux_matrix(labels: Sequence[int], edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0) -> list[list[float]]:
    basins = sorted({label for label in labels if label >= 0})
    index = {label: pos for pos, label in enumerate(basins)}
    matrix = [[0.0 for _ in basins] for _ in basins]
    for edge in edges:
        left, right = labels[edge.source], labels[edge.target]
        if left >= 0 and right >= 0 and left != right:
            matrix[index[left]][index[right]] += math.exp(-max(0.0, edge.work) / max(temperature, 1e-9))
    return matrix


def critical_edges(edges: Sequence[Edge], energies: Sequence[float], limit: int = 20) -> list[dict]:
    ranked = sorted(edges, key=lambda edge: (edge.work + max(0.0, energies[edge.target] - energies[edge.source])), reverse=True)
    return [{"source": edge.source, "target": edge.target, "work": edge.work, "barrier": max(0.0, energies[edge.target] - energies[edge.source]), "alignment": edge.alignment} for edge in ranked[:limit]]


def score_landscape(edges: Sequence[Edge], energies: Sequence[float], labels: Sequence[int]) -> float:
    if not edges:
        return 0.0
    downhill = sum(1 for edge in edges if energies[edge.target] <= energies[edge.source]) / len(edges)
    aligned = mean([edge.alignment for edge in edges])
    basin_count = len({label for label in labels if label >= 0})
    complexity = min(1.0, basin_count / max(1.0, len(labels) ** 0.5))
    return 0.45 * downhill + 0.4 * aligned + 0.15 * complexity


def report(edges: Sequence[Edge], energies: Sequence[float], velocities: Sequence[Sequence[float]], attractors: Sequence[int], temperature: float = 1.0, draws: int = 128) -> LandscapeReport:
    labels = basin_membership(edges, energies, attractors)
    groups = {label: labels.count(label) for label in sorted(set(labels)) if label >= 0}
    paths = [stochastic_paths(edges, energies, attractor, draws=draws, temperature=temperature, seed=101 + attractor) for attractor in attractors]
    flux = [{"source": attractor, "summary": trajectory_summary(samples)} for attractor, samples in zip(attractors, paths)]
    work = [edge.work for edge in edges]
    return LandscapeReport(score_landscape(edges, energies, labels), groups, flux, critical_edges(edges, energies), {"edge_work": bootstrap(work, draws=min(256, max(32, draws)), seed=79).as_dict()})
