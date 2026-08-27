from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean
from .graph import Edge


@dataclass(frozen=True)
class Intervention:
    source: int
    target: int
    score: float
    energy_drop: float
    alignment: float
    rationale: str

    def as_dict(self) -> dict:
        return {"source": self.source, "target": self.target, "score": self.score, "energy_drop": self.energy_drop, "alignment": self.alignment, "rationale": self.rationale}


def rank_interventions(edges: Sequence[Edge], energies: Sequence[float], labels: Sequence[str] | None = None, limit: int = 20) -> list[Intervention]:
    candidates = []
    for edge in edges:
        drop = energies[edge.source] - energies[edge.target]
        directional = max(0.0, edge.alignment)
        score = drop * (0.25 + directional) / (0.1 + edge.distance)
        rationale = "high downhill work with velocity-consistent direction" if score > 0 else "weak or uphill transition; retain as control"
        candidates.append(Intervention(edge.source, edge.target, score, drop, edge.alignment, rationale))
    return sorted(candidates, key=lambda item: item.score, reverse=True)[:max(0, limit)]


def transition_matrix(edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0) -> list[list[float]]:
    n = len(energies)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for edge in edges:
        delta = (energies[edge.target] - energies[edge.source]) / max(1e-9, temperature)
        matrix[edge.source][edge.target] = math.exp(max(-30.0, min(30.0, -delta))) * max(0.0, edge.alignment)
    for row in matrix:
        total = sum(row)
        if total:
            for index in range(n):
                row[index] /= total
    return matrix
