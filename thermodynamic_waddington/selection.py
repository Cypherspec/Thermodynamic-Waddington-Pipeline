from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, variance
from .graph import Edge


@dataclass(frozen=True)
class SelectionScore:
    node: int
    energy: float
    basin_mass: float
    accessibility: float
    uncertainty: float
    score: float
    reason: str

    def to_dict(self) -> dict[str, float | int | str]:
        return {"node": self.node, "energy": self.energy, "basin_mass": self.basin_mass, "accessibility": self.accessibility, "uncertainty": self.uncertainty, "score": self.score, "reason": self.reason}


def target_selection(energies: Sequence[float], edges: Sequence[Edge], uncertainty: Sequence[float] | None = None, top_k: int = 12) -> list[SelectionScore]:
    if not energies:
        return []
    uncertainty = uncertainty or [0.0] * len(energies)
    outgoing: dict[int, list[Edge]] = {}
    incoming: dict[int, int] = {}
    for edge in edges:
        outgoing.setdefault(edge.source, []).append(edge)
        incoming[edge.target] = incoming.get(edge.target, 0) + 1
    raw: list[SelectionScore] = []
    scale = max(1e-9, max(energies) - min(energies))
    for node, energy in enumerate(energies):
        downhill = sum(1 for edge in outgoing.get(node, []) if energies[edge.target] < energy)
        accessibility = 1.0 / (1.0 + downhill)
        mass = incoming.get(node, 0) / max(1, len(edges))
        confidence = 1.0 / (1.0 + uncertainty[node])
        basin = math_logistic(-(energy - min(energies)) / scale)
        score = 0.42 * accessibility + 0.28 * mass + 0.20 * basin + 0.10 * confidence
        reason = "low-energy, accessible, and recurrent" if score > 0.5 else "exploratory candidate requiring validation"
        raw.append(SelectionScore(node, energy, mass, accessibility, uncertainty[node], score, reason))
    return sorted(raw, key=lambda item: item.score, reverse=True)[:top_k]


def math_logistic(value: float) -> float:
    if value >= 0:
        exponent = pow(2.718281828, -value)
        return 1.0 / (1.0 + exponent)
    exponent = pow(2.718281828, value)
    return exponent / (1.0 + exponent)
