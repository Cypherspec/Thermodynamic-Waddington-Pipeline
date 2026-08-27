from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..arrays import dot, norm, sub
from ..graph import Edge


@dataclass(frozen=True)
class CurrentDecomposition:
    edge_count: int
    total_current: float
    conservative_component: float
    rotational_component: float
    housekeeping_heat: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "edge_count": self.edge_count,
            "total_current": self.total_current,
            "conservative_component": self.conservative_component,
            "rotational_component": self.rotational_component,
            "housekeeping_heat": self.housekeeping_heat,
        }


def decompose_current(points: Sequence[Sequence[float]], velocities: Sequence[Sequence[float]], energies: Sequence[float], edges: Sequence[Edge], diffusion: Sequence[float]) -> CurrentDecomposition:
    total = 0.0
    conservative = 0.0
    rotational = 0.0
    housekeeping = 0.0
    for edge in edges:
        displacement = sub(points[edge.target], points[edge.source])
        distance = norm(displacement)
        if distance == 0.0:
            continue
        projected_velocity = dot(velocities[edge.source], displacement) / distance
        energy_slope = (energies[edge.target] - energies[edge.source]) / distance
        local_current = projected_velocity / max(diffusion[edge.source], 1e-12)
        conservative_part = -energy_slope / max(diffusion[edge.source], 1e-12)
        total += abs(local_current)
        conservative += abs(conservative_part)
        rotational += abs(local_current - conservative_part)
        housekeeping += max(0.0, local_current * (projected_velocity + energy_slope))
    return CurrentDecomposition(len(edges), total, conservative, rotational, housekeeping)
