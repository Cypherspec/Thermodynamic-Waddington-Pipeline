from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..arrays import dot, mean, norm


@dataclass(frozen=True)
class FluxSummary:
    net_flux: tuple[float, ...]
    divergence: tuple[float, ...]
    circulation: float
    irreversible_fraction: float

    def to_dict(self) -> dict[str, object]:
        return {"net_flux": list(self.net_flux), "divergence": list(self.divergence), "circulation": self.circulation, "irreversible_fraction": self.irreversible_fraction}


def estimate_flux(points: Sequence[Sequence[float]], velocities: Sequence[Sequence[float]], edges: Sequence[tuple[int, int]]) -> FluxSummary:
    n = len(points)
    net = [0.0] * n
    divergence = [0.0] * n
    directed = 0.0
    absolute = 0.0
    for source, target in edges:
        displacement = [target_value - source_value for source_value, target_value in zip(points[source], points[target])]
        length = norm(displacement) or 1e-12
        projection = dot(velocities[source], displacement) / length
        directed += projection
        absolute += abs(projection)
        net[source] -= projection
        net[target] += projection
    for index in range(n):
        divergence[index] = net[index]
    circulation = abs(directed) / max(1, len(edges))
    irreversible = 0.5 * (1.0 + directed / max(absolute, 1e-12))
    return FluxSummary(tuple(net), tuple(divergence), circulation, max(0.0, min(1.0, irreversible)))
