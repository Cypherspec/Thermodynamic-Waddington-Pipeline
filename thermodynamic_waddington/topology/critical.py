from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class CriticalPoint:
    node: int
    kind: str
    energy: float
    degree: int
    persistence_proxy: float


def classify_critical_points(energies: Sequence[float], outgoing: Sequence[Sequence[int]]) -> list[CriticalPoint]:
    points: list[CriticalPoint] = []
    for node, energy in enumerate(energies):
        neighbors = list(outgoing[node]) if node < len(outgoing) else []
        lower = sum(1 for other in neighbors if energies[other] < energy)
        higher = sum(1 for other in neighbors if energies[other] > energy)
        if lower == 0:
            kind = "minimum"
        elif higher == 0:
            kind = "maximum"
        else:
            kind = "saddle" if lower > 0 and higher > 0 else "regular"
        contrast = min([abs(energy - energies[other]) for other in neighbors] or [0.0])
        points.append(CriticalPoint(node, kind, float(energy), len(neighbors), float(contrast)))
    return points
