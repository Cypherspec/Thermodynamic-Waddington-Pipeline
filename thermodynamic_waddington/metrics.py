from __future__ import annotations

import math
from typing import Sequence

from .arrays import mean, variance
from .graph import Edge


def edge_work(edges: Sequence[Edge]) -> list[float]:
    return [edge.work for edge in edges]


def mean_alignment(edges: Sequence[Edge]) -> float:
    return mean([edge.alignment for edge in edges]) if edges else 0.0


def irreversibility_index(edges: Sequence[Edge]) -> float:
    if not edges:
        return 0.0
    forward = sum(max(0.0, edge.work) for edge in edges)
    backward = sum(max(0.0, -edge.work) for edge in edges)
    return (forward - backward) / (forward + backward + 1e-12)


def landscape_roughness(energies: Sequence[float]) -> float:
    if len(energies) < 2:
        return 0.0
    differences = [energies[i + 1] - energies[i] for i in range(len(energies) - 1)]
    return math.sqrt(variance(differences))


def basin_contrast(energies: Sequence[float], attractors: Sequence[int]) -> float:
    if not energies or not attractors:
        return 0.0
    baseline = mean(energies)
    basin = mean([energies[i] for i in attractors if 0 <= i < len(energies)])
    return baseline - basin
