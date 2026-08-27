from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .landscape import barrier_heights


@dataclass(frozen=True)
class ContinuationStep:
    parameter: float
    attractor_count: int
    barrier_count: int
    mean_barrier: float


def continue_landscape(energies: Sequence[float], edges, parameter_grid: Sequence[float]) -> list[ContinuationStep]:
    base = list(energies)
    steps: list[ContinuationStep] = []
    for parameter in parameter_grid:
        transformed = [value * (1.0 + float(parameter)) for value in base]
        barriers = barrier_heights(edges, transformed)
        attractors = sum(1 for edge in edges if transformed[edge.source] <= transformed[edge.target])
        steps.append(ContinuationStep(float(parameter), attractors, len(barriers), sum(barriers) / len(barriers) if barriers else 0.0))
    return steps
