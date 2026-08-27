from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .arrays import mean
from .graph import Edge, adjacency


@dataclass(frozen=True)
class BifurcationPoint:
    parameter: float
    basin_count: int
    minimum_energy: float
    curvature_proxy: float


def continuation_scan(edges: Sequence[Edge], base_energies: Sequence[float], temperatures: Sequence[float]) -> list[BifurcationPoint]:
    outgoing = adjacency(edges, len(base_energies))
    results = []
    for temperature in temperatures:
        minima = []
        curvatures = []
        for node, neighbors in enumerate(outgoing):
            if neighbors and all(base_energies[node] <= base_energies[edge.target] for edge in neighbors):
                minima.append(node)
                curvatures.append(mean([base_energies[edge.target] - base_energies[node] for edge in neighbors]))
        results.append(BifurcationPoint(temperature, len(minima), min((base_energies[node] for node in minima), default=0.0), mean(curvatures) if curvatures else 0.0))
    return results


def detect_transitions(scan: Sequence[BifurcationPoint]) -> list[tuple[float, float]]:
    transitions = []
    for previous, current in zip(scan, scan[1:]):
        if previous.basin_count != current.basin_count:
            transitions.append((previous.parameter, current.parameter))
    return transitions
