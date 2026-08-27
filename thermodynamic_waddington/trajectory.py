from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import norm
from .graph import Edge, adjacency


@dataclass(frozen=True)
class Trajectory:
    start: int
    nodes: list[int]
    total_work: float
    confidence: float
    terminated: str

    def as_dict(self) -> dict:
        return {"start": self.start, "nodes": self.nodes, "total_work": self.total_work, "confidence": self.confidence, "terminated": self.terminated}


def _choose(outgoing: Sequence[Edge], energies: Sequence[float], visited: set[int]) -> Edge | None:
    candidates = [edge for edge in outgoing if edge.target not in visited]
    if not candidates:
        return None
    return max(candidates, key=lambda edge: edge.alignment - 0.15 * (energies[edge.target] - energies[edge.source]))


def integrate(edges: Sequence[Edge], energies: Sequence[float], starts: Sequence[int], max_steps: int = 64) -> list[Trajectory]:
    graph = adjacency(edges, len(energies))
    result = []
    for start in starts:
        node = start
        nodes = [node]
        visited = {node}
        work = 0.0
        alignments = []
        reason = "step_limit"
        for _ in range(max_steps):
            edge = _choose(graph[node], energies, visited)
            if edge is None:
                reason = "local_terminal"
                break
            work += energies[edge.target] - energies[edge.source]
            alignments.append(edge.alignment)
            node = edge.target
            nodes.append(node)
            visited.add(node)
        confidence = sum(alignments) / len(alignments) if alignments else 0.0
        result.append(Trajectory(start, nodes, work, confidence, reason))
    return result


def basin_membership(energies: Sequence[float], attractors: Sequence[int]) -> list[int]:
    if not energies or not attractors:
        return []
    return [min(attractors, key=lambda attractor: abs(energies[index] - energies[attractor])) for index in range(len(energies))]


def stochastic_paths(edges: Sequence[Edge], energies: Sequence[float], starts: Sequence[int], attractors: Sequence[int], temperature: float = 1.0, replicates: int = 8) -> list[Trajectory]:
    return integrate(edges, energies, starts, max_steps=max(1, len(energies)))


def trajectory_summary(trajectories: Sequence[Trajectory]) -> dict[str, float]:
    if not trajectories:
        return {"mean_length": 0.0, "mean_work": 0.0, "mean_confidence": 0.0}
    return {
        "mean_length": sum(len(item.nodes) for item in trajectories) / len(trajectories),
        "mean_work": sum(item.total_work for item in trajectories) / len(trajectories),
        "mean_confidence": sum(item.confidence for item in trajectories) / len(trajectories),
    }


def simulate_trajectories(edges: Sequence[Edge], energies: Sequence[float], starts: Sequence[int], attractors: Sequence[int] | None = None, temperature: float = 1.0, steps: int = 64, replicates: int = 8, seed: int = 7) -> list[Trajectory]:
    return integrate(edges, energies, starts, max_steps=steps)


def stochastic_paths(edges: Sequence[Edge], energies: Sequence[float], starts: Sequence[int], attractors: Sequence[int] | None = None, temperature: float = 1.0, replicates: int = 8) -> list[Trajectory]:
    return integrate(edges, energies, starts, max_steps=max(1, len(energies)))
