from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean
from .graph import Edge, adjacency


@dataclass(frozen=True)
class Trajectory:
    nodes: list[int]
    energies: list[float]
    work: float
    residence_time: float

    def to_dict(self) -> dict[str, object]:
        return {"nodes": self.nodes, "energies": self.energies, "work": self.work, "residence_time": self.residence_time}


def simulate_trajectories(edges: Sequence[Edge], energies: Sequence[float], starts: Sequence[int], steps: int = 24, temperature: float = 1.0, seed: int = 17) -> list[Trajectory]:
    rng = random.Random(seed)
    outgoing = adjacency(edges, len(energies))
    trajectories = []
    for start in starts:
        node = start
        nodes = [node]
        path_energies = [energies[node]]
        work = 0.0
        for _ in range(steps):
            candidates = outgoing[node]
            if not candidates:
                break
            weights = [math.exp(max(-50.0, min(50.0, -(energies[edge.target] - energies[node]) / max(temperature, 1e-9)))) for edge in candidates]
            threshold = rng.random() * sum(weights)
            selected = candidates[-1]
            for edge, weight in zip(candidates, weights):
                threshold -= weight
                if threshold <= 0:
                    selected = edge
                    break
            work += energies[selected.target] - energies[node]
            node = selected.target
            nodes.append(node)
            path_energies.append(energies[node])
        trajectories.append(Trajectory(nodes, path_energies, work, mean([1.0] * len(nodes))))
    return trajectories


def basin_occupancy(trajectories: Sequence[Trajectory], basins: Sequence[int]) -> dict[int, float]:
    counts = {basin: 0.0 for basin in basins}
    total = 0.0
    for trajectory in trajectories:
        for node in trajectory.nodes:
            if node in counts:
                counts[node] += 1.0
                total += 1.0
    if total:
        counts = {node: count / total for node, count in counts.items()}
    return counts
