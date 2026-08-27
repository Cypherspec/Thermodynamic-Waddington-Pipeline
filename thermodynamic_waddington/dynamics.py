from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, norm
from .graph import Edge


@dataclass(frozen=True)
class Trajectory:
    nodes: tuple[int, ...]
    energies: tuple[float, ...]
    probability: float
    work: float

    def to_dict(self) -> dict[str, object]:
        return {"nodes": list(self.nodes), "energies": list(self.energies), "probability": self.probability, "work": self.work}


def transition_probabilities(edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0) -> dict[int, list[tuple[int, float]]]:
    scale = max(1e-9, temperature)
    buckets: dict[int, list[Edge]] = {}
    for edge in edges:
        buckets.setdefault(edge.source, []).append(edge)
    output: dict[int, list[tuple[int, float]]] = {}
    for source, candidates in buckets.items():
        logits = [-(energies[edge.target] - energies[source]) / scale for edge in candidates]
        anchor = max(logits) if logits else 0.0
        scores = [math.exp(max(-60.0, min(60.0, value - anchor))) for value in logits]
        normalizer = sum(scores) or 1.0
        output[source] = [(edge.target, score / normalizer) for edge, score in zip(candidates, scores)]
    return output


def transition_risk(edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0) -> list[dict[str, float | int]]:
    probabilities = transition_probabilities(edges, energies, temperature)
    lookup = {(source, target): probability for source, items in probabilities.items() for target, probability in items}
    return [{"source": edge.source, "target": edge.target, "activation": max(0.0, energies[edge.target] - energies[edge.source]), "probability": lookup.get((edge.source, edge.target), 0.0)} for edge in edges]


def summarize_risk(risks: Sequence[dict[str, float | int]]) -> dict[str, float]:
    if not risks:
        return {"mean_activation": 0.0, "mean_probability": 0.0}
    return {"mean_activation": mean(float(item["activation"]) for item in risks), "mean_probability": mean(float(item["probability"]) for item in risks)}


def simulate_trajectories(edges: Sequence[Edge], energies: Sequence[float], starts: Sequence[int], steps: int = 25, replicates: int = 32, temperature: float = 1.0, seed: int = 3) -> tuple[Trajectory, ...]:
    rng = random.Random(seed)
    transitions = transition_probabilities(edges, energies, temperature)
    trajectories: list[Trajectory] = []
    for start in starts:
        for _ in range(replicates):
            node = start
            nodes = [node]
            probability = 1.0
            for _ in range(steps):
                choices = transitions.get(node)
                if not choices:
                    break
                threshold = rng.random()
                cumulative = 0.0
                selected = choices[-1]
                for choice in choices:
                    cumulative += choice[1]
                    if threshold <= cumulative:
                        selected = choice
                        break
                node, chance = selected
                probability *= chance
                nodes.append(node)
            path_energy = [energies[index] for index in nodes]
            work = sum(max(0.0, path_energy[i + 1] - path_energy[i]) for i in range(len(path_energy) - 1))
            trajectories.append(Trajectory(tuple(nodes), tuple(path_energy), probability, work))
    return tuple(trajectories)


def trajectory_summary(trajectories: Sequence[Trajectory]) -> dict[str, float]:
    if not trajectories:
        return {"mean_length": 0.0, "mean_work": 0.0, "mean_probability": 0.0, "terminal_energy": 0.0}
    return {"mean_length": mean([len(path.nodes) for path in trajectories]), "mean_work": mean([path.work for path in trajectories]), "mean_probability": mean([path.probability for path in trajectories]), "terminal_energy": mean([path.energies[-1] for path in trajectories])}
