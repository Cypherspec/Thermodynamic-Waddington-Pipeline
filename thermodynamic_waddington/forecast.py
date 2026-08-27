from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean
from .graph import Edge, adjacency


@dataclass(frozen=True)
class FateForecast:
    origin: int
    probabilities: list[float]
    expected_steps: float
    entropy: float
    terminal_nodes: list[int]


def _normalize(values: Sequence[float]) -> list[float]:
    total = sum(max(0.0, value) for value in values)
    return [max(0.0, value) / total for value in values] if total else [1.0 / len(values) for _ in values]


def transition_kernel(edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0, irreversibility: float = 0.0) -> list[list[tuple[int, float]]]:
    outgoing = adjacency(edges, len(energies))
    result: list[list[tuple[int, float]]] = []
    for source, candidates in enumerate(outgoing):
        logits = []
        targets = []
        for edge in candidates:
            downhill = (energies[source] - energies[edge.target]) / max(temperature, 1e-12)
            bias = irreversibility * edge.alignment
            logits.append(math.exp(max(-50.0, min(50.0, downhill + bias))))
            targets.append(edge.target)
        weights = _normalize(logits)
        result.append(list(zip(targets, weights)))
    return result


def propagate(origin: int, kernel: Sequence[Sequence[tuple[int, float]]], energies: Sequence[float], steps: int = 80) -> list[float]:
    state = [0.0] * len(energies)
    state[origin] = 1.0
    for _ in range(max(1, steps)):
        next_state = [0.0] * len(state)
        for source, mass in enumerate(state):
            if mass == 0.0 or not kernel[source]:
                next_state[source] += mass
                continue
            for target, probability in kernel[source]:
                next_state[target] += mass * probability
        state = next_state
    return state


def forecast(origins: Sequence[int], edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0, steps: int = 80) -> list[FateForecast]:
    kernel = transition_kernel(edges, energies, temperature)
    terminal = sorted(range(len(energies)), key=lambda index: energies[index])[:max(1, min(8, len(energies)))]
    reports: list[FateForecast] = []
    for origin in origins:
        probabilities = propagate(origin, kernel, energies, steps)
        terminal_probabilities = [probabilities[index] for index in terminal]
        total = sum(terminal_probabilities)
        normalized = [value / total for value in terminal_probabilities] if total else [0.0 for _ in terminal]
        entropy = -sum(value * math.log(value) for value in normalized if value > 0)
        reports.append(FateForecast(origin, normalized, float(steps), entropy, terminal))
    return reports


def sample_trajectory(start: int, edges: Sequence[Edge], energies: Sequence[float], steps: int = 30, seed: int = 7) -> list[int]:
    rng = random.Random(seed)
    kernel = transition_kernel(edges, energies)
    path = [start]
    current = start
    for _ in range(steps):
        if not kernel[current]:
            break
        threshold = rng.random()
        cumulative = 0.0
        next_node = current
        for target, probability in kernel[current]:
            cumulative += probability
            if threshold <= cumulative:
                next_node = target
                break
        path.append(next_node)
        current = next_node
    return path


def forecast_fates(energies: Sequence[float], labels: Sequence[str], temperature: float = 1.0) -> dict[str, object]:
    grouped: dict[str, list[float]] = {}
    for energy, label in zip(energies, labels):
        grouped.setdefault(label, []).append(energy)
    scores = {label: math.exp(-mean(values) / max(temperature, 1e-12)) for label, values in grouped.items() if values}
    total = sum(scores.values()) or 1.0
    probabilities = {label: value / total for label, value in scores.items()}
    return {"probabilities": probabilities, "entropy": -sum(value * math.log(value) for value in probabilities.values() if value > 0), "most_likely": max(probabilities, key=probabilities.get) if probabilities else None}
