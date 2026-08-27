from __future__ import annotations

import math
from typing import Sequence

from .arrays import mean, normalize, safe_log


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        return 0.0
    mx, my = mean(x), mean(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denominator_x = math.sqrt(sum((a - mx) ** 2 for a in x))
    denominator_y = math.sqrt(sum((b - my) ** 2 for b in y))
    denominator = denominator_x * denominator_y
    return numerator / denominator if denominator else 0.0


def fate_probabilities(energies: Sequence[float], fate_labels: Sequence[str], temperature: float) -> dict[str, float]:
    groups: dict[str, list[float]] = {}
    for energy, label in zip(energies, fate_labels):
        groups.setdefault(label, []).append(energy)
    scores = {label: math.exp(-mean(values) / temperature) for label, values in groups.items()}
    total = sum(scores.values()) or 1.0
    return {label: value / total for label, value in scores.items()}


def empirical_probabilities(labels: Sequence[str]) -> dict[str, float]:
    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    total = len(labels) or 1
    return {label: count / total for label, count in counts.items()}


def fate_calibration(energies: Sequence[float], labels: Sequence[str], outcomes: Sequence[str], temperature: float) -> dict[str, object]:
    predicted = fate_probabilities(energies, labels, temperature)
    observed = empirical_probabilities(outcomes)
    keys = sorted(set(predicted) | set(observed))
    return {"fates": keys, "predicted": [predicted.get(key, 0.0) for key in keys], "observed": [observed.get(key, 0.0) for key in keys], "pearson": pearson([predicted.get(key, 0.0) for key in keys], [observed.get(key, 0.0) for key in keys])}


def attractor_depth_order(energies: Sequence[float], labels: Sequence[str]) -> list[dict[str, float | str]]:
    groups: dict[str, list[float]] = {}
    for energy, label in zip(energies, labels):
        groups.setdefault(label, []).append(energy)
    result = [{"label": label, "mean_free_energy": mean(values), "n": len(values)} for label, values in groups.items()]
    return sorted(result, key=lambda item: float(item["mean_free_energy"]))


def velocity_alignment_score(edges: Sequence[dict[str, float | int]]) -> float:
    if not edges:
        return 0.0
    return mean(float(edge["alignment"]) for edge in edges)


def reversibility_gap(energies: Sequence[float], edges: Sequence[dict[str, float | int]]) -> float:
    pairs = {(int(edge["source"]), int(edge["target"])): edge for edge in edges}
    gaps = []
    for (source, target), edge in pairs.items():
        reverse = pairs.get((target, source))
        if reverse:
            observed = float(edge["work"]) + float(reverse["work"])
            expected = (energies[target] - energies[source]) + (energies[source] - energies[target])
            gaps.append(abs(observed - expected))
    return mean(gaps)
