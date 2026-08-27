from __future__ import annotations

import math
from typing import Sequence

from .arrays import mean


def _histogram(values: Sequence[float], bins: int) -> list[float]:
    if not values:
        return [0.0] * bins
    low, high = min(values), max(values)
    width = (high - low) / bins if high > low else 1.0
    result = [0.0] * bins
    for value in values:
        index = min(bins - 1, max(0, int((value - low) / width)))
        result[index] += 1.0
    total = sum(result)
    return [value / total for value in result] if total else result


def entropy(values: Sequence[float], bins: int = 16) -> float:
    return -sum(probability * math.log(probability) for probability in _histogram(values, bins) if probability > 0)


def mutual_information(first: Sequence[float], second: Sequence[float], bins: int = 16) -> float:
    if len(first) != len(second):
        raise ValueError("variables must have equal length")
    if not first:
        return 0.0
    one = _histogram(first, bins)
    two = _histogram(second, bins)
    joint = [[0.0] * bins for _ in range(bins)]
    low_a, high_a = min(first), max(first)
    low_b, high_b = min(second), max(second)
    width_a = (high_a - low_a) / bins if high_a > low_a else 1.0
    width_b = (high_b - low_b) / bins if high_b > low_b else 1.0
    for a, b in zip(first, second):
        ia = min(bins - 1, max(0, int((a - low_a) / width_a)))
        ib = min(bins - 1, max(0, int((b - low_b) / width_b)))
        joint[ia][ib] += 1.0 / len(first)
    return sum(joint[i][j] * math.log(joint[i][j] / max(one[i] * two[j], 1e-15)) for i in range(bins) for j in range(bins) if joint[i][j] > 0)


def transfer_entropy(source: Sequence[float], target: Sequence[float], bins: int = 8) -> float:
    if len(source) != len(target) or len(source) < 3:
        return 0.0
    target_now = target[1:]
    target_past = target[:-1]
    source_past = source[:-1]
    return max(0.0, mutual_information(source_past, target_now, bins) - mutual_information(target_past, target_now, bins))


def velocity_information(expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]]) -> list[float]:
    if len(expression) != len(velocity):
        raise ValueError("expression and velocity must have equal cell counts")
    result = []
    for row, vector in zip(expression, velocity):
        scale = math.sqrt(sum(value * value for value in row)) + 1e-12
        result.append(sum(abs(a * b) for a, b in zip(row, vector)) / scale)
    return result


def information_summary(points: Sequence[Sequence[float]], diffusions: Sequence[float], energies: Sequence[float]) -> dict[str, float]:
    scales = [math.log1p(max(0.0, value)) for value in diffusions]
    return {
        "energy_entropy": entropy(energies),
        "diffusion_entropy": entropy(diffusions),
        "energy_diffusion_mutual_information": mutual_information(energies, diffusions),
        "mean_coordinate_entropy": mean([entropy([row[index] for row in points]) for index in range(len(points[0]))]) if points else 0.0,
        "mean_signal_to_noise": mean([abs(energy) / max(1e-12, scale) for energy, scale in zip(energies, scales)]) if energies else 0.0,
    }
