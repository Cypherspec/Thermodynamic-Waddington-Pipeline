from __future__ import annotations

import random
from typing import Callable, Sequence

from ..arrays import mean


def permutation_p_value(statistic: Callable[[Sequence[float], Sequence[float]], float], first: Sequence[float], second: Sequence[float], permutations: int = 999, seed: int = 0) -> tuple[float, float]:
    observed = abs(statistic(first, second))
    pooled = list(first) + list(second)
    split = len(first)
    rng = random.Random(seed)
    exceedances = 0
    for _ in range(max(1, permutations)):
        shuffled = list(pooled)
        rng.shuffle(shuffled)
        score = abs(statistic(shuffled[:split], shuffled[split:]))
        if score >= observed:
            exceedances += 1
    return observed, (exceedances + 1) / (max(1, permutations) + 1)


def rank(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    result = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position + 1
        while end < len(order) and values[order[end]] == values[order[position]]:
            end += 1
        value = 0.5 * (position + end - 1) + 1.0
        for index in order[position:end]:
            result[index] = value
        position = end
    return result


def spearman(first: Sequence[float], second: Sequence[float]) -> float:
    if len(first) != len(second) or len(first) < 2:
        return 0.0
    left = rank(first)
    right = rank(second)
    left_mean = mean(left)
    right_mean = mean(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    denominator_left = sum((a - left_mean) ** 2 for a in left)
    denominator_right = sum((b - right_mean) ** 2 for b in right)
    denominator = (denominator_left * denominator_right) ** 0.5
    return numerator / denominator if denominator else 0.0
