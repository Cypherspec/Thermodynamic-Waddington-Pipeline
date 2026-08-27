from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Iterable, Sequence


def shape(matrix: Sequence[Sequence[float]]) -> tuple[int, int]:
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    if any(len(row) != cols for row in matrix):
        raise ValueError("matrix rows have inconsistent lengths")
    return rows, cols


def zeros(rows: int, cols: int) -> list[list[float]]:
    return [[0.0 for _ in range(cols)] for _ in range(rows)]


def copy_matrix(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    return [list(map(float, row)) for row in matrix]


def transpose(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    rows, cols = shape(matrix)
    return [[float(matrix[i][j]) for i in range(rows)] for j in range(cols)]


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def norm(a: Sequence[float]) -> float:
    return math.sqrt(max(0.0, dot(a, a)))


def add(a: Sequence[float], b: Sequence[float]) -> list[float]:
    return [x + y for x, y in zip(a, b)]


def sub(a: Sequence[float], b: Sequence[float]) -> list[float]:
    return [x - y for x, y in zip(a, b)]


def scale(a: Sequence[float], factor: float) -> list[float]:
    return [factor * x for x in a]


def mean(values: Iterable[float]) -> float:
    values = list(values)
    return sum(values) / len(values) if values else 0.0


def variance(values: Iterable[float]) -> float:
    values = list(values)
    if len(values) < 2:
        return 0.0
    center = mean(values)
    return sum((value - center) ** 2 for value in values) / (len(values) - 1)


def row_means(matrix: Sequence[Sequence[float]]) -> list[float]:
    return [mean(row) for row in matrix]


def column_means(matrix: Sequence[Sequence[float]]) -> list[float]:
    return row_means(transpose(matrix))


def center(matrix: Sequence[Sequence[float]]) -> tuple[list[list[float]], list[float]]:
    means = column_means(matrix)
    return [[value - means[j] for j, value in enumerate(row)] for row in matrix], means


def pairwise_squared_distances(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    count, _ = shape(matrix)
    distances = zeros(count, count)
    for i in range(count):
        for j in range(i + 1, count):
            distance = sum((a - b) ** 2 for a, b in zip(matrix[i], matrix[j]))
            distances[i][j] = distance
            distances[j][i] = distance
    return distances


def median(values: Iterable[float]) -> float:
    values = sorted(values)
    if not values:
        return 0.0
    middle = len(values) // 2
    if len(values) % 2:
        return values[middle]
    return 0.5 * (values[middle - 1] + values[middle])


def percentile(values: Iterable[float], probability: float) -> float:
    values = sorted(values)
    if not values:
        return 0.0
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    weight = position - lower
    return values[lower] * (1 - weight) + values[upper] * weight


def quantile(values: Iterable[float], probability: float) -> float:
    return percentile(values, probability)


def logsumexp(values: Iterable[float]) -> float:
    values = list(values)
    if not values:
        return float("-inf")
    maximum = max(values)
    if math.isinf(maximum):
        return maximum
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def safe_log(value: float, floor: float = 1e-12) -> float:
    return math.log(max(floor, value))


def safe_exp(value: float, ceiling: float = 60.0) -> float:
    return math.exp(max(-ceiling, min(ceiling, value)))


def normalize(values: Sequence[float]) -> list[float]:
    total = sum(values)
    if total <= 0:
        return [1.0 / len(values) for _ in values] if values else []
    return [value / total for value in values]


def read_matrix(path: str | Path) -> list[list[float]]:
    path = Path(path)
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text())
        if isinstance(payload, dict):
            payload = payload.get("data", payload.get("matrix"))
        return copy_matrix(payload)
    if path.suffix.lower() == ".jsonl":
        return [[float(value) for value in json.loads(line)] for line in path.read_text().splitlines() if line.strip()]
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open(newline="") as handle:
        return [[float(value) for value in row] for row in csv.reader(handle, delimiter=delimiter) if row]


def write_matrix(path: str | Path, matrix: Sequence[Sequence[float]]) -> None:
    path = Path(path)
    if path.suffix.lower() == ".json":
        path.write_text(json.dumps(matrix, separators=(",", ":")))
        return
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open("w", newline="") as handle:
        csv.writer(handle, delimiter=delimiter).writerows(matrix)


def solve_linear_system(matrix: Sequence[Sequence[float]], vector: Sequence[float], floor: float = 1e-12) -> list[float]:
    size = len(vector)
    augmented = [list(map(float, matrix[row])) + [float(vector[row])] for row in range(size)]
    for pivot in range(size):
        pivot_row = max(range(pivot, size), key=lambda row: abs(augmented[row][pivot]))
        if abs(augmented[pivot_row][pivot]) < floor:
            augmented[pivot_row][pivot] = floor
        augmented[pivot], augmented[pivot_row] = augmented[pivot_row], augmented[pivot]
        divisor = augmented[pivot][pivot]
        augmented[pivot] = [value / divisor for value in augmented[pivot]]
        for row in range(size):
            if row == pivot:
                continue
            factor = augmented[row][pivot]
            augmented[row] = [left - factor * right for left, right in zip(augmented[row], augmented[pivot])]
    return [augmented[row][-1] for row in range(size)]


def matrix_vector(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [dot(row, vector) for row in matrix]


def matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    """Multiply a rectangular matrix by a vector."""
    if not matrix:
        return []
    return [sum(float(value) * float(weight) for value, weight in zip(row, vector)) for row in matrix]


def matrix_multiply(first: Sequence[Sequence[float]], second: Sequence[Sequence[float]]) -> list[list[float]]:
    right = transpose(second)
    return [[dot(row, column) for column in right] for row in first]


def trace(matrix: Sequence[Sequence[float]]) -> float:
    return sum(matrix[index][index] for index in range(min(len(matrix), len(matrix[0]) if matrix else 0)))
