from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import dot, mean, norm, solve_linear_system


Matrix = list[list[float]]


@dataclass(frozen=True)
class TensorEstimate:
    tensor: Matrix
    eigenvalues: tuple[float, ...]
    trace: float
    determinant_proxy: float
    condition_proxy: float
    shrinkage: float

    def to_dict(self) -> dict[str, object]:
        return {
            "tensor": self.tensor,
            "eigenvalues": list(self.eigenvalues),
            "trace": self.trace,
            "determinant_proxy": self.determinant_proxy,
            "condition_proxy": self.condition_proxy,
            "shrinkage": self.shrinkage,
        }


def identity(dimensions: int, scale: float = 1.0) -> Matrix:
    return [[scale if row == column else 0.0 for column in range(dimensions)] for row in range(dimensions)]


def outer(vector: Sequence[float]) -> Matrix:
    return [[left * right for right in vector] for left in vector]


def add(first: Matrix, second: Matrix) -> Matrix:
    return [[a + b for a, b in zip(row_a, row_b)] for row_a, row_b in zip(first, second)]


def scale(matrix: Matrix, factor: float) -> Matrix:
    return [[factor * value for value in row] for row in matrix]


def covariance(samples: Sequence[Sequence[float]], floor: float = 1e-8, shrinkage: float = 0.12) -> TensorEstimate:
    if not samples:
        return TensorEstimate([], (), 0.0, 0.0, 0.0, shrinkage)
    dimensions = len(samples[0])
    centers = [mean(sample[column] for sample in samples) for column in range(dimensions)]
    matrix = [[0.0 for _ in range(dimensions)] for _ in range(dimensions)]
    denominator = max(1, len(samples) - 1)
    for sample in samples:
        centered = [value - center for value, center in zip(sample, centers)]
        matrix = add(matrix, outer(centered))
    matrix = scale(matrix, 1.0 / denominator)
    target = identity(dimensions, sum(matrix[column][column] for column in range(dimensions)) / max(1, dimensions))
    matrix = add(scale(matrix, 1.0 - shrinkage), scale(target, shrinkage))
    for index in range(dimensions):
        matrix[index][index] += floor
    eigenvalues = power_eigenvalues(matrix)
    positive = [max(floor, value) for value in eigenvalues]
    determinant = math.prod(positive)
    condition = max(positive) / min(positive)
    return TensorEstimate(matrix, tuple(positive), sum(positive), determinant, condition, shrinkage)


def power_eigenvalues(matrix: Matrix, iterations: int = 80) -> list[float]:
    if not matrix:
        return []
    working = [row[:] for row in matrix]
    values = []
    for component in range(len(matrix)):
        vector = [1.0 if index == component else 0.37 for index in range(len(matrix))]
        for _ in range(iterations):
            product = [sum(working[row][column] * vector[column] for column in range(len(matrix))) for row in range(len(matrix))]
            magnitude = norm(product) or 1.0
            vector = [value / magnitude for value in product]
        product = [sum(working[row][column] * vector[column] for column in range(len(matrix))) for row in range(len(matrix))]
        eigenvalue = dot(vector, product)
        values.append(eigenvalue)
        for row in range(len(matrix)):
            for column in range(len(matrix)):
                working[row][column] -= eigenvalue * vector[row] * vector[column]
    return values


def mahalanobis(vector: Sequence[float], tensor: Matrix, floor: float = 1e-8) -> float:
    regularized = [row[:] for row in tensor]
    for index in range(len(regularized)):
        regularized[index][index] += floor
    solution = solve_linear_system(regularized, list(vector))
    return math.sqrt(max(0.0, dot(vector, solution)))


def tensor_log_determinant(tensor: Matrix, floor: float = 1e-12) -> float:
    return sum(math.log(max(floor, value)) for value in power_eigenvalues(tensor))
