from __future__ import annotations

import math
import random
from typing import Sequence

from .arrays import dot, norm, zeros


def matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [dot(row, vector) for row in matrix]


def covariance(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    if not matrix:
        return []
    n = len(matrix)
    p = len(matrix[0])
    means = [sum(row[j] for row in matrix) / n for j in range(p)]
    result = zeros(p, p)
    divisor = max(1, n - 1)
    for row in matrix:
        centered = [row[j] - means[j] for j in range(p)]
        for i in range(p):
            for j in range(i, p):
                result[i][j] += centered[i] * centered[j] / divisor
    for i in range(p):
        for j in range(i):
            result[i][j] = result[j][i]
    return result


def rayleigh(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> float:
    transformed = matvec(matrix, vector)
    denominator = dot(vector, vector)
    return dot(vector, transformed) / denominator if denominator else 0.0


def normalize_vector(vector: Sequence[float], floor: float = 1e-14) -> list[float]:
    length = norm(vector)
    if length < floor:
        return [0.0 for _ in vector]
    return [value / length for value in vector]


def power_eigenvector(matrix: Sequence[Sequence[float]], seed: int, iterations: int = 250) -> tuple[float, list[float]]:
    size = len(matrix)
    rng = random.Random(seed)
    vector = normalize_vector([rng.uniform(-1.0, 1.0) for _ in range(size)])
    for _ in range(iterations):
        vector = normalize_vector(matvec(matrix, vector))
    return rayleigh(matrix, vector), vector


def deflate(matrix: Sequence[Sequence[float]], eigenvalue: float, eigenvector: Sequence[float]) -> list[list[float]]:
    result = [list(row) for row in matrix]
    for i in range(len(result)):
        for j in range(len(result)):
            result[i][j] -= eigenvalue * eigenvector[i] * eigenvector[j]
    return result


def pca(matrix: Sequence[Sequence[float]], components: int, seed: int) -> tuple[list[list[float]], list[float], list[list[float]]]:
    if not matrix:
        return [], [], []
    means = [sum(row[j] for row in matrix) / len(matrix) for j in range(len(matrix[0]))]
    centered = [[value - means[j] for j, value in enumerate(row)] for row in matrix]
    cov = covariance(matrix)
    eigenvalues: list[float] = []
    eigenvectors: list[list[float]] = []
    residual = cov
    for component in range(min(components, len(cov))):
        value, vector = power_eigenvector(residual, seed + component * 101)
        value = max(0.0, value)
        eigenvalues.append(value)
        eigenvectors.append(vector)
        residual = deflate(residual, value, vector)
    embedding = [[dot(row, vector) for vector in eigenvectors] for row in centered]
    return embedding, eigenvalues, eigenvectors


def project(rows: Sequence[Sequence[float]], means: Sequence[float], vectors: Sequence[Sequence[float]]) -> list[list[float]]:
    return [[dot([value - means[j] for j, value in enumerate(row)], vector) for vector in vectors] for row in rows]


def two_dimensional_projection(matrix: Sequence[Sequence[float]], seed: int) -> list[list[float]]:
    embedding, _, _ = pca(matrix, 2, seed)
    return [row + [0.0] * (2 - len(row)) for row in embedding]
