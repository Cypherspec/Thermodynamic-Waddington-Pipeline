from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import dot, norm
from .graph import Edge


@dataclass(frozen=True)
class MetricSummary:
    mean_speed: float
    mean_curvature: float
    anisotropy: float
    curl_proxy: float
    divergence_proxy: float
    entropy_rate: float = 0.0
    manifold_dimension: float = 0.0
    tangent_alignment: float = 0.0

    def to_dict(self) -> dict[str, float]:
        return self.__dict__.copy()


def local_curvature(points: Sequence[Sequence[float]], edges: Sequence[Edge]) -> list[float]:
    incident: list[list[int]] = [[] for _ in points]
    for edge in edges:
        incident[edge.source].append(edge.target)
    values: list[float] = []
    for index, point in enumerate(points):
        neighbor_indices = incident[index]
        if not neighbor_indices:
            values.append(0.0)
            continue
        center = [sum(points[target][dimension] for target in neighbor_indices) / len(neighbor_indices) for dimension in range(len(point))]
        displacement = [center[dimension] - point[dimension] for dimension in range(len(point))]
        values.append(norm(displacement) / max(1e-12, 1.0 + norm(point)))
    return values


def graph_entropy(edges: Sequence[Edge], n: int) -> float:
    outgoing = [[] for _ in range(n)]
    for edge in edges:
        outgoing[edge.source].append(max(1e-12, edge.alignment + 1.0))
    values = []
    for weights in outgoing:
        total = sum(weights)
        if total <= 0:
            continue
        values.append(-sum((weight / total) * math.log(weight / total) for weight in weights))
    return sum(values) / len(values) if values else 0.0


def effective_dimension(points: Sequence[Sequence[float]]) -> float:
    if not points:
        return 0.0
    centered = []
    means = [sum(row[d] for row in points) / len(points) for d in range(len(points[0]))]
    for row in points:
        centered.append([value - means[d] for d, value in enumerate(row)])
    variances = [sum(row[d] * row[d] for row in centered) / max(1, len(points) - 1) for d in range(len(means))]
    total = sum(variances)
    if total <= 1e-12:
        return 0.0
    probabilities = [variance / total for variance in variances if variance > 1e-12]
    return math.exp(-sum(probability * math.log(probability) for probability in probabilities))


def summarize(points: Sequence[Sequence[float]], velocities: Sequence[Sequence[float]], edges: Sequence[Edge]) -> MetricSummary:
    speeds = [norm(vector) for vector in velocities]
    curvature = local_curvature(points, edges)
    if not velocities:
        return MetricSummary(0.0, 0.0, 0.0, 0.0, 0.0)
    norms = [max(1e-12, norm(vector)) for vector in velocities]
    anisotropy = max(norms) / min(norms)
    divergence = 0.0
    curl = 0.0
    alignment: list[float] = []
    for edge in edges:
        delta = [velocities[edge.target][d] - velocities[edge.source][d] for d in range(len(velocities[edge.source]))]
        displacement = [points[edge.target][d] - points[edge.source][d] for d in range(len(points[edge.source]))]
        scale = max(1e-12, norm(displacement))
        divergence += dot(delta, displacement) / scale
        alignment.append(edge.alignment)
        if len(delta) >= 2:
            curl += delta[0] * displacement[1] - delta[1] * displacement[0]
    denominator = max(1, len(edges))
    return MetricSummary(sum(speeds) / len(speeds), sum(curvature) / len(curvature), anisotropy, curl / denominator, divergence / denominator, graph_entropy(edges, len(points)), effective_dimension(points), sum(alignment) / len(alignment) if alignment else 0.0)


def project_velocity(velocity: Sequence[float], basis: Sequence[Sequence[float]]) -> list[float]:
    return [dot(velocity, axis) for axis in basis]


def local_jacobian_proxy(points: Sequence[Sequence[float]], velocities: Sequence[Sequence[float]], edges: Sequence[Edge]) -> list[float]:
    values = []
    for edge in edges:
        displacement = [points[edge.target][d] - points[edge.source][d] for d in range(len(points[edge.source]))]
        delta = [velocities[edge.target][d] - velocities[edge.source][d] for d in range(len(velocities[edge.source]))]
        values.append(dot(delta, displacement) / max(1e-12, dot(displacement, displacement)))
    return values


def geodesic_edge_cost(edge: Edge, diffusion: float, curvature: float = 0.0) -> float:
    return edge.distance * (1.0 + curvature) / math.sqrt(max(1e-12, diffusion))
