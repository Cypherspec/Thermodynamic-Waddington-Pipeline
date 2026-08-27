from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Sequence

from .arrays import dot, mean, variance
from .graph import Edge, NeighborGraph


@dataclass(frozen=True)
class PathSample:
    source: int
    target: int
    work: float
    drift_work: float
    density_work: float
    noise_work: float
    displacement: float
    alignment: float
    protocol_weight: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "source": self.source,
            "target": self.target,
            "work": self.work,
            "drift_work": self.drift_work,
            "density_work": self.density_work,
            "noise_work": self.noise_work,
            "displacement": self.displacement,
            "alignment": self.alignment,
            "protocol_weight": self.protocol_weight,
        }


def _safe(value: float, default: float = 0.0) -> float:
    return value if math.isfinite(value) else default


def velocity_residuals(points: Sequence[Sequence[float]], velocities: Sequence[Sequence[float]], graph: NeighborGraph) -> list[list[float]]:
    residuals: list[list[float]] = []
    for index, velocity in enumerate(velocities):
        neighbours = graph.neighbors[index]
        if not neighbours:
            residuals.append([0.0 for _ in velocity])
            continue
        local_mean = [sum(velocities[j][feature] for j in neighbours) / len(neighbours) for feature in range(len(velocity))]
        residuals.append([float(value) - local_mean[feature] for feature, value in enumerate(velocity)])
    return residuals


def estimate_diffusion_tensor(residuals: Sequence[Sequence[float]], graph: NeighborGraph, floor: float) -> list[list[list[float]]]:
    tensors: list[list[list[float]]] = []
    for index, residual in enumerate(residuals):
        neighbours = graph.neighbors[index]
        if not neighbours:
            tensors.append([[floor if row == column else 0.0 for column in range(len(residual))] for row in range(len(residual))])
            continue
        scale = max(floor, sum(value * value for value in residual) / max(1, len(residual)))
        tensors.append([[scale if row == column else 0.0 for column in range(len(residual))] for row in range(len(residual))])
    return tensors


def tensor_summary(tensors: Sequence[Sequence[Sequence[float]]]) -> dict[str, object]:
    traces = [sum(row[index] for index, row in enumerate(tensor) if index < len(row)) for tensor in tensors]
    if not traces:
        return {"cells": 0, "trace_mean": 0.0, "trace_min": 0.0, "trace_max": 0.0, "model": "local diagonal residual covariance proxy"}
    return {"cells": len(traces), "trace_mean": mean(traces), "trace_min": min(traces), "trace_max": max(traces), "model": "local diagonal residual covariance proxy", "warning": "replicate-aware diffusion calibration is required for physical interpretation"}


def edge_path_samples(edges: Sequence[Edge], temperature: float) -> list[PathSample]:
    scale = max(temperature, 1e-12)
    return [PathSample(edge.source, edge.target, edge.work, edge.work * edge.alignment, edge.work * (1.0 - edge.alignment), max(0.0, edge.work * edge.work / scale), edge.distance, edge.alignment, math.exp(max(-50.0, min(50.0, -edge.work / scale)))) for edge in edges]


def path_sample_summary(samples: Sequence[PathSample], temperature: float) -> dict[str, object]:
    if not samples:
        return {"count": 0, "temperature": temperature, "work_mean": 0.0, "work_variance": 0.0, "effective_sample_size": 0.0}
    works = [sample.work for sample in samples]
    weights = [sample.protocol_weight for sample in samples]
    total = sum(weights)
    ess = total * total / max(1e-12, sum(weight * weight for weight in weights))
    return {"count": len(samples), "temperature": temperature, "work_mean": mean(works), "work_variance": variance(works), "work_min": min(works), "work_max": max(works), "effective_sample_size": ess, "negative_work_fraction": sum(1 for value in works if value < 0) / len(works), "components": {"drift_mean": mean([sample.drift_work for sample in samples]), "density_mean": mean([sample.density_work for sample in samples]), "noise_mean": mean([sample.noise_work for sample in samples])}}


def protocol_sanity_checks(samples: Sequence[PathSample]) -> dict[str, object]:
    if not samples:
        return {"status": "insufficient_data", "checks": []}
    works = [sample.work for sample in samples]
    alignments = [sample.alignment for sample in samples]
    finite = all(math.isfinite(value) for value in works)
    return {"status": "diagnostic_only", "finite_work": finite, "has_forward_and_reverse_work": any(value < 0 for value in works) and any(value >= 0 for value in works), "alignment_mean": mean(alignments), "protocol_defined": False, "reason": "an observational velocity field is not a controlled nonequilibrium protocol; Jarzynski equality is therefore not identified without additional experimental design"}


def bootstrap_path_work(edges: Sequence[Edge], replicates: int, seed: int) -> dict[str, object]:
    if not edges or replicates < 1:
        return {"replicates": 0, "work_mean_interval": [0.0, 0.0], "work_mean": 0.0}
    rng = random.Random(seed)
    estimates: list[float] = []
    for _ in range(replicates):
        sample = [edge.work for edge in edges if rng.random() < 0.8] or [edge.work for edge in edges]
        estimates.append(mean(sample))
    ordered = sorted(estimates)
    return {"replicates": replicates, "work_mean": mean(estimates), "work_mean_interval": [ordered[max(0, int(0.025 * len(ordered)) - 1)], ordered[min(len(ordered) - 1, int(0.975 * len(ordered)))]]}


def local_current(points: Sequence[Sequence[float]], velocities: Sequence[Sequence[float]], graph: NeighborGraph, densities: Sequence[float], diffusions: Sequence[float]) -> list[float]:
    currents: list[float] = []
    for index, velocity in enumerate(velocities):
        neighbours = graph.neighbors[index]
        if not neighbours:
            currents.append(0.0)
            continue
        divergence = 0.0
        for target in neighbours:
            displacement = [points[target][feature] - points[index][feature] for feature in range(len(points[index]))]
            distance_sq = max(1e-12, dot(displacement, displacement))
            neighbour_density = densities[target]
            local_density = max(1e-12, densities[index])
            drift = dot(velocity[:len(displacement)], displacement) / math.sqrt(distance_sq)
            diffusion_gradient = (neighbour_density - local_density) / math.sqrt(distance_sq)
            divergence += drift * local_density - diffusions[index] * diffusion_gradient
        currents.append(divergence / len(neighbours))
    return currents


def current_summary(currents: Sequence[float]) -> dict[str, object]:
    if not currents:
        return {"cells": 0, "mean_abs_current": 0.0, "max_abs_current": 0.0}
    absolute = [abs(_safe(value)) for value in currents]
    return {"cells": len(currents), "mean_current": mean(currents), "mean_abs_current": mean(absolute), "max_abs_current": max(absolute), "nonzero_fraction": sum(1 for value in absolute if value > 1e-12) / len(absolute), "interpretation": "probability-current proxy; active transcription and sampling effects are not separated"}
