from __future__ import annotations

"""Auditable local nonequilibrium diagnostics for effective fate fields.

The field is a coarse-grained diagnostic in the supplied coordinates. It is not
claimed to be a molecular thermodynamic state function. The implementation
separates conservative-looking gradient structure from rotational probability
current, estimates local entropy-production proxies, and exposes path action
rather than silently calling it physical work.
"""

import math
from dataclasses import asdict, dataclass
from typing import Any, Sequence


@dataclass(frozen=True)
class NonequilibriumConfig:
    regularization: float = 1e-8
    local_radius: float = 1.0
    temperature: float = 1.0
    max_neighbors: int = 12


@dataclass(frozen=True)
class NonequilibriumField:
    conservative: tuple[tuple[float, ...], ...]
    rotational: tuple[tuple[float, ...], ...]
    entropy_proxy: tuple[float, ...]
    divergence: tuple[float, ...]
    circulation: float
    gradient_alignment: float
    status: str
    assumptions: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _dot(left: Sequence[float], right: Sequence[float]) -> float:
    return sum(float(a) * float(b) for a, b in zip(left, right))


def _norm(row: Sequence[float]) -> float:
    return math.sqrt(max(0.0, _dot(row, row)))


def _sub(left: Sequence[float], right: Sequence[float]) -> list[float]:
    return [float(a) - float(b) for a, b in zip(left, right)]


def _scale(row: Sequence[float], factor: float) -> list[float]:
    return [float(x) * factor for x in row]


def _safe_mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _nearest(points: Sequence[Sequence[float]], index: int, limit: int) -> list[int]:
    base = points[index]
    candidates = [(sum((float(a) - float(b)) ** 2 for a, b in zip(base, row)), j) for j, row in enumerate(points) if j != index]
    candidates.sort(key=lambda item: (item[0], item[1]))
    return [j for _, j in candidates[:max(2, limit)]]


def _local_gradient(points: Sequence[Sequence[float]], values: Sequence[float], index: int, neighbors: Sequence[int], regularization: float) -> list[float]:
    dimension = len(points[index])
    gradient = [0.0] * dimension
    weights = 0.0
    for j in neighbors:
        delta = _sub(points[j], points[index])
        distance = _norm(delta)
        if distance <= regularization:
            continue
        weight = 1.0 / (distance + regularization)
        dv = float(values[j]) - float(values[index])
        for k in range(dimension):
            gradient[k] += weight * dv * delta[k] / (distance * distance + regularization)
        weights += weight
    return _scale(gradient, 1.0 / weights) if weights else gradient


def estimate_field(points: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], energies: Sequence[float], config: NonequilibriumConfig | None = None) -> NonequilibriumField:
    cfg = config or NonequilibriumConfig()
    if len(points) != len(velocity) or len(points) != len(energies):
        raise ValueError("points, velocity, and energies must be cell-aligned")
    if not points:
        return NonequilibriumField((), (), (), (), 0.0, 0.0, "empty", ("No cells supplied.",))
    dimension = len(points[0])
    gradients: list[list[float]] = []
    divergences: list[float] = []
    for i in range(len(points)):
        neighbors = _nearest(points, i, cfg.max_neighbors)
        gradients.append(_local_gradient(points, energies, i, neighbors, cfg.regularization))
        if neighbors:
            divergences.append(_safe_mean([_dot(_sub(velocity[j], velocity[i]), _sub(points[j], points[i])) / (sum(x * x for x in _sub(points[j], points[i])) + cfg.regularization) for j in neighbors]))
        else:
            divergences.append(0.0)
    conservative: list[list[float]] = []
    rotational: list[list[float]] = []
    entropy: list[float] = []
    alignment: list[float] = []
    for grad, vel in zip(gradients, velocity):
        force = _scale(grad, -1.0 / max(cfg.temperature, cfg.regularization))
        conservative.append(force)
        residual = _sub(vel, force)
        rotational.append(residual)
        entropy.append(max(0.0, _dot(residual, residual) / max(cfg.temperature, cfg.regularization)))
        alignment.append(_dot(vel, force) / max(cfg.regularization, _norm(vel) * _norm(force)))
    circulation = _safe_mean([_norm(row) for row in rotational])
    gradient_alignment = _safe_mean(alignment)
    assumptions = ("Coordinates are a coarse-grained embedding.", "Velocity is treated as an observed or supplied effective drift.", "Entropy production is a residual norm proxy, not a calibrated physical heat rate.")
    return NonequilibriumField(tuple(tuple(x) for x in conservative), tuple(tuple(x) for x in rotational), tuple(entropy), tuple(divergences), circulation, gradient_alignment, "diagnostic_only", assumptions)


def path_action(path: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], diffusion: float = 1.0, dt: float = 1.0) -> dict[str, Any]:
    if len(path) != len(velocity):
        raise ValueError("path and velocity must have equal length")
    scale = max(float(diffusion), 1e-12)
    increments = []
    for i in range(1, len(path)):
        observed = _scale(_sub(path[i], path[i - 1]), 1.0 / max(dt, 1e-12))
        residual = _sub(observed, velocity[i - 1])
        increments.append(0.5 * _dot(residual, residual) * dt / scale)
    return {"action": sum(increments), "increments": increments, "diffusion": diffusion, "dt": dt, "interpretation": "Onsager-Machlup-style effective path score; not molecular work."}


def compare_reversible_irreversible(points: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], energies: Sequence[float], config: NonequilibriumConfig | None = None) -> dict[str, Any]:
    field = estimate_field(points, velocity, energies, config)
    rotational = [_norm(row) for row in field.rotational]
    conservative = [_norm(row) for row in field.conservative]
    return {"field": field.to_dict(), "rotational_to_conservative": _safe_mean(rotational) / max(1e-12, _safe_mean(conservative)), "reversibility_warning": "A nonzero residual is evidence of model mismatch or nonequilibrium structure, not proof of biological irreversibility.", "status": "diagnostic_only"}
