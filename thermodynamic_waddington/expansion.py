from __future__ import annotations

import hashlib
import json
import math
import random
from dataclasses import dataclass, asdict
from typing import Callable, Iterable, Sequence

from .arrays import dot, mean, norm, percentile, variance
from .graph import Edge


@dataclass(frozen=True)
class UncertaintyEnvelope:
    center: float
    lower: float
    upper: float
    standard_error: float
    effective_samples: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class Counterfactual:
    source: int
    target: int
    cost: float
    probability_gain: float
    control_norm: float
    mechanism: str

    def to_dict(self) -> dict[str, float | int | str]:
        return asdict(self)


def stable_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    return hashlib.sha256(payload).hexdigest()[:24]


def weighted_mean(values: Sequence[float], weights: Sequence[float]) -> float:
    if len(values) != len(weights) or not values:
        return 0.0
    total = sum(max(0.0, float(w)) for w in weights)
    return sum(float(v) * max(0.0, float(w)) for v, w in zip(values, weights)) / total if total else mean(values)


def weighted_quantile(values: Sequence[float], weights: Sequence[float], q: float) -> float:
    if not values or len(values) != len(weights):
        return 0.0
    pairs = sorted((float(v), max(0.0, float(w))) for v, w in zip(values, weights))
    total = sum(weight for _, weight in pairs)
    if total <= 0:
        return percentile([value for value, _ in pairs], q)
    threshold = max(0.0, min(1.0, q)) * total
    cumulative = 0.0
    for value, weight in pairs:
        cumulative += weight
        if cumulative >= threshold:
            return value
    return pairs[-1][0]


def jackknife_envelope(values: Sequence[float]) -> UncertaintyEnvelope:
    clean = [float(value) for value in values if math.isfinite(float(value))]
    if not clean:
        return UncertaintyEnvelope(0.0, 0.0, 0.0, 0.0, 0.0)
    center = mean(clean)
    if len(clean) == 1:
        return UncertaintyEnvelope(center, center, center, 0.0, 1.0)
    leave_one_out = [mean(clean[:i] + clean[i + 1:]) for i in range(len(clean))]
    pseudo_variance = (len(clean) - 1) / len(clean) * sum((value - mean(leave_one_out)) ** 2 for value in leave_one_out)
    error = math.sqrt(max(0.0, pseudo_variance / len(clean)))
    return UncertaintyEnvelope(center, center - 1.96 * error, center + 1.96 * error, error, float(len(clean)))


def bootstrap_envelope(values: Sequence[float], replicates: int = 256, seed: int = 17, statistic: Callable[[list[float]], float] = mean) -> UncertaintyEnvelope:
    clean = [float(value) for value in values if math.isfinite(float(value))]
    if not clean:
        return UncertaintyEnvelope(0.0, 0.0, 0.0, 0.0, 0.0)
    rng = random.Random(seed)
    samples = [float(statistic([clean[rng.randrange(len(clean))] for _ in clean])) for _ in range(max(1, replicates))]
    center = float(statistic(clean))
    spread = math.sqrt(max(0.0, variance(samples)))
    return UncertaintyEnvelope(center, percentile(samples, 0.025), percentile(samples, 0.975), spread, float(len(clean)))


def edge_work_distribution(edges: Sequence[Edge], temperature: float = 1.0) -> dict[str, object]:
    works = [float(edge.work) for edge in edges if math.isfinite(float(edge.work))]
    envelope = bootstrap_envelope(works, seed=temperature.__hash__() & 0xffff)
    return {"count": len(works), "mean": mean(works), "variance": variance(works), "envelope": envelope.to_dict()}


def current_balance(edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0) -> dict[str, float]:
    if not edges:
        return {"forward_mass": 0.0, "reverse_mass": 0.0, "net_bias": 0.0, "entropy_production_proxy": 0.0}
    forward = 0.0
    reverse = 0.0
    for edge in edges:
        delta = (float(energies[edge.target]) - float(energies[edge.source])) / max(1e-12, temperature)
        forward += math.exp(max(-40.0, min(40.0, -delta)))
        reverse += math.exp(max(-40.0, min(40.0, delta)))
    total = forward + reverse
    bias = (forward - reverse) / total if total else 0.0
    return {"forward_mass": forward, "reverse_mass": reverse, "net_bias": bias, "entropy_production_proxy": bias * bias * total / max(1, len(edges))}


def committor_iteration(edges: Sequence[Edge], energies: Sequence[float], source_set: set[int], target_set: set[int], temperature: float = 1.0, iterations: int = 256, tolerance: float = 1e-9) -> list[float]:
    n = len(energies)
    outgoing: list[list[Edge]] = [[] for _ in range(n)]
    for edge in edges:
        outgoing[edge.source].append(edge)
    values = [0.0 if i in source_set else 1.0 if i in target_set else 0.5 for i in range(n)]
    for _ in range(max(1, iterations)):
        next_values = list(values)
        residual = 0.0
        for source in range(n):
            if source in source_set or source in target_set or not outgoing[source]:
                continue
            weights = [math.exp(max(-40.0, min(40.0, -(energies[edge.target] - energies[source]) / max(temperature, 1e-12)))) for edge in outgoing[source]]
            total = sum(weights) or 1.0
            estimate = sum(weight * values[edge.target] for weight, edge in zip(weights, outgoing[source])) / total
            next_values[source] = estimate
            residual = max(residual, abs(estimate - values[source]))
        values = next_values
        if residual <= tolerance:
            break
    return values


def reactive_flux(edges: Sequence[Edge], committor: Sequence[float], energies: Sequence[float], temperature: float = 1.0) -> list[dict[str, float | int]]:
    result = []
    for edge in edges:
        potential = math.exp(max(-40.0, min(40.0, -(energies[edge.target] - energies[edge.source]) / max(temperature, 1e-12))))
        delta_q = float(committor[edge.target]) - float(committor[edge.source])
        result.append({"source": edge.source, "target": edge.target, "flux": potential * max(0.0, delta_q), "committor_delta": delta_q, "alignment": edge.alignment})
    return sorted(result, key=lambda item: float(item["flux"]), reverse=True)


def minimum_action_path(edges: Sequence[Edge], energies: Sequence[float], source: int, target: int, temperature: float = 1.0) -> dict[str, object]:
    n = len(energies)
    outgoing: list[list[Edge]] = [[] for _ in range(n)]
    for edge in edges:
        outgoing[edge.source].append(edge)
    distances = [float("inf")] * n
    previous: list[int | None] = [None] * n
    distances[source] = 0.0
    remaining = set(range(n))
    while remaining:
        node = min(remaining, key=lambda index: distances[index])
        remaining.remove(node)
        if not math.isfinite(distances[node]) or node == target:
            break
        for edge in outgoing[node]:
            action = max(0.0, float(edge.work)) + max(0.0, float(energies[edge.target]) - float(energies[node])) / max(temperature, 1e-12)
            candidate = distances[node] + action
            if candidate < distances[edge.target]:
                distances[edge.target] = candidate
                previous[edge.target] = node
    if not math.isfinite(distances[target]):
        return {"source": source, "target": target, "path": [], "action": float("inf"), "reachable": False}
    path = [target]
    while path[-1] != source and previous[path[-1]] is not None:
        path.append(previous[path[-1]])
    path.reverse()
    return {"source": source, "target": target, "path": path, "action": distances[target], "reachable": path[0] == source}


def attractor_transition_matrix(edges: Sequence[Edge], energies: Sequence[float], attractors: Sequence[int], temperature: float = 1.0) -> list[list[float]]:
    result = [[0.0 for _ in attractors] for _ in attractors]
    if not attractors:
        return result
    for i, source in enumerate(attractors):
        for j, target in enumerate(attractors):
            if i != j:
                barrier = max(0.0, float(energies[target]) - float(energies[source]))
                result[i][j] = math.exp(max(-40.0, min(40.0, -barrier / max(temperature, 1e-12))))
        total = sum(result[i])
        if total:
            result[i] = [value / total for value in result[i]]
    return result


def graph_signature(edges: Sequence[Edge]) -> str:
    return stable_hash([{"source": edge.source, "target": edge.target, "work": round(float(edge.work), 12)} for edge in edges])


def summarize_frontier(values: Sequence[float], edges: Sequence[Edge]) -> dict[str, object]:
    by_node: dict[int, list[float]] = {}
    for edge in edges:
        by_node.setdefault(edge.source, []).append(float(values[edge.target]) - float(values[edge.source]))
    return {"nodes": len(by_node), "positive_fraction": sum(1 for delta in (d for values_ in by_node.values() for d in values_) if delta > 0) / max(1, sum(len(values_) for values_ in by_node.values())), "node_ranges": {str(node): [min(deltas), max(deltas)] for node, deltas in sorted(by_node.items())}}


def adaptive_temperature(energies: Sequence[float], target_acceptance: float = 0.234) -> float:
    if not energies:
        return 1.0
    scale = math.sqrt(max(1e-12, variance(energies)))
    return max(1e-6, scale * (1.0 - 0.5 * max(0.0, min(1.0, target_acceptance))))

def local_curvature(energies: Sequence[float], edges: Sequence[Edge]) -> list[float]:
    neighborhoods: dict[int, list[float]] = {}
    for edge in edges:
        neighborhoods.setdefault(edge.source, []).append(float(energies[edge.target]) - float(energies[edge.source]))
    return [variance(neighborhoods.get(i, [0.0])) for i in range(len(energies))]

def basin_entropy(labels: Sequence[int]) -> float:
    counts: dict[int, int] = {}
    for label in labels:
        counts[int(label)] = counts.get(int(label), 0) + 1
    total = sum(counts.values())
    return -sum((count / total) * math.log(count / total) for count in counts.values() if count and total)

def spectral_relaxation_proxy(matrix: Sequence[Sequence[float]], iterations: int = 256) -> dict[str, float]:
    if not matrix:
        return {"gap": 0.0, "relaxation_time": float("inf"), "residual": 0.0}
    n = len(matrix)
    state = [1.0 / n] * n
    for _ in range(max(1, iterations)):
        next_state = [sum(state[i] * matrix[i][j] for i in range(n)) for j in range(n)]
        total = sum(next_state) or 1.0
        next_state = [value / total for value in next_state]
        residual = sum(abs(a - b) for a, b in zip(next_state, state))
        state = next_state
        if residual < 1e-10:
            break
    diagonal = mean([matrix[i][i] for i in range(n)])
    gap = max(0.0, 1.0 - diagonal)
    return {"gap": gap, "relaxation_time": 1.0 / max(gap, 1e-12), "residual": residual}

def sensitivity_scan(values: Sequence[float], perturbations: Sequence[float] = (-0.2, -0.1, 0.1, 0.2)) -> list[dict[str, float]]:
    center = mean(values) if values else 0.0
    spread = math.sqrt(max(0.0, variance(values))) if values else 0.0
    return [{"perturbation": float(delta), "mean_shift": delta * center, "spread_shift": delta * spread, "relative_response": abs(delta) / max(1e-12, 1.0 + abs(center))} for delta in perturbations]

def path_action(path: Sequence[int], energies: Sequence[float], edges: Sequence[Edge], temperature: float = 1.0) -> float:
    lookup = {(edge.source, edge.target): edge for edge in edges}
    total = 0.0
    for source, target in zip(path, path[1:]):
        edge = lookup.get((source, target))
        if edge is None:
            total += abs(float(energies[target]) - float(energies[source]))
        else:
            total += max(0.0, float(edge.work)) + max(0.0, float(energies[target]) - float(energies[source])) / max(temperature, 1e-12)
    return total

def compare_landscapes(first: Sequence[float], second: Sequence[float]) -> dict[str, float]:
    if len(first) != len(second) or not first:
        return {"correlation": 0.0, "rmse": float("inf"), "mean_shift": float("nan")}
    a, b = mean(first), mean(second)
    numerator = sum((x - a) * (y - b) for x, y in zip(first, second))
    denominator = math.sqrt(sum((x - a) ** 2 for x in first) * sum((y - b) ** 2 for y in second))
    rmse = math.sqrt(mean([(x - y) ** 2 for x, y in zip(first, second)]))
    return {"correlation": numerator / denominator if denominator else 0.0, "rmse": rmse, "mean_shift": b - a}
