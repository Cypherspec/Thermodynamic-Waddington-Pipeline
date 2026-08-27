from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import logsumexp, percentile
from .graph import Edge, adjacency


@dataclass(frozen=True)
class PathEstimate:
    node: int
    free_energy: float
    path_count: int
    effective_sample_size: float
    uncertainty: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "node": self.node,
            "free_energy": self.free_energy,
            "path_count": self.path_count,
            "effective_sample_size": self.effective_sample_size,
            "uncertainty": self.uncertainty,
        }


@dataclass(frozen=True)
class JarzynskiDiagnostic:
    node: int
    path_count: int
    work_mean: float
    work_std: float
    free_energy: float
    dissipation: float
    effective_sample_size: float
    reliability: str

    def to_dict(self) -> dict[str, float | int | str]:
        return {
            "node": self.node,
            "path_count": self.path_count,
            "work_mean": self.work_mean,
            "work_std": self.work_std,
            "free_energy": self.free_energy,
            "dissipation": self.dissipation,
            "effective_sample_size": self.effective_sample_size,
            "reliability": self.reliability,
        }


def jarzynski_diagnostics(path_work: Sequence[Sequence[float]], energies: Sequence[float], temperature: float) -> list[JarzynskiDiagnostic]:
    scale = max(float(temperature), 1e-12)
    records: list[JarzynskiDiagnostic] = []
    for node, values in enumerate(path_work):
        finite = [float(value) for value in values if math.isfinite(float(value))]
        fallback = float(energies[node]) if node < len(energies) else 0.0
        if not finite:
            records.append(JarzynskiDiagnostic(node, 0, 0.0, 0.0, fallback, 0.0, 0.0, "unsupported"))
            continue
        center = sum(finite) / len(finite)
        spread = math.sqrt(sum((value - center) ** 2 for value in finite) / max(1, len(finite)))
        shifted = [-value / scale for value in finite]
        pivot = max(shifted)
        weights = [math.exp(value - pivot) for value in shifted]
        total = sum(weights)
        normalized = [value / max(total, 1e-300) for value in weights]
        ess = 1.0 / max(sum(value * value for value in normalized), 1e-300)
        free_energy = -scale * (pivot + math.log(total / len(finite)))
        ratio = ess / len(finite)
        reliability = "stable" if ratio >= 0.5 and spread <= 3.0 * scale else "high_variance"
        records.append(JarzynskiDiagnostic(node, len(finite), center, spread, free_energy, center - free_energy, ess, reliability))
    return records


def robust_log_mean_exp(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    return logsumexp(values) - math.log(len(values))


def propagate_free_energy(edges: Sequence[Edge], n: int, reference: int, temperature: float, max_paths: int) -> tuple[list[float], list[int], list[list[float]]]:
    outgoing = adjacency(edges, n)
    energies = [float("inf")] * n
    counts = [0] * n
    path_work: list[list[float]] = [[] for _ in range(n)]
    energies[reference] = 0.0
    counts[reference] = 1
    frontier = [reference]
    for _ in range(max(1, n * 2)):
        changed = False
        next_frontier: list[int] = []
        for source in frontier:
            if not math.isfinite(energies[source]):
                continue
            for edge in outgoing[source]:
                if edge.target == reference:
                    continue
                candidates = path_work[edge.target]
                incoming = path_work[source] or [energies[source]]
                for prior in incoming[:max_paths]:
                    if len(candidates) < max_paths:
                        candidates.append(prior + edge.work)
                if candidates:
                    estimate = -temperature * robust_log_mean_exp([-value / temperature for value in candidates])
                    if estimate < energies[edge.target] - 1e-10:
                        energies[edge.target] = estimate
                        counts[edge.target] = len(candidates)
                        next_frontier.append(edge.target)
                        changed = True
        if not changed:
            break
        frontier = list(dict.fromkeys(next_frontier))
    finite = [value for value in energies if math.isfinite(value)]
    fallback = max(finite) if finite else 0.0
    return [value if math.isfinite(value) else fallback for value in energies], counts, path_work


def bootstrap_uncertainty(edges: Sequence[Edge], n: int, reference: int, temperature: float, max_paths: int, replicates: int, seed: int) -> list[float]:
    if replicates <= 1:
        return [0.0] * n
    import random
    rng = random.Random(seed)
    samples = [[] for _ in range(n)]
    for _ in range(replicates):
        chosen = [edge for edge in edges if rng.random() < 0.8] or list(edges)
        values, _, _ = propagate_free_energy(chosen, n, reference, temperature, max_paths)
        for index, value in enumerate(values):
            samples[index].append(value)
    result: list[float] = []
    for values in samples:
        if len(values) < 2:
            result.append(0.0)
            continue
        center = sum(values) / len(values)
        result.append(math.sqrt(sum((value - center) ** 2 for value in values) / (len(values) - 1)))
    return result


def path_estimates(energies: Sequence[float], counts: Sequence[int], uncertainties: Sequence[float], temperature: float) -> list[PathEstimate]:
    scale = max(float(temperature), 1e-12)
    return [PathEstimate(node, float(energy), int(counts[node]), max(1.0, counts[node] / (1.0 + uncertainties[node] / scale)), float(uncertainties[node])) for node, energy in enumerate(energies)]


def barrier_heights(edges: Sequence[Edge], energies: Sequence[float]) -> list[dict[str, float | int]]:
    return [{"source": edge.source, "target": edge.target, "barrier": max(0.0, energies[edge.target] - energies[edge.source]), "work": edge.work, "alignment": edge.alignment} for edge in edges]


def attractor_candidates(edges: Sequence[Edge], energies: Sequence[float], quantile: float = 0.9) -> list[int]:
    outgoing = adjacency(edges, len(energies))
    local = [node for node, energy in enumerate(energies) if outgoing[node] and all(energy <= energies[edge.target] for edge in outgoing[node])]
    if not local:
        return []
    cutoff = percentile([energies[node] for node in local], quantile)
    return [node for node in local if energies[node] <= cutoff]
