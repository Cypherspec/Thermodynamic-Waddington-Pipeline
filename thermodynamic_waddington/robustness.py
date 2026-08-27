from __future__ import annotations

"""Multiverse robustness analysis for the effective landscape estimator.

This module does not turn specification sensitivity into biological validation.
It quantifies how much basin ordering, energy ordering, and velocity alignment
change across predeclared analyst choices, then exposes instability as a blocker.
"""

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

from .config import FitConfig
from .model import fit_landscape
from .synthetic import make_synthetic_dataset


@dataclass(frozen=True)
class RobustnessScenario:
    name: str
    neighbors: int
    dimensions: int
    density_bandwidth: float
    diffusion_floor: float
    velocity_scale: float
    edge_alignment_threshold: float
    estimator: str
    energies: tuple[float, ...]
    ordering: tuple[int, ...]
    alignment: float | None
    finite: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _rank(values: Sequence[float]) -> list[int]:
    return [index for index, _ in sorted(enumerate(values), key=lambda item: (float(item[1]), item[0]))]


def _spearman(left: Sequence[float], right: Sequence[float]) -> float:
    n = min(len(left), len(right))
    if n < 2:
        return float("nan")
    rank_left = [0] * n
    rank_right = [0] * n
    for rank, index in enumerate(sorted(range(n), key=lambda i: (float(left[i]), i))):
        rank_left[index] = rank + 1
    for rank, index in enumerate(sorted(range(n), key=lambda i: (float(right[i]), i))):
        rank_right[index] = rank + 1
    ml = sum(rank_left) / n
    mr = sum(rank_right) / n
    numerator = sum((a - ml) * (b - mr) for a, b in zip(rank_left, rank_right))
    denominator = math.sqrt(sum((a - ml) ** 2 for a in rank_left) * sum((b - mr) ** 2 for b in rank_right))
    return numerator / denominator if denominator else 0.0


def _alignment(fit: Any) -> float | None:
    value = (getattr(fit, "diagnostics", {}) or {}).get("velocity_alignment")
    try:
        return float(value) if value is not None and math.isfinite(float(value)) else None
    except (TypeError, ValueError):
        return None


def _scenario_config(base: FitConfig, values: dict[str, Any]) -> FitConfig:
    payload = asdict(base)
    payload.update(values)
    payload["bootstrap_replicates"] = 1
    payload["bootstrap"] = None
    payload["max_paths"] = min(int(payload["max_paths"]), 24)
    return FitConfig(**payload)


def evaluate_robustness(
    expression: Sequence[Sequence[float]],
    velocity: Sequence[Sequence[float]],
    *,
    embedding: Sequence[Sequence[float]] | None = None,
    labels: Sequence[str] | None = None,
    cell_ids: Sequence[str] | None = None,
    base_config: FitConfig | None = None,
    scenarios: Sequence[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    base = base_config or FitConfig(neighbors=min(12, len(expression) - 1), dimensions=min(4, len(expression[0])), bootstrap_replicates=1, max_paths=24)
    default_scenarios = [
        {"name": "baseline"},
        {"name": "local_neighbors_low", "neighbors": max(3, base.neighbors // 2)},
        {"name": "local_neighbors_high", "neighbors": min(len(expression) - 1, max(base.neighbors + 4, base.neighbors * 2))},
        {"name": "diffusion_floor_low", "diffusion_floor": max(1e-8, base.diffusion_floor / 10)},
        {"name": "diffusion_floor_high", "diffusion_floor": base.diffusion_floor * 10},
        {"name": "velocity_half", "velocity_scale": base.velocity_scale * 0.5},
        {"name": "velocity_double", "velocity_scale": base.velocity_scale * 2.0},
        {"name": "alignment_filtered", "edge_alignment_threshold": min(0.9, max(0.2, base.edge_alignment_threshold + 0.2))},
        {"name": "density_only", "estimator": "density_only"},
    ]
    definitions = list(scenarios or default_scenarios)
    results: list[RobustnessScenario] = []
    for offset, definition in enumerate(definitions):
        name = str(definition.get("name", f"scenario_{offset}"))
        values = {key: value for key, value in definition.items() if key != "name"}
        config = _scenario_config(base, values)
        try:
            fit = fit_landscape(expression, velocity, config=config, embedding=embedding, labels=labels, cell_ids=cell_ids, metadata={"robustness_scenario": name, "velocity_status": "observed"})
            energies = tuple(float(value) for value in fit.energies)
            finite = all(math.isfinite(value) for value in energies)
            alignment = _alignment(fit)
            results.append(RobustnessScenario(name, config.neighbors, config.dimensions, config.density_bandwidth, config.diffusion_floor, config.velocity_scale, config.edge_alignment_threshold, config.estimator, energies, tuple(_rank(energies)), alignment, finite))
        except Exception as error:
            results.append(RobustnessScenario(name, config.neighbors, config.dimensions, config.density_bandwidth, config.diffusion_floor, config.velocity_scale, config.edge_alignment_threshold, config.estimator, tuple(), tuple(), None, False))
    baseline = next((item for item in results if item.name == "baseline"), results[0] if results else None)
    correlations = [_spearman(baseline.energies, item.energies) for item in results if baseline and item is not baseline and item.energies]
    finite_fraction = sum(item.finite for item in results) / len(results) if results else 0.0
    min_correlation = min(correlations) if correlations else float("nan")
    alignments = [item.alignment for item in results if item.alignment is not None]
    alignment_range = max(alignments) - min(alignments) if alignments else float("nan")
    blockers: list[str] = []
    if finite_fraction < 1.0:
        blockers.append("nonfinite_or_failed_sensitivity_scenario")
    if correlations and min_correlation < 0.70:
        blockers.append("energy_ordering_is_specification_sensitive")
    if alignments and alignment_range > 0.35:
        blockers.append("velocity_alignment_is_specification_sensitive")
    payload: dict[str, Any] = {
        "suite": "TW-multiverse-robustness-v1",
        "scenarios": [item.to_dict() for item in results],
        "summary": {"scenario_count": len(results), "finite_fraction": finite_fraction, "minimum_energy_spearman": min_correlation, "alignment_range": alignment_range},
        "blockers": blockers,
        "status": "robust_under_declared_specifications" if not blockers else "specification_sensitive_review_required",
        "claim_boundary": "Robustness across analyst specifications is not biological validation, causal evidence, or proof of molecular free energy.",
    }
    payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
    return payload


def run_synthetic_robustness(seed: int = 17, cells: int = 64, genes: int = 8) -> dict[str, Any]:
    dataset = make_synthetic_dataset(cells=cells, genes=genes, seed=seed).to_dict()
    return evaluate_robustness(dataset["expression"], dataset["velocity"], embedding=dataset["embedding"], labels=dataset["labels"], cell_ids=dataset["cell_ids"])


def write_robustness_report(output: str | Path = "experiments/robustness_report.json", **kwargs: Any) -> dict[str, Any]:
    payload = run_synthetic_robustness(**kwargs)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return payload
