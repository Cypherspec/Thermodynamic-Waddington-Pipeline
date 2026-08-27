from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

from .config import FitConfig
from .larry_validation import validate_larry_fate_ground_truth
from .model import fit_landscape
from .real_data import benchmark_provenance, derive_velocity_from_spliced_unspliced


@dataclass(frozen=True)
class BenchmarkResult:
    name: str
    status: str
    cells: int
    features: int
    metrics: dict[str, float]
    provenance: dict[str, Any]
    limitations: list[str]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _finite(values: Sequence[float]) -> list[float]:
    return [float(value) for value in values if math.isfinite(float(value))]


def _metric(name: str, fit: Any) -> dict[str, float]:
    energies = _finite(fit.energies)
    alignment = fit.diagnostics.get("velocity_alignment")
    return {
        "cells": float(len(energies)),
        "features": float(len(fit.metadata.get("gene_names", [])) or 0),
        "edges": float(len(fit.edges)),
        "attractors": float(len(fit.attractors)),
        "finite_energy_fraction": len(energies) / max(1, len(fit.energies)),
        "energy_range": (max(energies) - min(energies)) if energies else float("nan"),
        "velocity_alignment": float(alignment) if isinstance(alignment, (int, float)) else float("nan"),
    }


def run_pancreas_benchmark(path: str | Path = "data/real/endocrinogenesis_day15.h5ad", max_cells: int = 128, max_genes: int = 64, label_key: str = "clusters_coarse", seed: int = 19) -> BenchmarkResult:
    source = Path(path)
    bundle = derive_velocity_from_spliced_unspliced(source, max_cells=max_cells, max_genes=max_genes, label_key=label_key)
    config = FitConfig(neighbors=min(12, max_cells - 1), dimensions=min(6, max_genes), bootstrap=3, seed=seed, trajectory_steps=5, trajectory_replicates=8)
    fit = fit_landscape(bundle["expression"], bundle["velocity"], config, labels=bundle["labels"], cell_ids=bundle["cell_ids"], metadata=bundle["metadata"])
    provenance = benchmark_provenance(source, "scvelo_pancreas_proxy", "https://github.com/theislab/scvelo_notebooks/raw/master/data/Pancreas/endocrinogenesis_day15.h5ad", False, expression_key="spliced", velocity_key="derived:unspliced-minus-spliced", label_key=label_key)
    return BenchmarkResult("pancreas_proxy", "passed_with_proxy_velocity", len(bundle["expression"]), len(bundle["expression"][0]), _metric("pancreas_proxy", fit), provenance, ["The file has spliced/unspliced layers but no precomputed velocity layer.", "Velocity is the transparent unspliced-minus-spliced proxy, not a scVelo dynamical estimate.", "This validates execution, provenance, and numerical behavior, not biological truth."])


def run_larry_ground_truth_benchmark(root: str | Path = "data/real/larry") -> BenchmarkResult:
    report = validate_larry_fate_ground_truth(root)
    metrics = {"cells_with_fate_bias": float(report["cells_with_fate_bias"]), "early_cells": float(report["early_cell_count"]), "fates": float(len(report["fates"]))}
    return BenchmarkResult("larry_fate_ground_truth", "loaded_reference_only", int(metrics["cells_with_fate_bias"]), int(metrics["fates"]), metrics, {"source": report["source"], "ground_truth_sha256": report["ground_truth_sha256"]}, list(report["limitations"]))


def run_external_validation_report(output: str | Path = "experiments/external_validation.json", pancreas_path: str | Path = "data/real/endocrinogenesis_day15.h5ad", larry_root: str | Path = "data/real/larry") -> dict[str, Any]:
    results = [run_pancreas_benchmark(pancreas_path).to_dict(), run_larry_ground_truth_benchmark(larry_root).to_dict()]
    report = {"protocol": "external_validation_v1", "status": "passed_with_explicit_boundaries", "results": results, "claim": "The implementation is reproducible against public benchmark artifacts and a public lineage-fate reference. It does not establish a molecular free-energy law or causal intervention effect."}
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    return report
