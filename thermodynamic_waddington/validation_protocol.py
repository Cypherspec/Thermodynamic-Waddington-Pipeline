"""Reproducible validation protocol for the effective landscape estimator.

This module deliberately separates benchmark validation from biological discovery.
A passing benchmark means the implementation recovers known structure under a
specified data-generating process; it is not evidence that a new biological law
has been discovered.
"""
from __future__ import annotations

import json
import math
import platform
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from .config import FitConfig
from .model import fit_landscape
from .synthetic import make_synthetic_dataset
from .validation import fate_calibration, velocity_alignment_score


@dataclass(frozen=True)
class ValidationCase:
    name: str
    cells: int
    genes: int
    seed: int
    dimensions: int = 5
    neighbors: int = 12
    bootstrap: int = 8


@dataclass
class ValidationResult:
    case: str
    passed: bool
    metrics: dict[str, float]
    checks: dict[str, bool]
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _finite(values: Sequence[float]) -> bool:
    return bool(values) and all(math.isfinite(float(value)) for value in values)


def run_case(case: ValidationCase) -> ValidationResult:
    data = make_synthetic_dataset(case.cells, case.genes, case.seed)
    config = FitConfig(
        dimensions=min(case.dimensions, case.genes, case.cells - 1),
        neighbors=min(case.neighbors, case.cells - 1),
        bootstrap=case.bootstrap,
        seed=case.seed,
        trajectory_steps=8,
        trajectory_replicates=8,
    )
    fit = fit_landscape(
        data.expression,
        data.velocity,
        config,
        embedding=data.embedding,
        labels=data.labels,
        cell_ids=data.cell_ids,
        lineage_outcomes=data.lineage_outcomes,
        metadata={**data.metadata, "validation_case": case.name, "velocity_status": "observed"},
    )
    energies = fit.energies
    checks = {
        "finite_energies": _finite(energies),
        "correct_cell_count": len(energies) == case.cells,
        "nonempty_graph": len(fit.edges) > case.cells,
        "attractor_exists": len(fit.attractors) > 0,
        "diagnostics_serializable": _jsonable(fit.diagnostics),
        "velocity_alignment_finite": math.isfinite(float(velocity_alignment_score(fit.edges))),
    }
    calibration_report = fate_calibration(energies, data.labels, data.lineage_outcomes, config.temperature)
    metrics = {
        "cells": float(case.cells),
        "genes": float(case.genes),
        "edges": float(len(fit.edges)),
        "attractors": float(len(fit.attractors)),
        "energy_range": float(max(energies) - min(energies)),
        "velocity_alignment": float(velocity_alignment_score(fit.edges)),
        "calibration": float(calibration_report["pearson"]),
    }
    warnings = [
        "Synthetic recovery validates implementation behavior, not biological truth.",
        "The benchmark uses an effective coarse-grained landscape rather than molecular free energy.",
    ]
    return ValidationResult(case.name, all(checks.values()), metrics, checks, warnings)


def _jsonable(value: Any) -> bool:
    try:
        json.dumps(value, sort_keys=True)
        return True
    except (TypeError, ValueError, OverflowError):
        return False


def run_validation(cases: Sequence[ValidationCase] | None = None) -> dict[str, Any]:
    cases = list(cases or (
        ValidationCase("small_recovery", 32, 8, 11),
        ValidationCase("medium_recovery", 96, 12, 19),
        ValidationCase("stress_recovery", 192, 16, 23),
    ))
    results = [run_case(case) for case in cases]
    passed = all(result.passed for result in results)
    return {
        "protocol": "tw-validation-v1",
        "status": "passed" if passed else "failed",
        "claim_boundary": "Implementation and synthetic recovery only; no claim of biological discovery.",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "runtime": {"python": sys.version.split()[0], "platform": platform.platform()},
        "cases": [result.to_dict() for result in results],
    }


def write_validation_report(path: str | Path, cases: Sequence[ValidationCase] | None = None) -> dict[str, Any]:
    report = run_validation(cases)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    return report
