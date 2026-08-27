from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Any, Sequence

from .arrays import mean, percentile
from .model import LandscapeFit


@dataclass(frozen=True)
class VisualDiagnostic:
    key: str
    value: float
    status: str
    explanation: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _status(value: float, good: float, warning: float, inverse: bool = False) -> str:
    if inverse:
        if value <= good:
            return "good"
        if value <= warning:
            return "warning"
        return "critical"
    if value >= good:
        return "good"
    if value >= warning:
        return "warning"
    return "critical"


def diagnose_fit(fit: LandscapeFit) -> list[VisualDiagnostic]:
    uncertainty = mean(fit.uncertainties) if fit.uncertainties else 0.0
    coverage = float(fit.diagnostics.get("coverage", 0.0))
    alignment = float(fit.diagnostics.get("mean_velocity_alignment", 0.0) or 0.0)
    barriers = [float(item.get("barrier", 0.0)) for item in fit.barriers if isinstance(item, dict)]
    return [
        VisualDiagnostic("coverage", coverage, _status(coverage, 0.9, 0.65), "Fraction of cells reachable from the reference state."),
        VisualDiagnostic("velocity_alignment", alignment, _status(alignment, 0.65, 0.35), "Agreement between observed velocity and inferred directed graph edges."),
        VisualDiagnostic("mean_energy_uncertainty", uncertainty, _status(uncertainty, percentile(fit.energies, 0.35) if fit.energies else 1.0, percentile(fit.energies, 0.7) if fit.energies else 2.0, inverse=True), "Average bootstrap uncertainty relative to the fitted energy scale."),
        VisualDiagnostic("barrier_p50", percentile(barriers, 0.5) if barriers else 0.0, "good" if barriers else "warning", "Median estimated barrier among candidate attractor transitions."),
        VisualDiagnostic("energy_dynamic_range", percentile(fit.energies, 0.95) - percentile(fit.energies, 0.05) if fit.energies else 0.0, "good" if fit.energies else "critical", "Robust dynamic range of the effective energy coordinate."),
    ]


def diagnostics_payload(fit: LandscapeFit) -> dict[str, Any]:
    diagnostics = diagnose_fit(fit)
    counts = {"good": 0, "warning": 0, "critical": 0}
    for diagnostic in diagnostics:
        counts[diagnostic.status] += 1
    return {"diagnostics": [item.to_dict() for item in diagnostics], "status_counts": counts, "overall": "critical" if counts["critical"] else "warning" if counts["warning"] else "good"}
