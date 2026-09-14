"""High-level one-call API.

`analyze` runs the pipeline end to end and returns a compact, serializable
report: the entropy-production test for irreversibility, the kT landscape depth,
attractors, and - when progenitor and terminal cell types are given - the
committor commitment coordinate and the commitment free-energy barrier. It is the
recommended entry point; the individual modules remain available for detail.

    from thermodynamic_waddington import analyze
    report = analyze(expression, velocity, labels=labels,
                     source_labels=["Ductal"], target_labels=["Beta"])
    print(report.is_irreversible, report.entropy_production_pvalue)
    print(report.commitment_label, report.commitment_barrier_kt)
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Sequence

from .calibration import calibrate_fit
from .config import FitConfig
from .developmental import (
    commitment_profile,
    committor_free_energy_profile,
    developmental_coordinate,
)
from .model import fit_landscape


@dataclass
class AnalysisReport:
    n_cells: int
    n_edges: int
    n_attractors: int
    entropy_production_rate: float | None
    entropy_production_pvalue: float | None
    is_irreversible: bool | None
    landscape_range_kt: float | None
    path_vs_density_r2: float | None
    committor_order: list[str] | None = None
    commitment_label: str | None = None
    commitment_barrier_kt: float | None = None
    transition_state_q: float | None = None
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2))


def analyze(
    expression: Sequence[Sequence[float]],
    velocity: Sequence[Sequence[float]],
    labels: Sequence[str] | None = None,
    source_labels: Sequence[str] | None = None,
    target_labels: Sequence[str] | None = None,
    config: FitConfig | None = None,
    significance: float = 0.05,
) -> AnalysisReport:
    """Fit the landscape and return the headline thermodynamic results.

    Pass `labels` plus `source_labels` (progenitor) and `target_labels`
    (terminal) to also get the committor commitment coordinate and barrier.
    """
    fit = fit_landscape(expression, velocity, config=config, labels=labels)
    ep = fit.diagnostics.get("entropy_production", {}) or {}
    ep_rate = ep.get("entropy_production_rate")
    ep_p = ep.get("permutation_p_value")

    try:
        cal = calibrate_fit(fit)
        range_kt = cal.boltzmann_energy_range_kt
        r2 = cal.r_squared
    except Exception as exc:  # calibration needs a spread of embedding points
        range_kt = None
        r2 = None
        cal = None

    report = AnalysisReport(
        n_cells=len(fit.energies),
        n_edges=len(fit.edges),
        n_attractors=len(fit.attractors),
        entropy_production_rate=ep_rate,
        entropy_production_pvalue=ep_p,
        is_irreversible=(ep_p <= significance) if isinstance(ep_p, (int, float)) else None,
        landscape_range_kt=range_kt,
        path_vs_density_r2=r2,
    )
    if cal is None:
        report.warnings.append("calibration unavailable for this embedding")

    if labels is not None and source_labels and target_labels:
        try:
            q = developmental_coordinate(fit, source_labels, target_labels)
            prof = commitment_profile(fit, source_labels, target_labels)
            barrier = committor_free_energy_profile(q, temperature=(config or FitConfig()).temperature)
            report.committor_order = prof.order
            report.commitment_label = prof.commitment_label
            report.commitment_barrier_kt = barrier.barrier_kt
            report.transition_state_q = barrier.barrier_q
        except ValueError as exc:
            report.warnings.append(f"commitment coordinate skipped: {exc}")
    elif source_labels or target_labels:
        report.warnings.append("source_labels and target_labels both required for the committor")

    return report
