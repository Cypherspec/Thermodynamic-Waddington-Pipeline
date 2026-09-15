"""High-level one-call API.

`analyze` runs the pipeline end to end on plain matrices and returns a compact,
serializable report. `analyze_adata` does the same directly on an AnnData object,
reading the velocity and expression layers and writing the per-cell committor and
energy back into `adata.obs`, so the pipeline drops into a scanpy/scVelo workflow.

    from thermodynamic_waddington import analyze, analyze_adata

    report = analyze(expression, velocity, labels=labels,
                     source_labels=["Ductal"], target_labels=["Beta"])

    report = analyze_adata(adata, label_key="clusters",
                           source=["Ductal"], target=["Beta"])
    adata.obs["tw_committor"]   # written back
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


def _report_from_fit(fit, config, labels, source_labels, target_labels, significance):
    """Build the report and the per-cell committor from a fitted landscape."""
    ep = fit.diagnostics.get("entropy_production", {}) or {}
    ep_rate = ep.get("entropy_production_rate")
    ep_p = ep.get("permutation_p_value")

    try:
        cal = calibrate_fit(fit)
        range_kt = cal.boltzmann_energy_range_kt
        r2 = cal.r_squared
    except Exception:  # calibration needs a spread of embedding points
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

    committor = None
    if labels is not None and source_labels and target_labels:
        try:
            committor = developmental_coordinate(fit, source_labels, target_labels)
            prof = commitment_profile(fit, source_labels, target_labels)
            barrier = committor_free_energy_profile(committor, temperature=(config or FitConfig()).temperature)
            report.committor_order = prof.order
            report.commitment_label = prof.commitment_label
            report.commitment_barrier_kt = barrier.barrier_kt
            report.transition_state_q = barrier.barrier_q
        except ValueError as exc:
            report.warnings.append(f"commitment coordinate skipped: {exc}")
    elif source_labels or target_labels:
        report.warnings.append("source_labels and target_labels both required for the committor")

    return report, committor


def analyze(
    expression: Sequence[Sequence[float]],
    velocity: Sequence[Sequence[float]],
    labels: Sequence[str] | None = None,
    source_labels: Sequence[str] | None = None,
    target_labels: Sequence[str] | None = None,
    config: FitConfig | None = None,
    significance: float = 0.05,
) -> AnalysisReport:
    """Fit the landscape and return the headline thermodynamic results."""
    fit = fit_landscape(expression, velocity, config=config, labels=labels)
    report, _ = _report_from_fit(fit, config, labels, source_labels, target_labels, significance)
    return report


def _dense(x):
    import numpy as np
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, dtype=float)


def analyze_adata(
    adata,
    velocity_layer: str = "velocity",
    expression_layer: str | None = None,
    label_key: str | None = None,
    source: Sequence[str] | None = None,
    target: Sequence[str] | None = None,
    config: FitConfig | None = None,
    n_top_genes: int = 50,
    significance: float = 0.05,
    key_added: str = "tw",
) -> AnalysisReport:
    """Run the pipeline on an AnnData and write results back into it.

    Expression is taken from `expression_layer` (or `Ms`, or `.X`). Velocity is
    taken from `velocity_layer`, or derived as an unspliced-minus-spliced proxy
    if a `velocity` layer is absent but `spliced`/`unspliced` are present. The
    top `n_top_genes` by variance are used. The per-cell committor and free
    energy are written to `adata.obs['{key_added}_committor']` and
    `adata.obs['{key_added}_energy']`, and the report to
    `adata.uns['{key_added}']`.
    """
    import numpy as np

    layers = getattr(adata, "layers", {})

    if expression_layer and expression_layer in layers:
        expr = _dense(layers[expression_layer])
    elif "Ms" in layers:
        expr = _dense(layers["Ms"])
    else:
        expr = _dense(adata.X)

    if velocity_layer in layers:
        vel = np.nan_to_num(_dense(layers[velocity_layer]))
    elif "spliced" in layers and "unspliced" in layers:
        vel = _dense(layers["unspliced"]) - _dense(layers["spliced"])
    else:
        raise ValueError(
            f"no '{velocity_layer}' layer and no spliced/unspliced to derive a proxy; "
            f"available layers: {list(layers)}"
        )

    top = np.argsort(expr.var(axis=0))[::-1][:n_top_genes]
    expr = expr[:, top]
    vel = vel[:, top]

    labels = None
    if label_key is not None and label_key in adata.obs:
        labels = [str(v) for v in adata.obs[label_key].to_numpy()]

    fit = fit_landscape(expr.tolist(), vel.tolist(), config=config, labels=labels)
    report, committor = _report_from_fit(fit, config, labels, source, target, significance)

    try:
        adata.obs[f"{key_added}_energy"] = np.asarray(fit.energies, dtype=float)
        if committor is not None:
            adata.obs[f"{key_added}_committor"] = np.asarray(committor, dtype=float)
        adata.uns[key_added] = report.to_dict()
    except Exception as exc:  # writing back is best-effort
        report.warnings.append(f"could not write results into adata: {exc}")

    return report
