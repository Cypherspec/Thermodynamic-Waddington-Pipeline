from __future__ import annotations

import math
from typing import Any, Sequence

from .arrays import mean, variance
from .graph import Edge, NeighborGraph, velocity_alignment


def _finite_fraction(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    return sum(1 for value in values if math.isfinite(float(value))) / len(values)


def scientific_audit(
    points: Sequence[Sequence[float]],
    velocities: Sequence[Sequence[float]],
    graph: NeighborGraph,
    edges: Sequence[Edge],
    diffusions: Sequence[float],
    diffusion_floor: float,
    temperature: float,
    metadata: dict[str, object] | None = None,
) -> dict[str, Any]:
    metadata = metadata or {}
    n_cells = len(points)
    all_components = [float(value) for row in velocities for value in row]
    speeds = [math.sqrt(sum(float(value) ** 2 for value in row)) for row in velocities]
    alignments = [edge.alignment for edge in edges]
    floor_fraction = sum(abs(value - diffusion_floor) <= max(1e-12, diffusion_floor * 1e-8) for value in diffusions) / max(1, len(diffusions))
    velocity_status = str(metadata.get("velocity_status", "observed"))
    observed = velocity_status not in {"not_observed", "missing", "placeholder", "imputed"}
    blockers: list[dict[str, str]] = []
    if not observed:
        blockers.append({"code": "velocity_not_observed", "severity": "critical", "message": "RNA velocity is absent, placeholder, imputed, or explicitly marked unobserved."})
    if floor_fraction > 0.25:
        blockers.append({"code": "diffusion_floor_saturation", "severity": "high", "message": "More than one quarter of cells use the diffusion floor; local noise is not identified in those cells."})
    if not edges:
        blockers.append({"code": "no_directed_edges", "severity": "critical", "message": "No velocity-aligned graph edges survived filtering."})
    elif mean(alignments) < 0.2:
        blockers.append({"code": "weak_directional_alignment", "severity": "high", "message": "The retained edge field is weakly aligned with the supplied velocity field."})
    if n_cells < 100:
        blockers.append({"code": "small_snapshot", "severity": "medium", "message": "The snapshot is small for estimating a cell-state landscape and bootstrap uncertainty."})
    if len(all_components) and variance(all_components) == 0:
        blockers.append({"code": "zero_velocity_variance", "severity": "high", "message": "Velocity components have zero variance; drift and diffusion cannot be separated."})
    if len(graph.edges) and len(edges) / len(graph.edges) < 0.15:
        blockers.append({"code": "sparse_directed_graph", "severity": "high", "message": "Few geometric neighbors survive directionality filtering, making path propagation fragile."})
    if not metadata.get("replicate_trajectories"):
        blockers.append({"code": "no_replicate_trajectories", "severity": "critical", "message": "No replicate time-resolved trajectories or explicit nonequilibrium work protocol were supplied."})
    if not metadata.get("calibrated_diffusion"):
        blockers.append({"code": "uncalibrated_diffusion", "severity": "critical", "message": "The diffusion scale is a local velocity-spread proxy, not a calibrated transcriptional noise coefficient."})
    if not metadata.get("held_out_lineage_validation"):
        blockers.append({"code": "no_held_out_validation", "severity": "high", "message": "No held-out lineage or perturbation validation was supplied."})
    protocol_declared = bool(metadata.get("protocol") or metadata.get("work_protocol") or metadata.get("protocol_defined"))
    path_work = metadata.get("path_work") if isinstance(metadata.get("path_work"), dict) else {"count": len(edges), "effective_sample_size": 0.0}
    path_protocol = metadata.get("path_protocol") if isinstance(metadata.get("path_protocol"), dict) else {"protocol_defined": protocol_declared}
    diffusion_tensor = metadata.get("diffusion_tensor") if isinstance(metadata.get("diffusion_tensor"), dict) else {"cells": n_cells}
    current = metadata.get("current") if isinstance(metadata.get("current"), dict) else {"cells": n_cells}
    if not protocol_declared:
        blockers.append({"code": "protocol_defined_work_ensemble_missing", "severity": "critical", "message": "The observational graph is not a protocol-defined work ensemble."})
    required = [
        "controlled repeated trajectories or an explicitly justified path measure",
        "calibrated diffusion/noise scale",
        "independent trajectories and reverse/control checks",
        "held-out lineage and perturbation validation",
    ]
    return {
        "version": "0.3",
        "status": "provisional" if n_cells < 30 else "audited_effective_landscape",
        "protocol_declared": protocol_declared,
        "velocity_observed": observed,
        "velocity_status": velocity_status,
        "path_work": path_work,
        "path_protocol": path_protocol,
        "diffusion_tensor": diffusion_tensor,
        "current": current,
        "warnings": [blocker["message"] for blocker in blockers],
        "required_for_physical_claim": required,
        "claims": {
            "effective_landscape": True,
            "velocity_informed_if_observed": observed,
            "physical_free_energy_in_kT": False,
            "equilibrium_state_function_identified": False,
            "jarzynski_equality_validated": False,
        },
        "input_summary": {
            "cells": n_cells,
            "dimensions": len(points[0]) if points else 0,
            "velocity_components": len(velocities[0]) if velocities else 0,
            "finite_velocity_fraction": _finite_fraction(all_components),
            "mean_speed": mean(speeds) if speeds else 0.0,
            "directed_edges": len(edges),
            "geometric_edges": len(graph.edges),
            "mean_alignment": mean(alignments) if alignments else 0.0,
            "diffusion_floor_fraction": floor_fraction,
            "temperature_convention": temperature,
        },
        "blockers": blockers + (["protocol-defined work ensemble"] if not protocol_declared else []),
        "required_controls": [
            "replicate time-resolved trajectories or an explicitly defined intervention protocol",
            "calibrated diffusion/noise from replicate measurements or an instrument model",
            "forward and reverse protocol checks when invoking Crooks/Jarzynski estimators",
            "held-out lineage outcomes and perturbation response validation",
            "comparison against density-only, scVelo, CellRank, and Palantir baselines",
        ],
        "interpretation": "Use this output to audit an effective stochastic landscape. Do not report it as a molecular free-energy measurement without clearing the blockers and preregistering the protocol.",
    }
