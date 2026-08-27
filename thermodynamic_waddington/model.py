from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from .arrays import percentile, shape
from .config import FitConfig
from .graph import Edge, annotate_edges, build_knn, estimate_local_diffusion, local_density
from .jarzynski import barrier_heights, bootstrap_uncertainty, propagate_free_energy, attractor_candidates, jarzynski_diagnostics
from .linalg import pca, two_dimensional_projection
from .validation import fate_calibration, reversibility_gap, velocity_alignment_score
from .observability import fingerprint, make_provenance
from .geometry import summarize
from .topology.persistence import lower_star_pairs
from .topology.landscape import watershed_basins
from .analysis import basin_statistics, flux_matrix, report as analyze_report
from .dynamics import summarize_risk, transition_risk
from .metrics import basin_contrast, irreversibility_index, landscape_roughness, mean_alignment
from .causal import intervention_scan, rank_targets
from .spectral import spectral_summary
from .trajectory import simulate_trajectories, trajectory_summary
from .experiment import analyze_experiment
from .forecast import forecast
from .information import information_summary
from .lineage import lineage_summary
from .thermodynamics import effective_thermodynamic_summary
from .topology.cycles import cycle_summary
from .expansion import edge_work_distribution, current_balance, committor_iteration, reactive_flux, minimum_action_path, attractor_transition_matrix, graph_signature, summarize_frontier, adaptive_temperature, local_curvature, sensitivity_scan
from .cell_prediction import predict_population, prediction_summary
from .phenotype import phenotype_population, phenotype_summary
from .cell_atlas import interpret_population
from .multimodal import normalize_modalities, modality_summary
from .exploration import cell_measurement_ledger, population_context, tour_scenes
from .scientific_audit import scientific_audit
from .scientific import edge_path_samples, path_sample_summary, estimate_diffusion_tensor, tensor_summary, velocity_residuals, local_current, current_summary, bootstrap_path_work, protocol_sanity_checks
from .entropy_production import EntropyProductionConfig, estimate_entropy_production, cross_validate_with_jarzynski


def _to_2d_vectors(vectors: Sequence[Sequence[float]], n: int) -> list[list[float]]:
    if not vectors:
        return [[0.0, 0.0] for _ in range(n)]
    if len(vectors[0]) == 2:
        return [list(row) for row in vectors]
    return [list(row[:2]) + [0.0] * max(0, 2 - len(row[:2])) for row in vectors]


@dataclass
class LandscapeFit:
    config: dict[str, object]
    embedding: list[list[float]]
    energies: list[float]
    uncertainties: list[float]
    diffusions: list[float]
    densities: list[float]
    labels: list[str]
    cell_ids: list[str]
    edges: list[dict[str, float | int | str]]
    attractors: list[int]
    barriers: list[dict[str, float | int]]
    diagnostics: dict[str, object]
    metadata: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), separators=(",", ":")))

    @classmethod
    def load(cls, path: str | Path) -> "LandscapeFit":
        return cls(**json.loads(Path(path).read_text()))


def _serialize_edges(edges: Sequence[Edge]) -> list[dict[str, float | int | str]]:
    return [{"source": edge.source, "target": edge.target, "distance": edge.distance, "alignment": edge.alignment, "work": edge.work, "action": edge.action} for edge in edges]


def fit_landscape(expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], config: FitConfig | None = None, embedding: Sequence[Sequence[float]] | None = None, labels: Sequence[str] | None = None, cell_ids: Sequence[str] | None = None, lineage_outcomes: Sequence[str] | None = None, metadata: dict[str, object] | None = None) -> LandscapeFit:
    config = config or FitConfig()
    if not isinstance(config, FitConfig):
        raise TypeError("config must be a FitConfig instance")
    n_cells, n_features = shape(expression)
    velocity_rows, velocity_features = shape(velocity)
    if velocity_rows != n_cells:
        raise ValueError("expression and velocity must have the same number of cells")
    if velocity_features < 2:
        raise ValueError("velocity needs at least two dimensions")
    source_metadata = dict(metadata or {})
    measured_modalities = normalize_modalities(source_metadata, n_cells)
    source_metadata["modalities"] = measured_modalities
    source_metadata["modality_summary"] = modality_summary(measured_modalities, n_cells)
    velocity_status = str(source_metadata.get("velocity_status", "observed"))
    velocity_observed = velocity_status not in {"not_observed", "missing", "placeholder", "imputed"}
    config.validate(n_cells, n_features)
    reduced, eigenvalues, eigenvectors = pca(expression, config.dimensions, config.seed)
    points = [list(row) for row in reduced]
    if embedding is None:
        display_embedding = two_dimensional_projection(expression, config.seed)
    else:
        display_embedding = [list(map(float, row)) for row in embedding]
    graph = build_knn(points, config.neighbors)
    densities = local_density(points, graph, config.density_bandwidth)
    diffusions = estimate_local_diffusion(graph, velocity, config.diffusion_floor)
    residuals = velocity_residuals(points, velocity, graph)
    diffusion_tensors = estimate_diffusion_tensor(residuals, graph, config.diffusion_floor)
    edges = annotate_edges(points, velocity, graph, densities, diffusions, config.temperature, config.velocity_scale, config.edge_alignment_threshold)
    energies, counts, paths = propagate_free_energy(edges, n_cells, min(n_cells - 1, max(0, config.reference_index)), config.temperature, config.max_paths)
    uncertainties = bootstrap_uncertainty(edges, n_cells, min(n_cells - 1, max(0, config.reference_index)), config.temperature, config.max_paths, config.effective_bootstrap_replicates(), config.seed + 991)
    serialized = _serialize_edges(edges)
    labels = list(labels) if labels is not None else ["unlabeled"] * n_cells
    cell_ids = list(cell_ids) if cell_ids is not None else [f"cell_{i:05d}" for i in range(n_cells)]
    diagnostics: dict[str, object] = {
        "pca_explained_proxy": eigenvalues,
        "n_cells": n_cells,
        "n_features": n_features,
        "n_directed_edges": len(edges),
        "path_count_min": min(counts) if counts else 0,
        "path_count_max": max(counts) if counts else 0,
        "reference_cell": min(n_cells - 1, max(0, config.reference_index)),
        "lineage_outcomes_provided": lineage_outcomes is not None,
        "velocity_status": velocity_status,
        "velocity_observed": velocity_observed,
        "warnings": ["effective landscape; not an equilibrium thermodynamic state function"] + ([] if velocity_observed else ["RNA velocity was not observed; velocity alignment, directed flow, path work, transition pressure, and velocity-dependent free-energy estimates are unavailable or provisional."]),
    }
    attractors = attractor_candidates(edges, energies, config.barrier_quantile)
    barriers = barrier_heights(edges, energies)
    edge_pairs = [(edge.source, edge.target) for edge in edges]
    topology_pairs = lower_star_pairs(energies, edge_pairs)
    basin_records = watershed_basins(energies, edge_pairs)
    geometric = summarize(_to_2d_vectors(display_embedding, n_cells), _to_2d_vectors(velocity, n_cells), edges)
    path_samples = edge_path_samples(edges, config.temperature)
    currents = local_current(points, velocity, graph, densities, diffusions)
    path_diagnostics = jarzynski_diagnostics(paths, energies, config.temperature)
    entropy_production_report = None
    if config.enable_entropy_production and velocity_observed:
        entropy_production_report = estimate_entropy_production(
            points, velocity, densities, diffusions, graph,
            EntropyProductionConfig(
                temperature=config.temperature,
                velocity_scale=config.velocity_scale,
                bootstrap_replicates=config.entropy_production_bootstrap_replicates,
                permutation_replicates=config.entropy_production_permutation_replicates,
                seed=config.seed + 5231,
            ),
        )
    audit_metadata = dict(source_metadata)
    audit_metadata.update({"path_work": path_sample_summary(path_samples, config.temperature), "path_protocol": protocol_sanity_checks(path_samples), "diffusion_tensor": tensor_summary(diffusion_tensors), "current": current_summary(currents)})
    diagnostics["scientific_audit"] = scientific_audit(points, velocity, graph, edges, diffusions, config.diffusion_floor, config.temperature, audit_metadata)
    diagnostics["scientific_audit"]["velocity_observed"] = velocity_observed
    diagnostics["scientific_audit"]["velocity_status"] = velocity_status
    diagnostics["jarzynski_diagnostics"] = [record.to_dict() if hasattr(record, "to_dict") else dict(record) for record in path_diagnostics]
    diagnostics["jarzynski_path_diagnostics"] = diagnostics["jarzynski_diagnostics"]
    diagnostics["jarzynski_degenerate_fraction"] = sum(float(record.get("effective_sample_size", 0.0)) < 2.0 for record in diagnostics["jarzynski_diagnostics"]) / max(1, len(diagnostics["jarzynski_diagnostics"]))
    diagnostics.update({
        "velocity_alignment": velocity_alignment_score(serialized) if velocity_observed else None,
        "reversibility_gap": reversibility_gap(energies, serialized) if velocity_observed else None,
        "persistence_pairs": [pair.to_dict() for pair in topology_pairs],
        "basins": [basin.to_dict() for basin in basin_records],
        "geometric_summary": geometric.to_dict(),
        "coverage": sum(1 for value in energies if value == value and abs(value) < float("inf")) / max(1, n_cells),
        "lineage_fate_calibration": fate_calibration(energies, labels, lineage_outcomes, config.temperature) if lineage_outcomes is not None else None,
        "provenance": make_provenance({"expression": expression, "velocity": velocity, "metadata": metadata or {}}, ["pca", "knn", "diffusion", "path_work", "uncertainty", "topology"]).to_dict(),
        "input_fingerprint": fingerprint({"expression": expression, "velocity": velocity}),
        "advanced_metrics": {
            "mean_alignment": mean_alignment(edges),
            "irreversibility_index": irreversibility_index(edges),
            "landscape_roughness": landscape_roughness(energies),
            "basin_contrast": basin_contrast(energies, attractors),
        },
        "transition_risk": summarize_risk(transition_risk(edges, energies, config.temperature)),
        "trajectory_ensemble": trajectory_summary(simulate_trajectories(edges, energies, [min(attractors) if attractors else 0], steps=config.trajectory_steps, replicates=config.trajectory_replicates, temperature=config.temperature, seed=config.seed + 4001)),
        "basin_statistics": basin_statistics([basin.representative for basin in basin_records for _ in basin.members], energies, _to_2d_vectors(velocity, n_cells)),
        "flux_matrix": flux_matrix([next((basin.representative for basin in basin_records if index in basin.members), -1) for index in range(n_cells)], edges, energies, config.temperature),
        "counterfactual_targets": rank_targets(intervention_scan(edges, energies, damping=0.25), limit=12) if config.enable_counterfactuals else [],
        "spectral": spectral_summary(edges, energies) if config.enable_spectral_analysis else {},
        "experiment_summary": analyze_experiment(energies, serialized, labels, config.temperature),
        "forecast": [item.__dict__ for item in forecast([min(n_cells - 1, max(0, config.reference_index))], edges, energies, config.temperature)],
        "information_geometry": information_summary(points, diffusions, energies),
        "lineage_summary": lineage_summary(labels, lineage_outcomes),
        "thermodynamic_summary": effective_thermodynamic_summary(energies, diffusions, graph.edges, config.temperature),
        "cycle_summary": cycle_summary(LandscapeFit(asdict(config), display_embedding, energies, uncertainties, diffusions, densities, labels, cell_ids, serialized, attractors, barriers, {}, metadata or {})),
        "work_distribution": edge_work_distribution(edges, config.temperature),
        "current_balance": current_balance(edges, energies, config.temperature),
        "graph_signature": graph_signature(edges),
        "frontier_summary": summarize_frontier(energies, edges),
        "attractor_transition_matrix": attractor_transition_matrix(edges, energies, attractors, config.temperature),
        "adaptive_temperature": adaptive_temperature(energies),
        "local_curvature": local_curvature(energies, edges),
        "sensitivity_scan": sensitivity_scan(energies),
        "scientific_path_details": {
            "path_work": path_sample_summary(path_samples, config.temperature),
            "path_protocol": protocol_sanity_checks(path_samples),
            "path_work_bootstrap": bootstrap_path_work(edges, config.effective_bootstrap_replicates(), config.seed + 771),
            "diffusion_tensor": tensor_summary(diffusion_tensors),
            "current": current_summary(currents),
            "interpretation": "effective stochastic path-work analysis; physical kT claims require calibration, a specified protocol, and held-out experimental validation",
        },
        "entropy_production": entropy_production_report.to_dict() if entropy_production_report is not None else {
            "status": "unavailable",
            "reason": "No measured RNA velocity was supplied (or enable_entropy_production=False); entropy production requires a directional signal to distinguish from equilibrium.",
        },
        "entropy_production_cross_validation": (
            cross_validate_with_jarzynski(entropy_production_report, path_diagnostics)
            if entropy_production_report is not None else None
        ),
    })
    fit_metadata = dict(source_metadata)
    fit_metadata["modalities"] = measured_modalities
    fit_metadata["modality_summary"] = modality_summary(measured_modalities, n_cells)
    fit_metadata["velocity_status"] = velocity_status
    fit_metadata["velocity_observed"] = velocity_observed
    if isinstance(expression, list):
        fit_metadata.setdefault("expression", [list(map(float, row)) for row in expression])
    if not isinstance(fit_metadata.get("gene_names"), list) or len(fit_metadata.get("gene_names", [])) != n_features:
        fit_metadata["gene_names"] = [f"feature_{index}" for index in range(n_features)]
        fit_metadata.setdefault("gene_names_status", "synthetic_or_unannotated")
    else:
        fit_metadata.setdefault("gene_names_status", "provided_annotation")
    fit_metadata.setdefault("velocity_features", n_features)
    fit_metadata.setdefault("expression", [list(map(float, row)) for row in expression])
    fit_metadata.setdefault("velocity", [list(map(float, row)) for row in velocity])
    fit_metadata.setdefault("modalities", source_metadata.get("modalities", {}))
    fit_metadata.setdefault("modality_summary", source_metadata.get("modality_summary", {}))
    gene_names = fit_metadata.get("gene_names")
    if isinstance(expression, list) and isinstance(gene_names, list):
        predictions = predict_population(expression, gene_names, cell_ids, labels, energies, uncertainties)
        diagnostics["cell_predictions"] = [prediction.to_dict() for prediction in predictions]
        diagnostics["prediction_summary"] = prediction_summary(predictions)
        phenotypes = phenotype_population(expression, gene_names, cell_ids, labels, energies, uncertainties, lineage_outcomes)
        diagnostics["cell_phenotypes"] = [report.to_dict() for report in phenotypes]
        diagnostics["phenotype_summary"] = phenotype_summary(phenotypes)
        diagnostics["cell_interpretations"] = [report.to_dict() for report in interpret_population(expression, gene_names, cell_ids, labels, energies, uncertainties, lineage_outcomes)]
        diagnostics["cell_measurement_ledgers"] = [cell_measurement_ledger(measured_modalities, index) for index in range(n_cells)]
        diagnostics["cell_population_context"] = [population_context(measured_modalities, index) for index in range(n_cells)]
        live_cells = []
        for index in range(n_cells):
            live_cells.append({
                "index": index,
                "id": cell_ids[index],
                "label": labels[index],
                "energy": energies[index],
                "state_class": "deep basin / low effective energy" if energies[index] <= percentile(energies, 0.24) else "ridge / activated state" if energies[index] >= percentile(energies, 0.76) else "transition corridor / metastable state",
                "basin": next((basin_index for basin_index, basin in enumerate(basin_records) if index in basin.members), None),
                "velocity_observed": velocity_observed,
            })
        diagnostics["live_scenes"] = tour_scenes(live_cells, limit=min(10, max(3, n_cells)))
        diagnostics["exploration_summary"] = {
            "scene_count": len(diagnostics["live_scenes"]),
            "available_modalities": sorted(measured_modalities),
            "measurement_policy": "Measured matrices are carried through exactly; absent modalities remain unavailable.",
            "rna_landscape": "The Waddington field is fit from expression and measured RNA velocity when supplied.",
        }
    if not velocity_observed:
        diagnostics["velocity_dependent_fields"] = {"status": "unavailable", "reason": "No measured RNA velocity was supplied; zero velocity is not a measurement.", "fields": ["flow alignment", "directed transition pressure", "RNA-velocity path work", "velocity-informed free-energy interpretation"]}
        diagnostics["velocity_alignment"] = None
        diagnostics["reversibility_gap"] = None
        diagnostics["thermodynamic_summary"] = {"status": "density/expression-only provisional landscape", "warning": "Do not interpret F as a velocity-informed free-energy landscape without measured RNA velocity."}
    return LandscapeFit(asdict(config), display_embedding, energies, uncertainties, diffusions, densities, labels, cell_ids, serialized, attractors, barriers, diagnostics, fit_metadata)
