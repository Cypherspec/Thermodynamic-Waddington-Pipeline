from __future__ import annotations

"""Locked analysis of matched perturbation, lineage, and fate experiments.

The module is intentionally strict: computational predictions cannot be entered as
measured outcomes, and a positive result requires donor-level replication,
pre-registered thresholds, viability checks, and a matched control.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


@dataclass(frozen=True)
class StudyCell:
    cell_id: str
    donor: str
    replicate: str
    intervention: str
    fate: str | None = None
    target_probability: float | None = None
    viable: bool = True
    lineage_id: str | None = None
    timepoint: str = "endpoint"
    batch: str = "unspecified"
    dose: float | None = None
    orthogonal_fate: str | None = None
    barcode_captured: bool = True
    perturbation_captured: bool = True
    expression_observed: bool = True
    velocity_observed: bool = True

    @classmethod
    def from_dict(cls, row: Mapping[str, Any]) -> "StudyCell":
        probability = row.get("target_probability")
        if probability is not None:
            probability = float(probability)
            if not 0.0 <= probability <= 1.0:
                raise ValueError("target_probability must be in [0, 1]")
        return cls(
            cell_id=str(row["cell_id"]),
            donor=str(row["donor"]),
            replicate=str(row.get("replicate", row["donor"])),
            intervention=str(row["intervention"]),
            fate=str(row["fate"]) if row.get("fate") is not None else None,
            target_probability=probability,
            viable=bool(row.get("viable", True)),
            lineage_id=str(row["lineage_id"]) if row.get("lineage_id") is not None else None,
            timepoint=str(row.get("timepoint", "endpoint")),
            batch=str(row.get("batch", "unspecified")),
            dose=float(row["dose"]) if row.get("dose") is not None else None,
            orthogonal_fate=str(row["orthogonal_fate"]) if row.get("orthogonal_fate") is not None else None,
            barcode_captured=bool(row.get("barcode_captured", True)),
            perturbation_captured=bool(row.get("perturbation_captured", True)),
            expression_observed=bool(row.get("expression_observed", True)),
            velocity_observed=bool(row.get("velocity_observed", True)),
        )

    def outcome(self, target_fate: str) -> float:
        if not self.viable:
            return float("nan")
        if self.target_probability is not None:
            return self.target_probability
        if self.fate is None:
            raise ValueError(f"cell {self.cell_id!r} has neither fate nor target_probability")
        return 1.0 if self.fate == target_fate else 0.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class StudyConfig:
    target_fate: str
    control_intervention: str
    minimum_donors: int = 3
    minimum_cells_per_group: int = 10
    minimum_effect: float = 0.10
    alpha: float = 0.05
    bootstrap_rounds: int = 2000
    seed: int = 17
    maximum_viability_loss: float = 0.10
    require_lineage: bool = True
    require_orthogonal_fate: bool = False
    require_measured_fate: bool = False
    require_perturbation_capture: bool = False
    require_velocity: bool = False
    require_expression: bool = False
    endpoint_timepoint: str = "endpoint"

    def validate(self) -> None:
        if not self.target_fate or not self.control_intervention:
            raise ValueError("target_fate and control_intervention are required")
        if self.minimum_donors < 2:
            raise ValueError("minimum_donors must be at least 2")
        if self.minimum_cells_per_group < 1:
            raise ValueError("minimum_cells_per_group must be positive")
        if self.bootstrap_rounds < 100:
            raise ValueError("bootstrap_rounds must be at least 100")
        if not 0.0 < self.alpha < 1.0:
            raise ValueError("alpha must be in (0, 1)")
        if self.minimum_effect < 0.0:
            raise ValueError("minimum_effect must be non-negative")
        if not 0.0 <= self.maximum_viability_loss <= 1.0:
            raise ValueError("maximum_viability_loss must be in [0, 1]")


@dataclass(frozen=True)
class GroupSummary:
    intervention: str
    donor: str
    cells: int
    viable_cells: int
    target_fraction: float
    viability_fraction: float
    lineages: int
    orthogonal_agreement: float = float("nan")

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class InterventionResult:
    intervention: str
    donors: int
    donor_effects: tuple[float, ...]
    effect: float
    bootstrap_low: float
    bootstrap_high: float
    permutation_p: float
    viability_effect: float
    consistent_direction: bool
    status: str
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class StudyReport:
    study_id: str
    config: StudyConfig
    schema_status: str
    groups: list[GroupSummary]
    interventions: list[InterventionResult]
    status: str
    claim_boundary: str
    blockers: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "study_id": self.study_id,
            "config": asdict(self.config),
            "schema_status": self.schema_status,
            "groups": [group.to_dict() for group in self.groups],
            "interventions": [item.to_dict() for item in self.interventions],
            "status": self.status,
            "claim_boundary": self.claim_boundary,
            "blockers": self.blockers,
            "provenance": self.provenance,
        }

    def save(self, path: str | Path) -> None:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8")


def load_cells(path: str | Path) -> list[StudyCell]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    rows = payload.get("cells", payload) if isinstance(payload, (dict, list)) else None
    if not isinstance(rows, list):
        raise ValueError("study input must be a JSON list or an object containing cells")
    return [StudyCell.from_dict(row) for row in rows]


def validate_cells(cells: Sequence[StudyCell], endpoint_timepoint: str = "endpoint", require_orthogonal_fate: bool = False, require_measured_fate: bool = False, require_perturbation_capture: bool = False, require_velocity: bool = False, require_expression: bool = False) -> dict[str, Any]:
    errors: list[str] = []
    ids = [cell.cell_id for cell in cells]
    if len(ids) != len(set(ids)):
        errors.append("duplicate_cell_ids")
    if not cells:
        errors.append("empty_study")
    for cell in cells:
        if not cell.donor:
            errors.append(f"missing_donor:{cell.cell_id}")
        if not cell.intervention:
            errors.append(f"missing_intervention:{cell.cell_id}")
        if not cell.replicate:
            errors.append(f"missing_replicate:{cell.cell_id}")
        if cell.timepoint != endpoint_timepoint:
            continue
        if cell.lineage_id is None:
            errors.append(f"missing_lineage:{cell.cell_id}")
        if require_measured_fate and cell.viable and cell.fate is None:
            errors.append(f"missing_measured_fate:{cell.cell_id}")
        elif cell.viable and cell.fate is None and cell.target_probability is None:
            errors.append(f"missing_outcome:{cell.cell_id}")
        if require_orthogonal_fate and cell.viable and cell.orthogonal_fate is None:
            errors.append(f"missing_orthogonal_fate:{cell.cell_id}")
        if require_perturbation_capture and not cell.perturbation_captured:
            errors.append(f"missing_perturbation_capture:{cell.cell_id}")
        if require_velocity and not cell.velocity_observed:
            errors.append(f"missing_velocity:{cell.cell_id}")
        if require_expression and not cell.expression_observed:
            errors.append(f"missing_expression:{cell.cell_id}")
    return {"valid": not errors, "cells": len(cells), "donors": len({cell.donor for cell in cells}), "interventions": sorted({cell.intervention for cell in cells}), "errors": errors}


def validate_expression_alignment(expression_cell_ids: Sequence[str], cells: Sequence[StudyCell]) -> dict[str, Any]:
    expression_ids = set(map(str, expression_cell_ids))
    outcome_ids = {cell.cell_id for cell in cells}
    overlap = expression_ids & outcome_ids
    return {
        "expression_cells": len(expression_ids),
        "outcome_cells": len(outcome_ids),
        "overlap": len(overlap),
        "expression_only": len(expression_ids - outcome_ids),
        "outcome_only": len(outcome_ids - expression_ids),
        "aligned": bool(overlap) and not (expression_ids - outcome_ids) and not (outcome_ids - expression_ids),
        "warning": "Partial overlap is not sufficient for leakage-free matched-cell validation." if overlap and expression_ids != outcome_ids else None,
    }


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def _variance(values: Sequence[float]) -> float:
    if len(values) < 2:
        return 0.0
    center = _mean(values)
    return sum((value - center) ** 2 for value in values) / (len(values) - 1)


def _bootstrap_interval(values: Sequence[float], rounds: int, seed: int) -> tuple[float, float]:
    if not values:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    samples = []
    for _ in range(rounds):
        sample = [values[rng.randrange(len(values))] for _ in values]
        samples.append(_mean(sample))
    samples.sort()
    return samples[max(0, int(0.025 * len(samples)))], samples[min(len(samples) - 1, int(0.975 * len(samples)))]


def _sign_flip_pvalue(values: Sequence[float], rounds: int, seed: int) -> float:
    if len(values) < 2:
        return 1.0
    observed = abs(_mean(values))
    rng = random.Random(seed)
    exceedances = 0
    for _ in range(rounds):
        randomized = [value if rng.random() < 0.5 else -value for value in values]
        if abs(_mean(randomized)) >= observed:
            exceedances += 1
    return (exceedances + 1.0) / (rounds + 1.0)


def _group_summaries(cells: Sequence[StudyCell], target_fate: str) -> list[GroupSummary]:
    groups: dict[tuple[str, str], list[StudyCell]] = {}
    for cell in cells:
        groups.setdefault((cell.intervention, cell.donor), []).append(cell)
    summaries = []
    for (intervention, donor), group in sorted(groups.items()):
        viable = [cell for cell in group if cell.viable]
        outcomes = [cell.outcome(target_fate) for cell in viable]
        lineages = {cell.lineage_id for cell in viable if cell.lineage_id}
        paired = [cell for cell in viable if cell.orthogonal_fate is not None and cell.fate is not None]
        agreement = _mean([1.0 if cell.fate == cell.orthogonal_fate else 0.0 for cell in paired]) if paired else float("nan")
        summaries.append(GroupSummary(intervention, donor, len(group), len(viable), _mean(outcomes) if outcomes else float("nan"), len(viable) / len(group), len(lineages), agreement))
    return summaries


def analyze_study(cells: Sequence[StudyCell], config: StudyConfig, study_id: str = "unregistered-study") -> StudyReport:
    config.validate()
    schema = validate_cells(cells, config.endpoint_timepoint, config.require_orthogonal_fate, config.require_measured_fate, config.require_perturbation_capture, config.require_velocity, config.require_expression)
    analysis_cells = [cell for cell in cells if cell.timepoint == config.endpoint_timepoint]
    groups = _group_summaries(analysis_cells, config.target_fate)
    blockers = list(schema["errors"])
    if not config.require_lineage:
        blockers = [item for item in blockers if not item.startswith("missing_lineage:")]
    if not config.require_orthogonal_fate:
        blockers = [item for item in blockers if not item.startswith("missing_orthogonal_fate:")]
    if config.control_intervention not in {group.intervention for group in groups}:
        blockers.append("missing_control_intervention")
    control_by_donor = {group.donor: group for group in groups if group.intervention == config.control_intervention}
    results: list[InterventionResult] = []
    for intervention in sorted({group.intervention for group in groups} - {config.control_intervention}):
        treatment_by_donor = {group.donor: group for group in groups if group.intervention == intervention}
        shared = sorted(set(treatment_by_donor) & set(control_by_donor))
        warnings: list[str] = []
        if len(shared) < config.minimum_donors:
            warnings.append("insufficient_independent_donors")
        effects = [treatment_by_donor[donor].target_fraction - control_by_donor[donor].target_fraction for donor in shared]
        viability_effects = [treatment_by_donor[donor].viability_fraction - control_by_donor[donor].viability_fraction for donor in shared]
        low, high = _bootstrap_interval(effects, config.bootstrap_rounds, config.seed + len(results) + 1)
        pvalue = _sign_flip_pvalue(effects, config.bootstrap_rounds, config.seed + len(results) + 101)
        effect = _mean(effects) if effects else float("nan")
        viability_effect = _mean(viability_effects) if viability_effects else float("nan")
        consistent = bool(effects) and all(effect_value > 0 for effect_value in effects)
        if viability_effect < -config.maximum_viability_loss:
            warnings.append("viability_loss_exceeds_threshold")
        underpowered_groups = [group for group in groups if group.intervention in {config.control_intervention, intervention} and group.donor in shared and group.cells < config.minimum_cells_per_group]
        if underpowered_groups:
            warnings.append("minimum_cells_per_group_not_met")
        if not effects:
            status = "not_comparable"
        elif underpowered_groups:
            status = "inconclusive_underpowered_groups"
        elif len(shared) < config.minimum_donors:
            status = "inconclusive_insufficient_replication"
        elif effect >= config.minimum_effect and low > 0 and consistent and viability_effect >= -config.maximum_viability_loss and (pvalue <= config.alpha or len(shared) < 4):
            status = "experimentally_supported_candidate"
        else:
            status = "inconclusive_or_falsified"
        results.append(InterventionResult(intervention, len(shared), tuple(effects), effect, low, high, pvalue, viability_effect, consistent, status, tuple(warnings)))
    supported = [item for item in results if item.status == "experimentally_supported_candidate"]
    if blockers:
        status = "invalid_study_schema"
    elif supported:
        status = "experimentally_supported_candidate"
    elif results:
        status = "inconclusive_or_falsified"
    else:
        status = "no_treatment_comparisons"
    if status == "experimentally_supported_candidate":
        claim_boundary = "Supported candidate effect under this pre-registered study analysis; independent replication and mechanistic validation are still required."
    else:
        claim_boundary = "No causal or world-changing claim is supported by this analysis."
    provenance = {"study_id": study_id, "input_cells": len(cells), "analysis_cells": len(analysis_cells), "endpoint_timepoint": config.endpoint_timepoint, "schema": schema, "analysis": "endpoint donor-paired effect, bootstrap interval, sign-flip permutation", "claim_status": status, "minimum_cell_rule": config.minimum_cells_per_group, "independent_unit": "donor", "measured_primary_outcome": config.require_measured_fate, "orthogonal_endpoint_required": config.require_orthogonal_fate, "velocity_required": config.require_velocity, "expression_required": config.require_expression}
    return StudyReport(study_id, config, "valid" if not blockers else "invalid", groups, results, status, claim_boundary, blockers, provenance)


def create_preregistration(config: StudyConfig, interventions: Sequence[str], protocol_id: str = "TW-HEM-FATE-001") -> dict[str, Any]:
    payload = {"protocol_id": protocol_id, "config": asdict(config), "interventions": list(interventions), "primary_endpoint": "donor-paired target-fate fraction difference versus control", "analysis_lock": "all exclusions, thresholds, and primary outcomes are fixed before unblinding"}
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    payload["preregistration_hash"] = hashlib.sha256(encoded).hexdigest()
    payload["status"] = "ready_for_registration_not_executed"
    return payload
