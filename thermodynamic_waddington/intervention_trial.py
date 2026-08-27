from __future__ import annotations

"""Pre-registered analysis of measured intervention studies.

This module is deliberately conservative: it accepts measured outcomes but
never turns observational or synthetic rows into causal evidence. Donor-level
aggregation is the inferential unit; cells are nested observations.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class TrialObservation:
    cell_id: str
    donor_id: str
    intervention: str
    fate_probability: float
    viability: float
    lineage_id: str | None = None
    barcode_id: str | None = None
    assay_batch: str | None = None
    orthogonal_fate: str | None = None
    endpoint: str = "endpoint"


@dataclass(frozen=True)
class TrialConfig:
    control_intervention: str = "control"
    treatment_intervention: str = "treatment"
    target_fate: str = "target"
    minimum_donors: int = 3
    minimum_cells_per_donor: int = 8
    minimum_viability: float = 0.7
    bootstrap_rounds: int = 2000
    seed: int = 41
    preregistration_id: str = "TW-INTERVENTION-001"


@dataclass
class TrialReport:
    status: str
    donor_effects: list[dict[str, Any]] = field(default_factory=list)
    aggregate: dict[str, Any] = field(default_factory=dict)
    exclusions: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    claim_boundary: str = "Measured intervention analysis is not a discovery unless preregistered endpoints, independent donors, matched controls, lineage readout, viability, and orthogonal fate validation pass."
    preregistration_hash: str = ""

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
        return payload


def _finite(value: float) -> bool:
    return math.isfinite(float(value))


def _mean(values: Sequence[float]) -> float | None:
    return sum(values) / len(values) if values else None


def _interval(values: Sequence[float], rounds: int, seed: int) -> dict[str, Any]:
    if not values:
        return {"estimate": None, "lower": None, "upper": None, "rounds": 0}
    rng = random.Random(seed)
    draws = [sum(rng.choice(values) for _ in values) / len(values) for _ in range(max(1, rounds))]
    draws.sort()
    return {"estimate": _mean(values), "lower": draws[int(0.025 * (len(draws) - 1))], "upper": draws[int(0.975 * (len(draws) - 1))], "rounds": len(draws)}


def _as_observation(row: Mapping[str, Any], index: int) -> TrialObservation:
    return TrialObservation(
        cell_id=str(row.get("cell_id", f"cell-{index}")),
        donor_id=str(row.get("donor_id", row.get("donor", "missing_donor"))),
        intervention=str(row.get("intervention", "missing_intervention")),
        fate_probability=float(row.get("fate_probability", row.get("prediction", row.get("target_probability", float("nan"))))),
        viability=float(row.get("viability", 1.0)),
        lineage_id=None if row.get("lineage_id") is None else str(row["lineage_id"]),
        barcode_id=None if row.get("barcode_id") is None else str(row["barcode_id"]),
        assay_batch=None if row.get("assay_batch") is None else str(row["assay_batch"]),
        orthogonal_fate=None if row.get("orthogonal_fate") is None else str(row["orthogonal_fate"]),
        endpoint=str(row.get("endpoint", "endpoint")),
    )


def load_observations(path: str | Path) -> list[TrialObservation]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    rows = payload.get("observations", payload.get("cells", payload)) if isinstance(payload, (dict, list)) else []
    if not isinstance(rows, list):
        raise ValueError("trial input must be a list or contain observations/cells")
    return [_as_observation(row, i) for i, row in enumerate(rows)]


def preregistration(config: TrialConfig) -> dict[str, Any]:
    endpoints = {
        "primary": "donor-level treatment minus matched-control change in target-fate probability",
        "secondary": ["viability", "lineage/barcode concordance", "orthogonal endpoint concordance"],
        "exclusions": ["nonviable cells", "missing donor identity", "missing intervention", "failed assay QC"],
        "inferential_unit": "independent donor, not cell",
    }
    payload = {"protocol_id": config.preregistration_id, "config": asdict(config), "endpoints": endpoints, "status": "preregistered_analysis_specification"}
    payload["preregistration_hash"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return payload


def write_preregistration(output: str | Path, config: TrialConfig = TrialConfig()) -> dict[str, Any]:
    result = preregistration(config)
    target = Path(output); target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def analyze_trial(observations: Sequence[TrialObservation], config: TrialConfig = TrialConfig(), preregistration_hash: str | None = None) -> TrialReport:
    blockers: list[str] = []
    exclusions: list[str] = []
    valid = [row for row in observations if row.donor_id != "missing_donor" and row.intervention != "missing_intervention" and _finite(row.fate_probability) and _finite(row.viability)]
    exclusions.append(f"excluded {len(observations) - len(valid)} rows with missing identity, intervention, or non-finite measurement")
    if not preregistration_hash:
        blockers.append("preregistration hash is missing")
    if len({row.donor_id for row in valid}) < config.minimum_donors:
        blockers.append(f"fewer than {config.minimum_donors} independent donors")
    control = [row for row in valid if row.intervention == config.control_intervention]
    treatment = [row for row in valid if row.intervention == config.treatment_intervention]
    if not control or not treatment:
        blockers.append("matched control and treatment observations are required")
    donor_effects: list[dict[str, Any]] = []
    for donor in sorted({row.donor_id for row in valid}):
        c = [row for row in control if row.donor_id == donor and row.viability >= config.minimum_viability]
        t = [row for row in treatment if row.donor_id == donor and row.viability >= config.minimum_viability]
        if len(c) < config.minimum_cells_per_donor or len(t) < config.minimum_cells_per_donor:
            continue
        cmean = _mean([row.fate_probability for row in c]); tmean = _mean([row.fate_probability for row in t])
        lineage = sum(row.lineage_id is not None or row.barcode_id is not None for row in c + t) / len(c + t)
        orthogonal = sum(row.orthogonal_fate is not None for row in c + t) / len(c + t)
        donor_effects.append({"donor_id": donor, "control_cells": len(c), "treatment_cells": len(t), "control_mean": cmean, "treatment_mean": tmean, "effect": tmean - cmean, "lineage_or_barcode_fraction": lineage, "orthogonal_fate_fraction": orthogonal, "control_viability": _mean([r.viability for r in c]), "treatment_viability": _mean([r.viability for r in t])})
    if not donor_effects:
        blockers.append("no donor has the preregistered minimum cells in both arms")
    effects = [float(row["effect"]) for row in donor_effects]
    if donor_effects and min(float(row["lineage_or_barcode_fraction"]) for row in donor_effects) < 1.0:
        blockers.append("lineage or barcode measurements are incomplete")
    if donor_effects and min(float(row["orthogonal_fate_fraction"]) for row in donor_effects) < 1.0:
        blockers.append("orthogonal fate measurements are incomplete")
    aggregate = {"donor_count": len(donor_effects), "effect": _interval(effects, config.bootstrap_rounds, config.seed), "viability": {"control": _mean([float(row["control_viability"]) for row in donor_effects]), "treatment": _mean([float(row["treatment_viability"]) for row in donor_effects])}, "inferential_unit": "donor", "cell_pooling_used": False}
    status = "measured_intervention_result_ready_for_replication" if not blockers else "blocked_pending_measured_experiment"
    return TrialReport(status, donor_effects, aggregate, exclusions, blockers, preregistration_hash=preregistration_hash or "")


def write_trial_report(observations: Sequence[TrialObservation], output: str | Path, config: TrialConfig = TrialConfig(), preregistration_hash: str | None = None) -> dict[str, Any]:
    report = analyze_trial(observations, config, preregistration_hash)
    target = Path(output); target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report.to_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return report.to_dict()
