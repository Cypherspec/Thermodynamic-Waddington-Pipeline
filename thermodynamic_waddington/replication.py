from __future__ import annotations

"""Independent-study replication analysis.

This module makes the strongest useful computational upgrade without confusing
more cells with more biological experiments. The unit of replication is a
study or donor-level experimental run, never an individual cell. Reports retain
study fingerprints, exclusions, heterogeneity, and preregistered claim gates.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .flagship_benchmark import Observation, _metrics


@dataclass(frozen=True)
class ReplicationStudy:
    study_id: str
    provenance: str
    observations: tuple[Observation, ...]
    assay_version: str = "unspecified"
    preregistration_hash: str | None = None

    @property
    def fingerprint(self) -> str:
        payload = {
            "provenance": self.provenance,
            "assay_version": self.assay_version,
            "preregistration_hash": self.preregistration_hash,
            "observations": [asdict(row) for row in self.observations],
        }
        return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        return {
            "study_id": self.study_id,
            "provenance": self.provenance,
            "assay_version": self.assay_version,
            "preregistration_hash": self.preregistration_hash,
            "fingerprint": self.fingerprint,
            "observations": [asdict(row) for row in self.observations],
        }


@dataclass
class ReplicationReport:
    status: str
    independent_studies: int
    eligible_studies: int
    study_reports: list[dict[str, Any]] = field(default_factory=list)
    meta_analysis: dict[str, Any] = field(default_factory=dict)
    heterogeneity: dict[str, Any] = field(default_factory=dict)
    duplicate_fingerprints: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    claim_boundary: str = "Independent predictive replication is not causal evidence, molecular thermodynamics, or clinical validation."
    fingerprint: str = ""

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
        return payload


def _as_observation(row: Mapping[str, Any], prefix: str, index: int) -> Observation:
    return Observation(
        cell_id=str(row.get("cell_id", f"{prefix}-cell-{index}")),
        donor=str(row.get("donor", row.get("donor_id", "missing_donor"))),
        clone=str(row.get("clone", row.get("lineage_id", "missing_clone"))),
        timepoint=str(row.get("timepoint", "endpoint")),
        intervention=str(row.get("intervention", "unknown")),
        prediction=float(row.get("prediction", row.get("score", row.get("target_probability", float("nan"))))),
        outcome=float(row.get("outcome", row.get("target", row.get("fate", float("nan"))))),
        viable=bool(row.get("viable", True)),
    )


def load_studies(path: str | Path) -> list[ReplicationStudy]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    raw_studies = payload.get("studies", payload if isinstance(payload, list) else [])
    if not isinstance(raw_studies, list):
        raise ValueError("replication input must be a list or contain a studies list")
    studies: list[ReplicationStudy] = []
    for index, raw in enumerate(raw_studies):
        if not isinstance(raw, Mapping):
            raise ValueError(f"study {index} must be an object")
        study_id = str(raw.get("study_id", raw.get("id", f"study-{index}")))
        rows = raw.get("observations", raw.get("cells", []))
        if not isinstance(rows, list):
            raise ValueError(f"study {study_id} observations must be a list")
        observations = tuple(_as_observation(row, study_id, i) for i, row in enumerate(rows))
        studies.append(ReplicationStudy(study_id, str(raw.get("provenance", "unspecified")), observations, str(raw.get("assay_version", "unspecified")), raw.get("preregistration_hash")))
    return studies


def _bootstrap_mean(values: Sequence[float], rounds: int, seed: int) -> dict[str, float | int | None]:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if not finite:
        return {"estimate": None, "lower": None, "upper": None, "rounds": 0}
    rng = random.Random(seed)
    estimates = []
    for _ in range(max(1, int(rounds))):
        estimates.append(sum(rng.choice(finite) for _ in finite) / len(finite))
    estimates.sort()
    low = estimates[int(0.025 * (len(estimates) - 1))]
    high = estimates[int(0.975 * (len(estimates) - 1))]
    return {"estimate": sum(finite) / len(finite), "lower": low, "upper": high, "rounds": len(estimates)}


def _study_effect(study: ReplicationStudy) -> dict[str, Any]:
    eligible = [row for row in study.observations if row.viable and math.isfinite(row.prediction) and math.isfinite(row.outcome)]
    metrics = _metrics(eligible)
    return {
        "study_id": study.study_id,
        "provenance": study.provenance,
        "assay_version": study.assay_version,
        "fingerprint": study.fingerprint,
        "eligible_cells": len(eligible),
        "donors": len({row.donor for row in eligible if row.donor != "missing_donor"}),
        "clones": len({row.clone for row in eligible if row.clone != "missing_clone"}),
        "metrics": metrics,
    }


def evaluate_replication(studies: Sequence[ReplicationStudy], bootstrap_rounds: int = 256, min_cells: int = 8, min_studies: int = 2) -> ReplicationReport:
    fingerprints: dict[str, list[str]] = {}
    for study in studies:
        fingerprints.setdefault(study.fingerprint, []).append(study.study_id)
    duplicate_fingerprints = [fingerprint for fingerprint, ids in fingerprints.items() if len(ids) > 1]
    unique: list[ReplicationStudy] = []
    seen: set[str] = set()
    for study in studies:
        if study.fingerprint not in seen:
            unique.append(study)
            seen.add(study.fingerprint)
    reports = [_study_effect(study) for study in unique]
    eligible = [report for report in reports if int(report["eligible_cells"]) >= min_cells]
    effects = [float(report["metrics"]["spearman"]) for report in eligible if math.isfinite(float(report["metrics"]["spearman"]))]
    blockers: list[str] = []
    if duplicate_fingerprints:
        blockers.append("duplicate study fingerprints were deduplicated and cannot count as independent replication")
    if len(eligible) < min_studies:
        blockers.append(f"fewer than {min_studies} independent studies have at least {min_cells} eligible cells")
    if any(report["donors"] < 2 for report in eligible):
        blockers.append("each replication study should include at least two independent donors")
    if any(report["assay_version"] == "unspecified" for report in eligible):
        blockers.append("assay version or locked analysis provenance is missing for at least one study")
    mean_report = _bootstrap_mean(effects, bootstrap_rounds, 991)
    spread = None
    if len(effects) > 1:
        mean = sum(effects) / len(effects)
        spread = sum((effect - mean) ** 2 for effect in effects) / (len(effects) - 1)
    meta = {
        "unit_of_replication": "independent study",
        "effect": "held-out Spearman correlation per study",
        "study_level_effect": mean_report,
        "study_level_effects": effects,
        "cell_pooling_used": False,
    }
    heterogeneity = {
        "study_count": len(effects),
        "sample_variance": spread,
        "interpretation": "Heterogeneity is descriptive; it is not evidence of a shared causal mechanism.",
    }
    status = "replicated_predictive_result" if len(blockers) == 0 else "independent_replication_not_established"
    return ReplicationReport(status, len(unique), len(eligible), reports, meta, heterogeneity, duplicate_fingerprints, blockers)


def write_replication_report(studies: Sequence[ReplicationStudy], output: str | Path = "experiments/replication_report.json", bootstrap_rounds: int = 256) -> dict[str, Any]:
    report = evaluate_replication(studies, bootstrap_rounds=bootstrap_rounds)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report.to_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return report.to_dict()
