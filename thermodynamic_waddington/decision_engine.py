from __future__ import annotations

"""Evidence-gated decision support for cell-fate intervention programs.

This module connects landscape predictions, protocol diagnostics, robustness,
and measured-trial evidence without allowing any one computational score to
masquerade as causal biology. It produces auditable decisions: observe,
replicate, intervene, or stop. The decision is deliberately conservative when
required evidence is absent or contradictory.
"""

import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class DecisionConfig:
    target_fate: str = "target"
    minimum_effect: float = 0.10
    minimum_donor_count: int = 3
    minimum_viability: float = 0.70
    minimum_replication_studies: int = 2
    minimum_transportability_spearman: float = 0.30
    minimum_robustness_spearman: float = 0.70
    maximum_negative_control_abs_spearman: float = 0.20
    require_preregistration: bool = True
    require_lineage: bool = True
    require_orthogonal_fate: bool = True
    require_reverse_protocol: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class EvidenceRecord:
    name: str
    status: str
    value: float | None = None
    threshold: float | None = None
    provenance: str = ""
    detail: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class DecisionReport:
    decision: str
    attained_level: str
    confidence: str
    records: tuple[EvidenceRecord, ...]
    blockers: tuple[str, ...]
    warnings: tuple[str, ...]
    next_actions: tuple[str, ...]
    claim_boundary: str
    fingerprint: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "decision": self.decision,
            "attained_level": self.attained_level,
            "confidence": self.confidence,
            "records": [record.to_dict() for record in self.records],
            "blockers": list(self.blockers),
            "warnings": list(self.warnings),
            "next_actions": list(self.next_actions),
            "claim_boundary": self.claim_boundary,
            "fingerprint": self.fingerprint,
        }


def _as_mapping(value: Any) -> Mapping[str, Any]:
    if hasattr(value, "to_dict"):
        value = value.to_dict()
    return value if isinstance(value, Mapping) else {}


def _number(value: Any) -> float | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _record(name: str, status: str, value: Any = None, threshold: Any = None, provenance: str = "", detail: str = "") -> EvidenceRecord:
    return EvidenceRecord(name, status, _number(value), _number(threshold), provenance, detail)


def _find(mapping: Mapping[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in mapping:
            return mapping[key]
    return None


def _records_from_trial(trial: Mapping[str, Any], config: DecisionConfig) -> tuple[list[EvidenceRecord], list[str], list[str]]:
    records: list[EvidenceRecord] = []
    blockers: list[str] = []
    warnings: list[str] = []
    status = str(trial.get("status", "missing"))
    aggregate = _as_mapping(trial.get("aggregate"))
    effect = _as_mapping(aggregate.get("effect"))
    estimate = _number(effect.get("estimate"))
    donors = _number(aggregate.get("donor_count"))
    records.append(_record("measured_intervention_status", "pass" if status.startswith("measured_") else "fail", provenance="intervention_trial", detail=status))
    records.append(_record("donor_count", "pass" if donors is not None and donors >= config.minimum_donor_count else "fail", donors, config.minimum_donor_count, "intervention_trial"))
    records.append(_record("donor_level_effect", "pass" if estimate is not None and estimate >= config.minimum_effect else "fail", estimate, config.minimum_effect, "intervention_trial"))
    viability = _as_mapping(aggregate.get("viability"))
    for arm in ("control", "treatment"):
        value = _number(viability.get(arm))
        records.append(_record(f"{arm}_viability", "pass" if value is not None and value >= config.minimum_viability else "fail", value, config.minimum_viability, "intervention_trial"))
    for required, label in (("preregistration hash is missing", "preregistration"), ("lineage or barcode measurements are incomplete", "lineage_or_barcode"), ("orthogonal fate measurements are incomplete", "orthogonal_fate")):
        found = required in [str(item) for item in trial.get("blockers", [])]
        if label == "preregistration":
            found = not bool(trial.get("preregistration_hash"))
        required_enabled = (label == "preregistration" and config.require_preregistration) or (label == "lineage_or_barcode" and config.require_lineage) or (label == "orthogonal_fate" and config.require_orthogonal_fate)
        if required_enabled:
            records.append(_record(label, "fail" if found else "pass", provenance="intervention_trial"))
            if found:
                blockers.append(required)
    for blocker in trial.get("blockers", []):
        text = str(blocker)
        if text not in blockers and text != "preregistration hash is missing":
            blockers.append(text)
    for warning in trial.get("warnings", []):
        warnings.append(str(warning))
    return records, blockers, warnings


def _records_from_transportability(report: Mapping[str, Any], config: DecisionConfig) -> tuple[list[EvidenceRecord], list[str], list[str]]:
    records: list[EvidenceRecord] = []
    blockers: list[str] = []
    warnings: list[str] = []
    metrics = report.get("metrics", [])
    eligible = [item for item in metrics if isinstance(item, Mapping) and not bool(item.get("negative_control")) and _number(item.get("spearman")) is not None]
    count = len(eligible)
    pooled = _as_mapping(report.get("pooled_equal_study"))
    score = _number(pooled.get("spearman"))
    records.append(_record("independent_transportability_studies", "pass" if count >= config.minimum_replication_studies else "fail", count, config.minimum_replication_studies, "transportability"))
    records.append(_record("equal_study_spearman", "pass" if score is not None and score >= config.minimum_transportability_spearman else "fail", score, config.minimum_transportability_spearman, "transportability"))
    if count < config.minimum_replication_studies:
        blockers.append("independent replication studies are insufficient")
    if score is None or score < config.minimum_transportability_spearman:
        blockers.append("transportability predictive association is below threshold")
    negative = _as_mapping(report.get("negative_controls"))
    negative_score = _number(_find(negative, "max_abs_spearman", "maximum_abs_spearman"))
    if negative_score is not None:
        records.append(_record("negative_control_association", "pass" if negative_score <= config.maximum_negative_control_abs_spearman else "fail", negative_score, config.maximum_negative_control_abs_spearman, "transportability"))
        if negative_score > config.maximum_negative_control_abs_spearman:
            blockers.append("negative control association exceeds threshold")
    else:
        warnings.append("negative control score is not reported")
    for item in report.get("blockers", []):
        text = str(item)
        if text not in blockers:
            blockers.append(text)
    return records, blockers, warnings


def _records_from_robustness(report: Mapping[str, Any], config: DecisionConfig) -> tuple[list[EvidenceRecord], list[str], list[str]]:
    summary = _as_mapping(report.get("summary"))
    score = _number(summary.get("minimum_energy_spearman"))
    finite = _number(summary.get("finite_fraction"))
    records = [
        _record("robustness_finite_fraction", "pass" if finite is not None and finite == 1.0 else "fail", finite, 1.0, "robustness"),
        _record("minimum_energy_ordering_spearman", "pass" if score is not None and score >= config.minimum_robustness_spearman else "fail", score, config.minimum_robustness_spearman, "robustness"),
    ]
    blockers = [str(item) for item in report.get("blockers", [])]
    if score is None or score < config.minimum_robustness_spearman:
        blockers.append("landscape ordering is specification-sensitive")
    if finite is None or finite < 1.0:
        blockers.append("robustness scenarios contain failed or nonfinite fits")
    return records, sorted(set(blockers)), []


def _records_from_path_ensemble(report: Mapping[str, Any], config: DecisionConfig) -> tuple[list[EvidenceRecord], list[str], list[str]]:
    records: list[EvidenceRecord] = []
    blockers: list[str] = []
    warnings: list[str] = []
    protocols = report.get("protocols", [])
    records.append(_record("protocol_path_ensemble", "pass" if report.get("status") == "protocol_bound_ensemble_diagnostic" else "fail", provenance="path_ensemble"))
    if not protocols:
        blockers.append("no protocol-defined path ensemble")
    for protocol in protocols:
        if not isinstance(protocol, Mapping):
            continue
        forward = _number(protocol.get("forward_replicates")) or 0.0
        reverse = _number(protocol.get("reverse_replicates")) or 0.0
        records.append(_record(f"{protocol.get('protocol_id', 'protocol')}_forward_replicates", "pass" if forward >= 4 else "fail", forward, 4.0, "path_ensemble"))
        records.append(_record(f"{protocol.get('protocol_id', 'protocol')}_reverse_replicates", "pass" if reverse >= 4 else "fail", reverse, 4.0, "path_ensemble"))
    for item in report.get("blockers", []):
        blockers.append(str(item))
    warnings.extend(str(item) for item in report.get("warnings", []))
    if config.require_reverse_protocol and not any((_number(item.get("reverse_replicates")) or 0.0) >= 4 for item in protocols if isinstance(item, Mapping)):
        blockers.append("reverse protocol evidence is absent or underpowered")
    return records, sorted(set(blockers)), warnings


def decide(
    *,
    trial: Mapping[str, Any] | Any | None = None,
    transportability: Mapping[str, Any] | Any | None = None,
    robustness: Mapping[str, Any] | Any | None = None,
    path_ensemble: Mapping[str, Any] | Any | None = None,
    config: DecisionConfig | None = None,
) -> DecisionReport:
    config = config or DecisionConfig()
    records: list[EvidenceRecord] = []
    blockers: list[str] = []
    warnings: list[str] = []
    if trial is not None:
        new_records, new_blockers, new_warnings = _records_from_trial(_as_mapping(trial), config)
        records.extend(new_records); blockers.extend(new_blockers); warnings.extend(new_warnings)
    else:
        blockers.append("measured intervention trial report is missing")
        records.append(_record("measured_intervention_status", "fail", provenance="decision_engine"))
    if transportability is not None:
        new_records, new_blockers, new_warnings = _records_from_transportability(_as_mapping(transportability), config)
        records.extend(new_records); blockers.extend(new_blockers); warnings.extend(new_warnings)
    else:
        blockers.append("independent transportability report is missing")
    if robustness is not None:
        new_records, new_blockers, new_warnings = _records_from_robustness(_as_mapping(robustness), config)
        records.extend(new_records); blockers.extend(new_blockers); warnings.extend(new_warnings)
    else:
        blockers.append("robustness report is missing")
    if path_ensemble is not None:
        new_records, new_blockers, new_warnings = _records_from_path_ensemble(_as_mapping(path_ensemble), config)
        records.extend(new_records); blockers.extend(new_blockers); warnings.extend(new_warnings)
    else:
        warnings.append("protocol path ensemble was not supplied; computational landscape outputs remain effective diagnostics")
    blockers = sorted(set(blockers))
    warnings = sorted(set(warnings))
    passed = sum(record.status == "pass" for record in records)
    total = len(records)
    if not blockers:
        decision, level, confidence = "advance_to_replicated_intervention", "replicated_predictive_intervention_result", "high" if passed == total else "moderate"
        next_actions = ("pre-register confirmatory study", "replicate across independent donors and laboratories", "release raw data and locked analysis")
    elif trial is not None and not any("measured intervention" in blocker or "donor" in blocker for blocker in blockers):
        decision, level, confidence = "replicate_measured_signal", "measured_result_pending_replication", "moderate"
        next_actions = ("repeat the intervention with independent donors", "retain matched controls, lineage/barcode, viability, and orthogonal fate endpoints", "do not claim molecular free energy or causal generalization")
    else:
        decision, level, confidence = "do_not_upgrade_claim", "computational_or_incomplete_evidence", "low"
        next_actions = ("acquire measured randomized intervention data", "lock preregistration before unblinding", "report negative controls and failed endpoints")
    claim_boundary = "The decision engine ranks evidence states; it cannot establish a world-changing biological discovery. Causal and thermodynamic claims require measured, preregistered, independently replicated experiments."
    provisional = DecisionReport(decision, level, confidence, tuple(records), tuple(blockers), tuple(warnings), tuple(next_actions), claim_boundary, "")
    fingerprint = hashlib.sha256(json.dumps(provisional.to_dict(), sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
    return DecisionReport(decision, level, confidence, tuple(records), tuple(blockers), tuple(warnings), tuple(next_actions), claim_boundary, fingerprint)


def decide_from_files(*, trial_path: str | Path | None = None, transportability_path: str | Path | None = None, robustness_path: str | Path | None = None, path_ensemble_path: str | Path | None = None, config: DecisionConfig | None = None) -> dict[str, Any]:
    def load(path: str | Path | None) -> Any:
        return json.loads(Path(path).read_text(encoding="utf-8")) if path else None
    return decide(trial=load(trial_path), transportability=load(transportability_path), robustness=load(robustness_path), path_ensemble=load(path_ensemble_path), config=config).to_dict()


def write_decision_report(output: str | Path = "experiments/decision_report.json", **kwargs: Any) -> dict[str, Any]:
    payload = decide_from_files(**kwargs)
    target = Path(output); target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
