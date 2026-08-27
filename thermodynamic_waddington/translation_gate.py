from __future__ import annotations

"""Evidence gate for translating a computational result into a biological claim.

This module deliberately separates four levels of evidence: computation, public
observations, measured interventions, and independent replication. A stronger
claim can never be inferred from a weaker level merely because a score is high.
"""

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping


_LEVELS = {
    "computational": 1,
    "observational": 2,
    "measured_intervention": 3,
    "independent_replication": 4,
}


@dataclass(frozen=True)
class EvidenceRecord:
    record_id: str
    level: str
    dataset: str
    provenance: str
    independent_donors: int = 0
    matched_controls: bool = False
    randomized_intervention: bool = False
    lineage_or_barcode: bool = False
    viability: bool = False
    orthogonal_fate_assay: bool = False
    preregistered: bool = False
    held_out: bool = False
    effect_size: float | None = None
    uncertainty: float | None = None
    notes: tuple[str, ...] = ()

    def validate(self) -> list[str]:
        errors: list[str] = []
        if self.level not in _LEVELS:
            errors.append(f"unknown evidence level: {self.level}")
        if not self.dataset.strip():
            errors.append(f"{self.record_id}: dataset is empty")
        if not self.provenance.strip():
            errors.append(f"{self.record_id}: provenance is empty")
        if self.independent_donors < 0:
            errors.append(f"{self.record_id}: independent_donors must be nonnegative")
        if self.level == "measured_intervention":
            for field_name in (
                "matched_controls",
                "randomized_intervention",
                "lineage_or_barcode",
                "viability",
                "orthogonal_fate_assay",
                "preregistered",
            ):
                if not getattr(self, field_name):
                    errors.append(f"{self.record_id}: measured intervention missing {field_name}")
            if self.independent_donors < 2:
                errors.append(f"{self.record_id}: at least two independent donors are required")
        if self.level == "independent_replication" and self.independent_donors < 2:
            errors.append(f"{self.record_id}: replication requires at least two independent donors")
        if self.effect_size is not None and not (-1e6 < self.effect_size < 1e6):
            errors.append(f"{self.record_id}: effect_size is not finite-range")
        if self.uncertainty is not None and self.uncertainty < 0:
            errors.append(f"{self.record_id}: uncertainty must be nonnegative")
        return errors

    def to_dict(self) -> dict[str, Any]:
        return asdict(self) | {"notes": list(self.notes)}


@dataclass(frozen=True)
class TranslationDecision:
    status: str
    attained_level: str
    target_level: str
    blockers: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    records: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self) | {"blockers": list(self.blockers), "warnings": list(self.warnings), "records": list(self.records)}


def _target_level(target: str) -> str:
    aliases = {"discovery": "measured_intervention", "causal": "measured_intervention", "replicated": "independent_replication"}
    return aliases.get(target, target)


def evaluate_translation(records: Iterable[EvidenceRecord], target: str = "measured_intervention") -> TranslationDecision:
    items = list(records)
    target_level = _target_level(target)
    blockers: list[str] = []
    warnings: list[str] = []
    valid: list[EvidenceRecord] = []
    for record in items:
        errors = record.validate()
        if errors:
            blockers.extend(errors)
        else:
            valid.append(record)
    attained_value = max((_LEVELS.get(record.level, 0) for record in valid), default=0)
    attained = next((level for level, value in sorted(_LEVELS.items(), key=lambda pair: pair[1], reverse=True) if value == attained_value), "none")
    required_value = _LEVELS.get(target_level, 99)
    if not items:
        blockers.append("no evidence records supplied")
    if attained_value < required_value:
        blockers.append(f"target level {target_level} is not established by supplied evidence")
    if target_level == "independent_replication":
        eligible = [record for record in valid if record.level == target_level and record.held_out]
        if len(eligible) < 2:
            blockers.append("independent replication requires at least two held-out records")
    if any(record.level in {"computational", "observational"} for record in valid):
        warnings.append("computational or observational evidence cannot establish causality without measured randomized intervention")
    status = "translation_supported" if not blockers else "translation_blocked"
    return TranslationDecision(status, attained, target_level, tuple(dict.fromkeys(blockers)), tuple(dict.fromkeys(warnings)), tuple(record.record_id for record in valid))


def build_translation_report(root: str | Path = ".") -> dict[str, Any]:
    base = Path(root).resolve()
    records: list[EvidenceRecord] = []
    capsule_path = base / "experiments" / "research_capsule.json"
    if capsule_path.is_file():
        records.append(EvidenceRecord("computational-capsule", "computational", "research_capsule", "workspace artifact", held_out=True))
    for relative, record_id in (("experiments/larry_validation.json", "public-larry"), ("experiments/larry_figure5_benchmark.json", "public-larry-heldout")):
        if (base / relative).is_file():
            records.append(EvidenceRecord(record_id, "observational", relative, "public dataset artifact", held_out="heldout" in record_id))
    manifest = base / "experiments" / "measured_study_manifest.json"
    if manifest.is_file():
        records.append(EvidenceRecord("measured-study-manifest", "measured_intervention", "measured_study_manifest", "user-supplied measured study manifest"))
    decision = evaluate_translation(records, target="independent_replication")
    payload = {
        "schema": "thermodynamic-waddington/translation-gate-v1",
        "decision": decision.to_dict(),
        "records": [record.to_dict() for record in records],
        "scientific_boundary": "A public benchmark can establish predictive computational evidence; causal biological discovery requires measured randomized intervention data and independent replication.",
    }
    payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return payload


def write_translation_report(path: str | Path = "experiments/translation_gate.json", root: str | Path = ".") -> dict[str, Any]:
    payload = build_translation_report(root)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
