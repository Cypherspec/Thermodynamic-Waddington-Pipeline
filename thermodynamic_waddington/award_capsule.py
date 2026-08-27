from __future__ import annotations

"""Reproducible research capsule for external review.

An award-worthy project needs an auditable artifact: exact claims, data hashes,
locked metrics, negative results, environment, and a one-command reproduction
recipe. This capsule packages that evidence without overstating it.
"""

import hashlib
import json
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


@dataclass(frozen=True)
class CapsuleClaim:
    claim_id: str
    statement: str
    status: str
    evidence: tuple[str, ...]
    falsifier: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ResearchCapsule:
    capsule_id: str
    title: str
    created_at: str
    runtime: dict[str, str]
    claims: list[CapsuleClaim]
    artifacts: list[dict[str, Any]]
    commands: list[str]
    limitations: list[str]
    status: str = "review_ready_not_discovery"

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["claims"] = [claim.to_dict() for claim in self.claims]
        payload["fingerprint"] = fingerprint(payload)
        return payload


def fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def _file_record(root: Path, relative: str) -> dict[str, Any]:
    path = root / relative
    if not path.exists():
        return {"path": relative, "exists": False}
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024): digest.update(block)
    return {"path": relative, "exists": True, "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def build_capsule(root: str | Path = ".") -> ResearchCapsule:
    root = Path(root)
    artifacts = [_file_record(root, relative) for relative in (
        "experiments/larry_figure5_championship.json",
        "experiments/larry_tw_fate_model.json",
        "experiments/larry_tw_fate_model_fixed.json",
        "experiments/larry_tw_enriched_train_only.json",
        "experiments/larry_figure5_championship.json",
        "experiments/external_validation_full.json",
        "experiments/preregistration.json",
        "experiments/intervention_registry.json",
        "experiments/measured_study_manifest.template.json",
        "README.md",
    )]
    claims = [
        CapsuleClaim("C1", "The pipeline runs reproducibly against public LARRY lineage/fate artifacts.", "supported_by_public_artifact", ("experiments/external_validation_full.json",), "schema or hash mismatch"),
        CapsuleClaim("C2", "The fixed-panel real LARRY predictor generalizes strongly to held-out fate outcomes.", "falsified_or_not_supported", ("experiments/larry_tw_fate_model.json",), "held-out Spearman remains near zero"),
        CapsuleClaim("C3", "The project has a locked path from prediction to preregistered intervention analysis.", "implemented_not_executed", ("experiments/preregistration.json", "experiments/measured_study_manifest.template.json"), "missing measured randomized study"),
        CapsuleClaim("C4", "The inferred field is physical molecular free energy in kT.", "not_established", ("thermodynamic_waddington/scientific_audit.py",), "requires calibrated controlled work ensemble"),
    ]
    commands = ["python -m unittest discover -s tests -q", "python -m compileall -q thermodynamic_waddington", "python -m thermodynamic_waddington.real_data_cli validate-external", "python -m thermodynamic_waddington.discovery_cli --out experiments/discovery_ledger.json",
        "python -m thermodynamic_waddington.capsule_cli --verify experiments/research_capsule.json"]
    limitations = ["No wet-lab intervention outcomes are present in this workspace.", "Public LARRY expression and fate arrays do not constitute randomized intervention evidence.", "Fixed-panel and train-only enriched-panel negative predictive results are retained rather than hidden.",
        "No wet-lab data are synthesized or presented as biological validation.", "Clinical or therapeutic relevance is not established."]
    return ResearchCapsule("TW-CAPSULE-v1", "Thermodynamic Waddington: auditable cell-fate landscape and intervention validation", datetime.now(timezone.utc).isoformat(), {"python": sys.version.split()[0], "platform": platform.platform()}, claims, artifacts, commands, limitations)


def write_capsule(path: str | Path = "experiments/research_capsule.json", root: str | Path = ".") -> dict[str, Any]:
    capsule = build_capsule(root)
    payload = capsule.to_dict()
    target = Path(path); target.parent.mkdir(parents=True, exist_ok=True); target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload
