from __future__ import annotations

"""Cross-artifact integrity audit for evidence-bounded research releases."""

import hashlib
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class AuditCheck:
    name: str
    passed: bool
    observed: Any
    required: Any
    interpretation: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _load(root: Path, relative: str) -> dict[str, Any] | None:
    path = root / relative
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_evidence_audit(root: str | Path = ".") -> dict[str, Any]:
    root = Path(root).resolve()
    review = _load(root, "experiments/review_packet.json") or {}
    submission = _load(root, "experiments/submission_dossier.json") or {}
    capsule = _load(root, "experiments/research_capsule.json") or {}
    falsification = _load(root, "experiments/falsification_report.json") or {}
    reproduction = _load(root, "experiments/reproduction_run.json") or {}
    metrics = review.get("validated_observations", {}).get("benchmark_metrics", [])
    checks = [
        AuditCheck("software_tests", all(step.get("passed", False) for step in reproduction.get("steps", [])), reproduction.get("steps", []), "all reproduction steps pass", "Software integrity only."),
        AuditCheck("falsification_suite", falsification.get("status") == "passed", falsification.get("status"), "passed", "Synthetic adversarial checks only."),
        AuditCheck("heldout_benchmark_present", bool(metrics) and all(int(row.get("n", 0)) > 0 for row in metrics), [row.get("n") for row in metrics], "nonempty held-out metrics", "Public predictive comparison is available."),
        AuditCheck("negative_results_retained", len(review.get("negative_results", [])) >= 3, len(review.get("negative_results", [])), ">=3", "Negative evidence is not hidden."),
        AuditCheck("wetlab_gate_blocked", "wetlab_validation" in review.get("claim_gate", {}).get("blocked", []), review.get("claim_gate", {}).get("blocked", []), "wetlab_validation blocked", "No software artifact can substitute for measured intervention biology."),
        AuditCheck("causal_gate_blocked", "causal_biological_discovery" in review.get("claim_gate", {}).get("blocked", []), review.get("claim_gate", {}).get("blocked", []), "causal_biological_discovery blocked", "Observational/public benchmarks are not causal evidence."),
        AuditCheck("capsule_fingerprint_present", len(capsule.get("fingerprint", "")) == 64, capsule.get("fingerprint", ""), "64 hex characters", "Release provenance is hash-addressed."),
        AuditCheck("submission_is_bounded", submission.get("submission_status") == "computational_review_ready_biological_discovery_not_established", submission.get("submission_status"), "bounded review status", "The dossier does not overclaim."),
    ]
    artifact_hashes = {}
    for relative in ("experiments/review_packet.json", "experiments/submission_dossier.json", "experiments/research_capsule.json", "experiments/falsification_report.json", "experiments/reproduction_run.json"):
        path = root / relative
        if path.exists():
            artifact_hashes[relative] = {"bytes": path.stat().st_size, "sha256": _sha256(path)}
    payload: dict[str, Any] = {
        "audit_id": "TW-EVIDENCE-AUDIT-v1",
        "status": "computationally_reproducible_evidence_bounded" if all(check.passed for check in checks) else "audit_incomplete",
        "checks": [check.to_dict() for check in checks],
        "passed": sum(check.passed for check in checks),
        "total": len(checks),
        "artifact_hashes": artifact_hashes,
        "claim_boundary": "This audit verifies software integrity, artifact consistency, adversarial synthetic checks, and public predictive evidence. It cannot establish a biological discovery without measured randomized interventions, independent donors, matched controls, lineage/barcode data, viability, and orthogonal fate assays.",
    }
    payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return payload


def write_evidence_audit(output: str | Path = "experiments/evidence_audit.json", root: str | Path = ".") -> dict[str, Any]:
    payload = build_evidence_audit(root)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
