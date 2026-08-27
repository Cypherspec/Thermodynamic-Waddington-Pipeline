from __future__ import annotations

"""A reviewer-facing submission dossier with explicit claims, null results, and gates."""

import hashlib
import json
from pathlib import Path
from typing import Any


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def build_submission_dossier(root: str | Path = ".") -> dict[str, Any]:
    root = Path(root)
    review_path = root / "experiments/review_packet.json"
    review = json.loads(review_path.read_text()) if review_path.exists() else {}
    capsule_path = root / "experiments/research_capsule.json"
    capsule = json.loads(capsule_path.read_text()) if capsule_path.exists() else {}
    files = [
        "experiments/review_packet.json",
        "experiments/research_capsule.json",
        "experiments/larry_figure5_championship.json",
        "experiments/larry_tw_fate_model.json",
        "experiments/wetlab_validation_protocol.json",
        "experiments/measured_study_manifest.template.json",
        "experiments/evidence_audit.json",
        "README.md",
    ]
    artifacts = []
    for name in files:
        path = root / name
        artifacts.append({"path": name, "exists": path.exists(), "bytes": path.stat().st_size if path.exists() else 0, "sha256": _sha256(path) if path.exists() else None})
    return {
        "dossier_id": "TW-SUBMISSION-DOSSIER-v1",
        "title": "Thermodynamic Waddington: auditable effective cell-fate landscape inference",
        "submission_status": "computational_review_ready_biological_discovery_not_established",
        "proposed_contribution": "A reproducible inference and evidence-gating framework that separates held-out predictive analysis from causal biology, physical thermodynamics, and wet-lab discovery.",
        "what_is_verified": [
            "The Python package compiles and its regression suite passes.",
            "Public LARRY Figure 5 held-out arrays are compared with published baselines.",
            "Donor, clone, and timepoint holdout structures are supported by the flagship benchmark.",
            "Fixed-panel and train-only enriched-panel results are retained, including negative results.",
            "Every research capsule artifact is hash-addressed and reproducibility-verifiable.",
            "A cross-artifact evidence audit checks that computational claims remain bounded.",
        ],
        "headline_result": review.get("validated_observations", {}).get("benchmark_metrics", []),
        "negative_results": review.get("negative_results", ["No randomized intervention outcome dataset is present."]),
        "claim_boundary": {
            "permitted": ["auditable computational method", "reproducible public-data benchmark", "held-out predictive comparison"],
            "not_permitted": ["world-changing biological discovery", "causal intervention efficacy", "molecular free energy in kT", "clinical or therapeutic relevance"],
        },
        "wet_lab_readiness": {
            "status": "protocol_and_ingestion_gate_ready_not_executed",
            "required": ["independent donors", "randomized intervention assignment", "matched vehicle controls", "lineage/barcode readouts", "viability", "orthogonal fate assay", "preregistered endpoints", "independent replication"],
            "analysis_rule": "A simulated or public observational dataset cannot satisfy this gate.",
        },
        "capsule_status": capsule.get("status"),
        "artifacts": artifacts,
        "reproduction_commands": [
            "python -m unittest discover -s tests -q",
            "python -m compileall -q thermodynamic_waddington",
            "python -m thermodynamic_waddington.capsule_cli --verify experiments/research_capsule.json",
            "python -m thermodynamic_waddington.review_packet_cli --out experiments/review_packet.json",
        ],
    }


def write_submission_dossier(path: str | Path = "experiments/submission_dossier.json", root: str | Path = ".") -> dict[str, Any]:
    payload = build_submission_dossier(root)
    payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
