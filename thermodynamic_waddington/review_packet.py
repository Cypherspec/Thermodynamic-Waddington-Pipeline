from __future__ import annotations

"""Generate a reviewer-facing, evidence-bounded research packet.

This is deliberately a claim audit, not a marketing generator: it records what
was measured, what was benchmarked, what failed, and which experiment is still
required before a causal biological claim is permitted.
"""

import hashlib
import json
from pathlib import Path
from typing import Any


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while block := f.read(1024 * 1024):
            h.update(block)
    return h.hexdigest()


def build_review_packet(root: str | Path = ".") -> dict[str, Any]:
    root = Path(root)
    files = [
        "experiments/research_capsule.json",
        "experiments/larry_figure5_championship.json",
        "experiments/larry_tw_fate_model.json",
        "experiments/larry_tw_fate_model_fixed.json",
        "experiments/external_validation_full.json",
        "experiments/wetlab_validation_protocol.json",
        "experiments/evidence_audit.json",
        "experiments/measured_study_manifest.template.json",
        "README.md",
    ]
    artifacts = []
    for name in files:
        path = root / name
        record: dict[str, Any] = {"path": name, "exists": path.exists()}
        if path.exists():
            record.update({"bytes": path.stat().st_size, "sha256": _sha256(path)})
        artifacts.append(record)
    championship = {}
    champ_path = root / "experiments/larry_figure5_championship.json"
    if champ_path.exists():
        championship = json.loads(champ_path.read_text())
    metrics = championship.get("metrics", [])
    return {
        "packet_id": "TW-REVIEW-PACKET-v1",
        "title": "Thermodynamic Waddington: evidence-bounded computational review packet",
        "status": "review_ready_not_biological_discovery",
        "executive_conclusion": "The public LARRY artifact supports a reproducible held-out predictive benchmark. It does not establish causal intervention efficacy, molecular free energy, or wet-lab validation.",
        "validated_observations": {
            "public_benchmark": "Weinreb2020 LARRY Figure 5 held-out outcome arrays",
            "comparative_methods": [row.get("method") for row in metrics],
            "heldout_cells": metrics[0].get("n") if metrics else None,
            "best_published_baseline": max(metrics, key=lambda row: row.get("spearman", float("-inf"))).get("method") if metrics else None,
            "benchmark_metrics": metrics,
        },
        "negative_results": [
            "The Thermodynamic Waddington fixed-panel LARRY predictor remains near-zero Spearman on the held-out fate benchmark.",
            "The enriched train-only panel is retained as a negative result rather than selected post hoc.",
            "No randomized intervention outcome dataset is present.",
        ],
        "claim_gate": {
            "allowed": ["reproducible_public_artifact", "heldout_predictive_comparison"],
            "blocked": ["causal_biological_discovery", "wetlab_validation", "molecular_free_energy_in_kT", "clinical_relevance"],
            "reason": "independent donors, matched randomized controls, lineage/barcode readouts, viability, and orthogonal fate assays are not yet measured in a study bundle",
        },
        "required_next_experiment": {
            "design": "independent donors with randomized intervention assignment and matched vehicle controls",
            "measurements": ["cell-resolved expression", "lineage/barcode identity", "viability", "orthogonal fate assay", "timepoint", "batch and donor metadata"],
            "analysis_lock": "preregister endpoints, exclusions, donor-paired estimator, null permutation, and replication threshold before unblinding",
        },
        "artifacts": artifacts,
        "evidence_audit": {"status": "available" if (root / "experiments/evidence_audit.json").exists() else "missing", "path": "experiments/evidence_audit.json"},
        "reproduction": [
            "python -m unittest discover -s tests -q",
            "python -m compileall -q thermodynamic_waddington",
            "python -m thermodynamic_waddington.capsule_cli --verify experiments/research_capsule.json",
        ],
    }


def write_review_packet(path: str | Path = "experiments/review_packet.json", root: str | Path = ".") -> dict[str, Any]:
    payload = build_review_packet(root)
    payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
