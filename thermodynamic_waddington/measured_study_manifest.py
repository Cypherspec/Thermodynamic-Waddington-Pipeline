from __future__ import annotations

"""Schema and gates for ingesting a real measured intervention study."""

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

REQUIRED_FILES = ("cell_metadata", "expression", "lineage", "perturbation", "viability", "orthogonal_fate")
REQUIRED_METADATA = ("cell_id", "donor", "replicate", "intervention", "timepoint")


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def empty_manifest(protocol_id: str = "TW-HEM-FATE-001", preregistration_hash: str = "") -> dict[str, Any]:
    return {
        "manifest_version": "TW-MEASURED-STUDY-v1",
        "status": "awaiting_data",
        "protocol_id": protocol_id,
        "preregistration_hash": preregistration_hash,
        "blinding": {"analysis_locked": False, "intervention_unblinded": False},
        "study": {"target_fate": "Monocyte", "control_intervention": "vehicle", "donors": [], "biological_replicates": [], "timepoints": ["baseline", "early", "commitment", "endpoint"]},
        "files": {name: {"path": "", "sha256": "", "verified": False} for name in REQUIRED_FILES},
        "metadata_schema": {"required_columns": list(REQUIRED_METADATA), "lineage_column": "lineage_id", "viability_column": "viable", "orthogonal_fate_column": "fate"},
        "analysis": {"status": "not_run", "locked_code_fingerprint": "", "output_path": ""},
        "claim_gate": {"minimum_donors": 3, "minimum_cells_per_group": 10, "requires_lineage": True, "requires_orthogonal_fate": True, "requires_viability": True},
        "claim_boundary": "A manifest is a data-integrity record, not evidence of an effect.",
    }


def validate_manifest(manifest: Mapping[str, Any], root: str | Path = ".", require_files: bool = False) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    if not manifest.get("protocol_id"):
        errors.append("missing_protocol_id")
    if not manifest.get("preregistration_hash"):
        errors.append("missing_preregistration_hash")
    blinding = manifest.get("blinding", {})
    if blinding.get("intervention_unblinded") and not blinding.get("analysis_locked"):
        errors.append("intervention_unblinded_before_analysis_lock")
    files = manifest.get("files", {})
    for name in REQUIRED_FILES:
        entry = files.get(name, {})
        if not entry.get("path"):
            errors.append(f"missing_file:{name}")
            continue
        path = Path(root) / str(entry["path"])
        if require_files and not path.exists():
            errors.append(f"file_not_found:{name}")
            continue
        if path.exists() and entry.get("sha256"):
            observed = sha256_file(path)
            if observed != entry["sha256"]:
                errors.append(f"sha256_mismatch:{name}")
        elif path.exists():
            warnings.append(f"sha256_not_recorded:{name}")
    gate = manifest.get("claim_gate", {})
    if int(gate.get("minimum_donors", 0)) < 3:
        errors.append("minimum_donors_below_three")
    if not gate.get("requires_lineage"):
        errors.append("lineage_gate_disabled")
    if not gate.get("requires_viability"):
        errors.append("viability_gate_disabled")
    return {"valid": not errors, "errors": errors, "warnings": warnings, "status": "ready_for_measured_analysis" if not errors else "blocked"}


def write_template(path: str | Path = "experiments/measured_study_manifest.json") -> dict[str, Any]:
    payload = empty_manifest()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload
