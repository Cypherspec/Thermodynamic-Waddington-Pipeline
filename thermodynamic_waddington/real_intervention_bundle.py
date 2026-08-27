from __future__ import annotations

"""Strict ingestion of measured perturbation studies.

The loader accepts a cell-level JSONL/JSON table and refuses to call it real
measured evidence unless the study manifest, preregistration hash, required
assays, and cell-level fields all agree.
"""

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

from .causal_validation import StudyCell, StudyConfig, StudyReport, analyze_study
from .measured_study_manifest import validate_manifest


@dataclass(frozen=True)
class IngestResult:
    status: str
    cells: int
    donors: int
    interventions: tuple[str, ...]
    blockers: tuple[str, ...]
    warnings: tuple[str, ...]
    manifest_status: str
    claim_status: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _rows(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() == ".jsonl":
        return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload.get("cells", payload.get("observations", payload)) if isinstance(payload, dict) else payload
    if not isinstance(rows, list):
        raise ValueError("data file must contain a JSON list, cells, or observations")
    return rows


def ingest_measured_study(data_path: str | Path, manifest_path: str | Path, output_path: str | Path | None = None, *, require_files: bool = False) -> dict[str, Any]:
    data_path = Path(data_path)
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    manifest_result = validate_manifest(manifest, Path(manifest_path).parent, require_files=require_files)
    blockers = list(manifest_result["errors"])
    warnings = list(manifest_result["warnings"])
    if manifest.get("status") in {"awaiting_data", "synthetic", "simulated"}:
        blockers.append("manifest_not_marked_as_measured")
    rows = _rows(data_path)
    cells = [StudyCell.from_dict(row) for row in rows]
    study = manifest.get("study", {})
    gate = manifest.get("claim_gate", {})
    config = StudyConfig(
        target_fate=str(study.get("target_fate", "Monocyte")),
        control_intervention=str(study.get("control_intervention", "vehicle")),
        minimum_donors=int(gate.get("minimum_donors", 3)),
        minimum_cells_per_group=int(gate.get("minimum_cells_per_group", 10)),
        require_lineage=bool(gate.get("requires_lineage", True)),
        require_orthogonal_fate=bool(gate.get("requires_orthogonal_fate", True)),
        require_measured_fate=True,
        require_perturbation_capture=True,
        require_velocity=bool(gate.get("requires_velocity", False)),
        require_expression=bool(gate.get("requires_expression", True)),
        endpoint_timepoint=str(study.get("endpoint_timepoint", "endpoint")),
    )
    report = analyze_study(cells, config, manifest.get("protocol_id", "unregistered"))
    blockers.extend(report.blockers)
    status = "accepted_for_locked_analysis" if not blockers and report.status != "invalid_study_schema" else "blocked"
    result = IngestResult(status, len(cells), len({cell.donor for cell in cells}), tuple(sorted({cell.intervention for cell in cells})), tuple(sorted(set(blockers))), tuple(sorted(set(warnings))), manifest_result["status"], report.status)
    payload = {"ingest": result.to_dict(), "study_report": report.to_dict(), "claim_boundary": "Measured-data ingestion is not proof of an effect; only the locked analysis of a valid randomized study can produce a candidate result."}
    if output_path:
        output_path = Path(output_path); output_path.parent.mkdir(parents=True, exist_ok=True); output_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload
