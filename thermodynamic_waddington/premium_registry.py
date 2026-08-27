from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from .premium_schema import DatasetRecord, EvidenceRecord, ProjectManifest, RunRecord, utc_now

SCHEMA = "thermodynamic-waddington/premium-registry-v1"
DEFAULT_MANIFEST = "premium/project.json"
DEFAULT_DATASETS = "premium/datasets.json"
DEFAULT_RUNS = "premium/runs.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid registry file {path}: {error}") from error


def _save(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def default_manifest() -> ProjectManifest:
    return ProjectManifest(project_id="tw-premium", name="Thermodynamic Waddington Premium")


def load_manifest(root: str | Path) -> ProjectManifest:
    base = Path(root).resolve()
    payload = _load(base / DEFAULT_MANIFEST, None)
    return default_manifest() if payload is None else ProjectManifest.from_dict(payload)


def save_manifest(root: str | Path, manifest: ProjectManifest) -> Path:
    target = Path(root).resolve() / DEFAULT_MANIFEST
    _save(target, manifest.to_dict())
    return target


def register_dataset(root: str | Path, dataset_path: str | Path, *, dataset_id: str, name: str | None = None, velocity_status: str = "unknown", labels_present: bool = False, metadata: dict[str, Any] | None = None) -> DatasetRecord:
    base = Path(root).resolve()
    path = Path(dataset_path).resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"dataset is missing or empty: {path}")
    payload = _load(path, None) if path.suffix.lower() == ".json" else None
    cells = len(payload.get("expression", [])) if isinstance(payload, dict) else 0
    features = len(payload.get("expression", [[]])[0]) if isinstance(payload, dict) and payload.get("expression") else 0
    if isinstance(payload, dict):
        labels_present = labels_present or bool(payload.get("labels"))
        velocity_status = str((payload.get("metadata") or {}).get("velocity_status", velocity_status))
    record = DatasetRecord(dataset_id, name or path.stem, _relative(path, base), _sha256(path), cells, features, velocity_status, labels_present, metadata=dict(metadata or {}))
    registry_path = base / DEFAULT_DATASETS
    registry = _load(registry_path, {"schema": "thermodynamic-waddington/dataset-registry-v1", "datasets": []})
    rows = [item for item in registry.get("datasets", []) if item.get("dataset_id") != dataset_id]
    rows.append(record.to_dict())
    _save(registry_path, {"schema": registry.get("schema", "thermodynamic-waddington/dataset-registry-v1"), "updated_at": utc_now(), "datasets": sorted(rows, key=lambda item: item["dataset_id"])})
    return record


def load_datasets(root: str | Path) -> list[DatasetRecord]:
    payload = _load(Path(root).resolve() / DEFAULT_DATASETS, {"datasets": []})
    return [DatasetRecord.from_dict(item) for item in payload.get("datasets", [])]


def register_run(root: str | Path, run: RunRecord) -> Path:
    base = Path(root).resolve()
    target = base / DEFAULT_RUNS
    payload = _load(target, {"schema": "thermodynamic-waddington/run-registry-v1", "runs": []})
    rows = [item for item in payload.get("runs", []) if item.get("run_id") != run.run_id]
    rows.append(run.to_dict())
    _save(target, {"schema": payload.get("schema", "thermodynamic-waddington/run-registry-v1"), "updated_at": utc_now(), "runs": sorted(rows, key=lambda item: item["requested_at"])})
    return target


def load_runs(root: str | Path) -> list[RunRecord]:
    payload = _load(Path(root).resolve() / DEFAULT_RUNS, {"runs": []})
    return [RunRecord.from_dict(item) for item in payload.get("runs", [])]


def audit_artifacts(root: str | Path, paths: Iterable[str | Path] | None = None) -> dict[str, Any]:
    base = Path(root).resolve()
    selected = [Path(item) for item in paths] if paths is not None else [Path("README.md"), Path("pyproject.toml"), Path("experiments/evidence_release.json"), Path("experiments/review_packet.json"), Path("experiments/research_capsule.json")]
    records: list[EvidenceRecord] = []
    missing: list[str] = []
    for relative in selected:
        path = relative if relative.is_absolute() else base / relative
        if not path.is_file() or path.stat().st_size == 0:
            missing.append(_relative(path, base))
            continue
        classification = "evidence" if "experiment" in path.parts else "source"
        keys: tuple[str, ...] = ()
        if path.suffix.lower() == ".json":
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                keys = tuple(sorted(payload.keys())) if isinstance(payload, dict) else ()
            except json.JSONDecodeError:
                classification = "malformed"
        modified = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
        records.append(EvidenceRecord(_relative(path, base), "present", path.stat().st_size, _sha256(path), modified, classification, keys))
    digest = hashlib.sha256()
    for record in records:
        digest.update(json.dumps(record.to_dict(), sort_keys=True).encode())
    return {
        "schema": SCHEMA,
        "status": "complete" if not missing else "incomplete",
        "root": str(base),
        "artifacts": [record.to_dict() for record in records],
        "missing": missing,
        "fingerprint": digest.hexdigest(),
        "runtime": {"python": sys.version.split()[0], "platform": platform.platform()},
        "scientific_boundary": "Artifact integrity and reproducibility are not evidence of a causal biological effect.",
    }
