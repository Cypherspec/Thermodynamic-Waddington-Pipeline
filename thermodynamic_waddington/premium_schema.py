from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _tuple(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    return tuple(str(item) for item in value)


@dataclass(frozen=True)
class ProjectManifest:
    project_id: str
    name: str
    version: str = "1.0.0"
    description: str = "Auditable effective cell-fate landscape research platform"
    scientific_boundary: str = "Computational evidence is not a causal biological discovery without measured intervention replication."
    created_at: str = field(default_factory=utc_now)
    capabilities: tuple[str, ...] = (
        "effective_landscape_inference",
        "protocol_aware_work_diagnostics",
        "lineage_and_perturbation_gates",
        "reproducibility_and_provenance",
        "evidence_release_management",
    )
    policies: dict[str, Any] = field(default_factory=lambda: {
        "require_measured_intervention_for_causal_claims": True,
        "require_independent_replication_for_upgrade": True,
        "retain_negative_results": True,
        "allow_visualization": False,
    })

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["capabilities"] = list(self.capabilities)
        return value

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "ProjectManifest":
        return cls(
            project_id=str(value.get("project_id", "tw-project")),
            name=str(value.get("name", "Thermodynamic Waddington")),
            version=str(value.get("version", "1.0.0")),
            description=str(value.get("description", cls.description)),
            scientific_boundary=str(value.get("scientific_boundary", cls.scientific_boundary)),
            created_at=str(value.get("created_at", utc_now())),
            capabilities=_tuple(value.get("capabilities", cls.capabilities)),
            policies=dict(value.get("policies", {})),
        )


@dataclass(frozen=True)
class DatasetRecord:
    dataset_id: str
    name: str
    path: str
    sha256: str
    cells: int
    features: int
    velocity_status: str
    labels_present: bool
    registered_at: str = field(default_factory=utc_now)
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "DatasetRecord":
        return cls(
            dataset_id=str(value["dataset_id"]),
            name=str(value["name"]),
            path=str(value["path"]),
            sha256=str(value["sha256"]),
            cells=int(value["cells"]),
            features=int(value["features"]),
            velocity_status=str(value.get("velocity_status", "unknown")),
            labels_present=bool(value.get("labels_present", False)),
            registered_at=str(value.get("registered_at", utc_now())),
            metadata=dict(value.get("metadata", {})),
        )


@dataclass(frozen=True)
class RunRecord:
    run_id: str
    dataset_id: str
    status: str
    requested_at: str = field(default_factory=utc_now)
    started_at: str | None = None
    finished_at: str | None = None
    output_path: str | None = None
    error: str | None = None
    config: dict[str, Any] = field(default_factory=dict)
    summary: dict[str, Any] = field(default_factory=dict)
    fingerprint: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "RunRecord":
        return cls(
            run_id=str(value["run_id"]),
            dataset_id=str(value["dataset_id"]),
            status=str(value["status"]),
            requested_at=str(value.get("requested_at", utc_now())),
            started_at=value.get("started_at"),
            finished_at=value.get("finished_at"),
            output_path=value.get("output_path"),
            error=value.get("error"),
            config=dict(value.get("config", {})),
            summary=dict(value.get("summary", {})),
            fingerprint=value.get("fingerprint"),
        )


@dataclass(frozen=True)
class EvidenceRecord:
    artifact: str
    status: str
    size_bytes: int
    sha256: str
    modified_at: str
    classification: str
    keys: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["keys"] = list(self.keys)
        return value
