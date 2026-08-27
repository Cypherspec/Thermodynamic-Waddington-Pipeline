from __future__ import annotations

import hashlib
import json
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any


@dataclass(frozen=True)
class Provenance:
    run_id: str
    created_at: str
    python: str
    platform: str
    input_fingerprint: str
    stages: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {"run_id": self.run_id, "created_at": self.created_at, "python": self.python, "platform": self.platform, "input_fingerprint": self.input_fingerprint, "stages": list(self.stages)}


def fingerprint(payload: Any) -> str:
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()
    return hashlib.sha256(serialized).hexdigest()[:20]


def make_provenance(payload: Any, stages: list[str] | tuple[str, ...]) -> Provenance:
    digest = fingerprint(payload)
    return Provenance(digest, datetime.now(timezone.utc).isoformat(), sys.version.split()[0], platform.platform(), digest, tuple(stages))


def merge_diagnostics(*records: dict[str, Any]) -> dict[str, Any]:
    merged: dict[str, Any] = {}
    for index, record in enumerate(records):
        for key, value in record.items():
            merged[f"stage_{index}.{key}"] = value
    return merged
