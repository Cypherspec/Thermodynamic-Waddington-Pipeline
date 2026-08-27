from __future__ import annotations

import hashlib
import json
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any


@dataclass(frozen=True)
class RunManifest:
    run_id: str
    created_at: str
    python: str
    platform: str
    input_digest: str
    config: dict[str, Any]
    assumptions: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "created_at": self.created_at,
            "python": self.python,
            "platform": self.platform,
            "input_digest": self.input_digest,
            "config": self.config,
            "assumptions": list(self.assumptions),
        }


def digest_payload(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def make_manifest(payload: Any, config: dict[str, Any], assumptions: tuple[str, ...] = ()) -> RunManifest:
    digest = digest_payload(payload)
    run_id = digest[:16]
    return RunManifest(run_id, datetime.now(timezone.utc).isoformat(), sys.version.split()[0], platform.platform(), digest, config, assumptions)


def audit_record(event: str, payload: Any, manifest: RunManifest) -> dict[str, Any]:
    return {"event": event, "run_id": manifest.run_id, "input_digest": manifest.input_digest, "payload": payload, "timestamp": datetime.now(timezone.utc).isoformat()}
