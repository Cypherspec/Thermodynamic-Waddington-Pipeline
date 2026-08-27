from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from typing import Any


def content_hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def create_manifest(command: str, inputs: dict[str, Any], config: dict[str, Any], seed: int | None) -> dict[str, Any]:
    return {"schema": "twaddington.manifest.v1", "created_at": datetime.now(timezone.utc).isoformat(), "command": command, "input_hash": content_hash(inputs), "config": config, "seed": seed, "runtime": {"python": sys.version, "platform": platform.platform()}, "reproducibility": {"deterministic": seed is not None, "notes": "Randomness is isolated to bootstrap and synthetic generation paths."}}
