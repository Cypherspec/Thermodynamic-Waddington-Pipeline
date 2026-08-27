from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable


@dataclass
class ExperimentLedger:
    name: str
    seed: int
    started_at: float = field(default_factory=time.time)
    events: list[dict[str, Any]] = field(default_factory=list)

    def record(self, stage: str, payload: dict[str, Any] | None = None) -> None:
        self.events.append({"stage": stage, "timestamp": time.time(), "payload": payload or {}})

    def fingerprint(self) -> str:
        encoded = json.dumps({"name": self.name, "seed": self.seed, "events": self.events}, sort_keys=True).encode()
        return hashlib.sha256(encoded).hexdigest()

    def manifest(self) -> dict[str, Any]:
        return {"name": self.name, "seed": self.seed, "started_at": self.started_at, "fingerprint": self.fingerprint(), "runtime": {"python": sys.version, "platform": platform.platform()}, "events": self.events}

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.manifest(), indent=2))


def tracked(ledger: ExperimentLedger, stage: str, function: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    ledger.record(stage + ":started")
    try:
        result = function(*args, **kwargs)
    except Exception as error:
        ledger.record(stage + ":failed", {"error": str(error)})
        raise
    ledger.record(stage + ":completed")
    return result
