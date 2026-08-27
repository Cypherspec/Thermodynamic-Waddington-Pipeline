from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable


@dataclass
class ExperimentLedger:
    name: str
    started_at: float = field(default_factory=time.time)
    events: list[dict[str, Any]] = field(default_factory=list)

    def record(self, stage: str, **payload: Any) -> None:
        self.events.append({"stage": stage, "elapsed_seconds": time.time() - self.started_at, **payload})

    def run(self, stage: str, operation: Callable[[], Any]) -> Any:
        start = time.perf_counter()
        self.record(stage, status="started")
        try:
            value = operation()
        except Exception as error:
            self.record(stage, status="failed", error=type(error).__name__, message=str(error))
            raise
        self.record(stage, status="complete", duration_seconds=time.perf_counter() - start)
        return value

    def export(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps({"name": self.name, "started_at": self.started_at, "events": self.events}, indent=2))


def analyze_experiment(energies: list[float], edges: list[dict], labels: list[str], temperature: float) -> dict[str, Any]:
    finite = [value for value in energies if value == value]
    grouped: dict[str, list[float]] = {}
    for value, label in zip(energies, labels):
        grouped.setdefault(label, []).append(value)
    return {
        "n_cells": len(energies),
        "n_edges": len(edges),
        "energy_range": [min(finite), max(finite)] if finite else [0.0, 0.0],
        "temperature": temperature,
        "label_groups": {label: {"n": len(values), "mean_energy": sum(values) / len(values)} for label, values in grouped.items() if values},
    }
