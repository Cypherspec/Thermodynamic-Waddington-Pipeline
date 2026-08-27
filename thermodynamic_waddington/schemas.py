from __future__ import annotations

from typing import Any

BUNDLE_SCHEMA: dict[str, Any] = {"schema": "twaddington.bundle.v1", "required": ["expression", "velocity"], "optional": ["labels", "cell_ids", "metadata"]}
FIT_SCHEMA: dict[str, Any] = {"schema": "twaddington.fit.v1", "required": ["free_energy", "uncertainty", "attractors", "barriers", "diagnostics"], "optional": ["config", "metadata"]}


def validate_schema(payload: dict[str, Any], schema: dict[str, Any]) -> None:
    missing = [key for key in schema["required"] if key not in payload]
    if missing:
        raise ValueError(f"missing required fields for {schema['schema']}: {', '.join(missing)}")
