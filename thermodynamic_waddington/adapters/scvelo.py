from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..synthetic import SyntheticDataset


def _as_rows(value: Any) -> list[list[float]]:
    if hasattr(value, "tolist"):
        value = value.tolist()
    return [[float(item) for item in row] for row in value]


def load_scvelo_json(path: str | Path) -> SyntheticDataset:
    payload = json.loads(Path(path).read_text())
    expression = payload.get("expression") or payload.get("X") or payload.get("spliced")
    velocity = payload.get("velocity") or payload.get("velocity_vectors")
    if expression is None or velocity is None:
        raise ValueError("scVelo JSON must contain expression/X/spliced and velocity arrays")
    labels = payload.get("labels") or payload.get("cell_types") or ["unknown"] * len(expression)
    ids = payload.get("cell_ids") or [f"cell_{index:05d}" for index in range(len(expression))]
    return SyntheticDataset(_as_rows(expression), _as_rows(velocity), list(labels), list(ids))


def describe_velocity_payload(dataset: SyntheticDataset) -> dict[str, object]:
    cells = len(dataset.expression)
    genes = len(dataset.expression[0]) if dataset.expression else 0
    velocity_dim = len(dataset.velocity[0]) if dataset.velocity else 0
    return {"cells": cells, "expression_features": genes, "velocity_features": velocity_dim, "velocity_space": "caller-defined; align to expression features before fitting", "warning": "This adapter intentionally does not reproduce scVelo preprocessing. Record filtering, moments, and latent-time settings in the manifest."}
