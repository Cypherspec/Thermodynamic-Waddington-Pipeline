from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .arrays import read_matrix
from .synthetic import SyntheticDataset


def read_bundle(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text())
        if "expression" not in payload or "velocity" not in payload:
            raise ValueError("bundle JSON must contain expression and velocity")
        return payload
    raise ValueError("supported bundle input is JSON; use an adapter to convert AnnData/scVelo outputs")


def write_bundle(path: str | Path, payload: dict[str, Any]) -> None:
    Path(path).write_text(json.dumps(payload, separators=(",", ":")))


def save_csv_bundle(directory: str | Path, expression: list[list[float]], velocity: list[list[float]]) -> None:
    from .arrays import write_matrix
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    write_matrix(directory / "expression.csv", expression)
    write_matrix(directory / "velocity.csv", velocity)


def load_expression_velocity(expression_path: str | Path, velocity_path: str | Path) -> tuple[list[list[float]], list[list[float]]]:
    return read_matrix(expression_path), read_matrix(velocity_path)
