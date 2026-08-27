from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Sequence


def read_lineage_table(path: str | Path) -> dict[str, list[str]]:
    path = Path(path)
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {}
    keys = list(rows[0])
    cell_key = next((key for key in keys if key.lower() in {"cell", "cell_id", "barcode", "id"}), keys[0])
    fate_key = next((key for key in keys if "fate" in key.lower() or "lineage" in key.lower() or "type" in key.lower()), keys[1] if len(keys) > 1 else keys[0])
    result: dict[str, list[str]] = {}
    for row in rows:
        result.setdefault(row[cell_key], []).append(row[fate_key])
    return result


def fate_probabilities(lineage: dict[str, Sequence[str]]) -> dict[str, dict[str, float]]:
    output = {}
    for cell, fates in lineage.items():
        counts: dict[str, int] = {}
        for fate in fates:
            counts[str(fate)] = counts.get(str(fate), 0) + 1
        total = max(1, sum(counts.values()))
        output[cell] = {fate: count / total for fate, count in counts.items()}
    return output


def write_lineage_summary(path: str | Path, lineage: dict[str, Sequence[str]]) -> None:
    Path(path).write_text(json.dumps(fate_probabilities(lineage), indent=2, sort_keys=True))
