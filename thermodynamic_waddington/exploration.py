from __future__ import annotations

import math
from typing import Any, Sequence


def _numbers(values: Any) -> list[float]:
    if not isinstance(values, (list, tuple)):
        return []
    result: list[float] = []
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number):
            result.append(number)
    return result


def _stats(values: Any) -> dict[str, float | int]:
    numbers = _numbers(values)
    if not numbers:
        return {"count": 0, "mean": 0.0, "sd": 0.0, "minimum": 0.0, "maximum": 0.0, "nonzero_fraction": 0.0}
    average = sum(numbers) / len(numbers)
    variance = sum((value - average) ** 2 for value in numbers) / max(1, len(numbers) - 1)
    return {
        "count": len(numbers),
        "mean": average,
        "sd": math.sqrt(max(0.0, variance)),
        "minimum": min(numbers),
        "maximum": max(numbers),
        "nonzero_fraction": sum(abs(value) > 1e-12 for value in numbers) / len(numbers),
    }


def _top_features(values: Any, names: Any, limit: int = 8) -> list[dict[str, float | str]]:
    numbers = _numbers(values)
    labels = [str(name) for name in names] if isinstance(names, list) else []
    ranked = sorted(enumerate(numbers), key=lambda item: abs(item[1]), reverse=True)[:limit]
    return [{"feature": labels[index] if index < len(labels) else f"feature_{index}", "value": value} for index, value in ranked]


def cell_measurement_ledger(modalities: Any, index: int) -> list[dict[str, Any]]:
    from .multimodal import SUPPORTED_MODALITIES

    supplied = modalities if isinstance(modalities, dict) else {}
    ledger: list[dict[str, Any]] = []
    for name in SUPPORTED_MODALITIES:
        payload = supplied.get(name)
        if not isinstance(payload, dict):
            ledger.append({
                "modality": name,
                "status": "unavailable",
                "source": "not supplied",
                "feature_count": 0,
                "observed_count": 0,
                "summary": _stats(None),
                "top_features": [],
                "caveat": "No measurement matrix was supplied for this modality; it is not inferred from RNA.",
            })
            continue
        rows = payload.get("values", [])
        values = rows[index] if isinstance(rows, list) and index < len(rows) else None
        stats = _stats(values)
        features = payload.get("feature_names", [])
        ledger.append({
            "modality": str(name),
            "status": "measured" if stats["count"] else "unavailable",
            "source": str(payload.get("source", "supplied measurement")),
            "feature_count": len(values) if isinstance(values, list) else 0,
            "observed_count": stats["count"],
            "summary": stats,
            "top_features": _top_features(values, features),
            "caveat": "Supplied measurement; interpretation depends on assay design and normalization." if stats["count"] else "No finite measurement supplied for this cell.",
        })
    return ledger


def _population_z(values: Sequence[float], index: int) -> tuple[float, list[float]]:
    if not values or index >= len(values):
        return 0.0, []
    width = len(values[index])
    columns = [[float(row[column]) for row in values if column < len(row) and math.isfinite(float(row[column]))] for column in range(width)]
    row = [float(value) for value in values[index]]
    z_scores: list[float] = []
    for column, value in enumerate(row):
        column_values = columns[column] if column < len(columns) else []
        if not column_values:
            z_scores.append(0.0)
            continue
        average = sum(column_values) / len(column_values)
        variance = sum((item - average) ** 2 for item in column_values) / max(1, len(column_values) - 1)
        z_scores.append((value - average) / max(1e-9, math.sqrt(variance)))
    return math.sqrt(sum(score * score for score in z_scores) / max(1, len(z_scores))), z_scores


def population_context(modalities: Any, index: int) -> list[dict[str, Any]]:
    if not isinstance(modalities, dict):
        return []
    context: list[dict[str, Any]] = []
    for name in sorted(modalities):
        payload = modalities.get(name)
        if not isinstance(payload, dict) or not isinstance(payload.get("values"), list):
            continue
        score, z_scores = _population_z(payload["values"], index)
        names = payload.get("feature_names", [])
        ranked = sorted(enumerate(z_scores), key=lambda item: abs(item[1]), reverse=True)[:8]
        context.append({
            "modality": str(name),
            "population_distance": score,
            "outliers": [{"feature": names[column] if isinstance(names, list) and column < len(names) else f"feature_{column}", "z": z} for column, z in ranked],
            "interpretation": "cell-level deviation from the supplied population for this modality; not a diagnosis",
        })
    return context


def tour_scenes(cells: Sequence[dict[str, Any]], selected_index: int | None = None, limit: int = 7) -> list[dict[str, Any]]:
    if not cells:
        return []
    selected = next((cell for cell in cells if int(cell.get("index", -1)) == selected_index), None) if selected_index is not None else None
    ordered = [selected] if selected else []
    remaining = [cell for cell in cells if cell is not selected]
    ordered.extend(sorted(remaining, key=lambda cell: float(cell.get("energy", 0.0)))[: max(0, limit - len(ordered))])
    scenes: list[dict[str, Any]] = []
    for number, cell in enumerate(ordered[:limit], 1):
        scenes.append({
            "number": number,
            "title": "Selected cell" if number == 1 and selected else ("Lowest sampled basin" if number == 1 else "Landscape waypoint"),
            "cell_index": int(cell.get("index", 0)),
            "cell_id": str(cell.get("id", "cell")),
            "narration": f"{cell.get('id', 'Cell')} occupies the {cell.get('state_class', 'sampled state')} at effective F / kT {float(cell.get('energy', 0.0)):.3f}.",
            "evidence": ["expression", "RNA velocity" if cell.get("velocity_observed", False) else "velocity unavailable", "basin geometry"],
            "grounding": "Computed from the fitted report payload; no external narrative or unobserved measurement is added.",
        })
    return scenes
