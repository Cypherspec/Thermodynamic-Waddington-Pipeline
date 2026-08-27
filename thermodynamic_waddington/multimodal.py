from __future__ import annotations

from math import isfinite
from typing import Any


SUPPORTED_MODALITIES = (
    "proteomics",
    "imaging",
    "metabolomics",
    "viability_assays",
    "genomics",
    "epigenomics",
    "atac_seq",
    "methylation",
    "lipidomics",
    "phosphoproteomics",
    "spatial_transcriptomics",
    "morphology",
    "drug_response",
    "perturbation",
)


def _finite_row(row: Any) -> list[float] | None:
    if not isinstance(row, (list, tuple)):
        return None
    values: list[float] = []
    for value in row:
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        if not isfinite(number):
            return None
        values.append(number)
    return values


def _feature_names(value: Any, width: int) -> list[str]:
    if isinstance(value, list) and len(value) == width:
        return [str(item) for item in value]
    return [f"feature_{index}" for index in range(width)]


def normalize_modalities(metadata: dict[str, Any], n_cells: int) -> dict[str, dict[str, Any]]:
    raw = metadata.get("modalities", {})
    if not isinstance(raw, dict):
        raw = {}
    result: dict[str, dict[str, Any]] = {}
    for modality in SUPPORTED_MODALITIES:
        candidate = raw.get(modality, metadata.get(modality))
        if candidate is None:
            continue
        if isinstance(candidate, dict):
            matrix = candidate.get("values", candidate.get("matrix", candidate.get("data")))
            feature_names = candidate.get("feature_names", candidate.get("features", []))
            source = str(candidate.get("source", "supplied measurement"))
        else:
            matrix = candidate
            feature_names = []
            source = "supplied measurement"
        if not isinstance(matrix, list) or len(matrix) != n_cells:
            continue
        rows = [_finite_row(row) for row in matrix]
        if any(row is None for row in rows):
            continue
        width = len(rows[0] or [])
        if width == 0 or any(len(row or []) != width for row in rows):
            continue
        result[modality] = {
            "values": rows,
            "feature_names": _feature_names(feature_names, width),
            "source": source,
            "measured": True,
            "cell_count": n_cells,
            "feature_count": width,
        }
    return result


def modality_summary(modalities: dict[str, dict[str, Any]], n_cells: int) -> dict[str, Any]:
    return {
        "available": sorted(modalities),
        "missing": [name for name in SUPPORTED_MODALITIES if name not in modalities],
        "cell_count": n_cells,
        "measurement_policy": "Only supplied finite matrices are displayed as measured; absent modalities remain unavailable.",
    }


def cell_modalities(modalities: dict[str, dict[str, Any]], index: int) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for name, payload in modalities.items():
        rows = payload.get("values", [])
        if index < len(rows):
            result[name] = {
                "values": rows[index],
                "feature_names": payload.get("feature_names", []),
                "source": payload.get("source", "supplied measurement"),
                "measured": True,
            }
    return result
