from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Sequence


def _load_npy(path: Path):
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("numpy is required for LARRY validation") from exc
    return np.load(path, allow_pickle=True)


def _rank(values: Sequence[float]) -> list[int]:
    return sorted(range(len(values)), key=lambda i: (float(values[i]), i))


def spearman(first: Sequence[float], second: Sequence[float]) -> float:
    n = min(len(first), len(second))
    if n < 2:
        return 0.0
    a = [0] * n
    b = [0] * n
    for rank, index in enumerate(_rank(first[:n])):
        a[index] = rank
    for rank, index in enumerate(_rank(second[:n])):
        b[index] = rank
    ma = sum(a) / n
    mb = sum(b) / n
    numerator = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    denominator = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return numerator / denominator if denominator else 0.0


def validate_larry_fate_ground_truth(root: str | Path = "data/real/larry", predicted: Sequence[Sequence[float]] | None = None, early_indices: Sequence[int] | None = None) -> dict[str, object]:
    root = Path(root)
    fate = _load_npy(root / "clonal_fate_matrix.npy")
    times = _load_npy(root / "timepoints.npy")
    early = _load_npy(root / "early_cells.npy")
    names = [line.strip() for line in (root / "FINAL_fate_names.txt").read_text().splitlines() if line.strip()]
    if fate.ndim != 2 or fate.shape[1] != len(names):
        raise ValueError("LARRY fate matrix and fate-name list disagree")
    dominant = fate.argmax(axis=1).tolist()
    selected = [i for i, flag in enumerate(early.tolist()) if bool(flag)]
    if early_indices is not None:
        selected = [int(i) for i in early_indices]
    result: dict[str, object] = {
        "dataset": "Weinreb2020_LARRY_in_vitro",
        "source": "https://github.com/scDiffEq/LARRY-dataset",
        "cells_with_fate_bias": int(fate.shape[0]),
        "fates": names,
        "early_cell_count": len(selected),
        "timepoint_range": [float(times.min()), float(times.max())],
        "ground_truth_sha256": _sha256(root),
        "validation_status": "ground_truth_loaded",
        "limitations": ["lineage fate bias is an outcome reference, not a measured RNA velocity field", "expression-to-fate alignment requires the matching H5AD and cell IDs", "correlation is not causal validation"],
    }
    if predicted is not None:
        pred = [list(map(float, row)) for row in predicted]
        truth = [fate[i].tolist() for i in selected[:len(pred)]]
        if len(pred) and len(pred[0]) == len(names):
            result["prediction_validation"] = {"cells": len(pred), "per_fate_spearman": {names[j]: spearman([row[j] for row in pred], [row[j] for row in truth]) for j in range(len(names))}, "mean_spearman": sum(spearman([row[j] for row in pred], [row[j] for row in truth]) for j in range(len(names))) / len(names)}
        else:
            result["prediction_validation"] = {"status": "not_comparable", "reason": "prediction width does not equal ten LARRY fates"}
    return result



def validate_larry_heldout_benchmark(root: str | Path = "data/real/larry/figure5") -> dict[str, object]:
    """Score the official held-out neutrophil/monocyte fate benchmark.

    This is outcome prediction validation: it is deliberately separate from
    molecular free-energy validation and uses the authors' published masks and
    smoothed lineage ground truth.
    """
    root = Path(root)
    import numpy as np
    truth = _load_npy(root / "smoothed_groundtruth_from_heldout.npy").astype(float)
    heldout = _load_npy(root / "heldout_mask.npy").astype(bool)
    outliers = _load_npy(root / "outliers_in_SPRING_plot.npy").astype(bool)
    scores: dict[str, dict[str, float | int]] = {}
    for method in ("PBA", "FateID", "WOT"):
        pred = _load_npy(root / f"{method}_predictions.npy").astype(float)
        valid = np.isfinite(pred) & np.isfinite(truth) & (~outliers)
        if len(pred) != len(truth):
            raise ValueError(f"{method} prediction and ground truth lengths disagree")
        residual = pred[valid] - truth[valid]
        scores[method] = {
            "cells": int(valid.sum()),
            "spearman": float(spearman(pred[valid].tolist(), truth[valid].tolist())),
            "mae": float(np.mean(np.abs(residual))),
            "rmse": float(np.sqrt(np.mean(residual ** 2))),
            "heldout_fraction": float(heldout.mean()),
        }
    best = max(scores, key=lambda name: float(scores[name]["spearman"]))
    return {
        "dataset": "Weinreb2020_LARRY_Figure5_heldout",
        "source": "https://github.com/scDiffEq/LARRY-dataset",
        "benchmark": "neutrophil-monocyte fate prediction",
        "methods": scores,
        "best_baseline_by_spearman": best,
        "validation_status": "published_outcome_benchmark_loaded",
        "interpretation": "This validates fate-score prediction against a held-out lineage reference; it does not validate physical free energy or establish causality.",
    }


def _sha256(root: Path) -> str:
    import hashlib
    digest = hashlib.sha256()
    for path in sorted(root.glob("*")):
        if path.is_file():
            digest.update(path.name.encode())
            digest.update(path.read_bytes())
    return digest.hexdigest()


def write_larry_report(path: str | Path = "experiments/larry_validation.json") -> dict[str, object]:
    report = validate_larry_fate_ground_truth()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    return report
