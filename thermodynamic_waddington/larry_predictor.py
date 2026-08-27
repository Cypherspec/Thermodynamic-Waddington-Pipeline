from __future__ import annotations

"""Leakage-resistant LARRY day-2 fate prediction from aligned sparse expression.

This baseline deliberately uses only day-2 expression and the published training
portion of LARRY fate bias. It is an outcome predictor, not a thermodynamic proof.
"""

import json
from pathlib import Path
from typing import Any

import numpy as np

from .larry_expression import _read_metadata, _read_gene_names, stream_larry_matrix
from .benchmark_championship import evaluate_methods


def _load_masks(root: Path) -> dict[str, np.ndarray]:
    figure = root / "figure5"
    return {
        "timepoints": np.load(figure / "timepoints.npy").astype(float),
        "trajectory": np.load(figure / "neutrophil_monocyte_trajectory_mask.npy").astype(bool),
        "early": np.load(figure / "early_cells.npy").astype(bool),
        "heldout": np.load(figure / "heldout_mask.npy").astype(bool),
        "outliers": np.load(figure / "outliers_in_SPRING_plot.npy").astype(bool),
    }


def _selected_day2_indices(root: Path) -> tuple[np.ndarray, np.ndarray]:
    masks = _load_masks(root)
    day2 = masks["timepoints"] == 2.0
    nm = masks["trajectory"] & day2
    early = masks["early"]
    heldout = masks["heldout"]
    outliers = masks["outliers"]
    early_nm = early
    eligible = ~outliers
    return np.flatnonzero(nm)[eligible], heldout[eligible]


def _ground_truth(root: Path) -> np.ndarray:
    masks = _load_masks(root)
    truth = np.load(root / "figure5" / "smoothed_groundtruth_from_heldout.npy").astype(float)
    nm = masks["trajectory"] & (masks["timepoints"] == 2.0)
    eligible = ~masks["outliers"]
    return truth[eligible]


def _expression_features(root: Path, selected_cells: np.ndarray, genes: list[str], feature_names: tuple[str, ...] = ("Spi1", "Irf8", "Cebpa", "Gata1", "Cebpe", "Mpo", "Csf1r", "Elane", "Gata2")) -> tuple[np.ndarray, list[str]]:
    wanted = [genes.index(name) for name in feature_names if name in genes]
    names = [genes[index] for index in wanted]
    if not wanted:
        raise ValueError("none of the requested marker genes were found in the LARRY gene list")
    index_map = {int(cell): i for i, cell in enumerate(selected_cells.tolist())}
    values = np.zeros((len(selected_cells), len(wanted)), dtype=float)
    for cell, gene, value in stream_larry_matrix(root, selected_cells=selected_cells.tolist(), selected_genes=wanted):
        values[index_map[cell], wanted.index(gene)] = value
    return values, names


def _safe_zscore(values: np.ndarray) -> np.ndarray:
    center = values.mean(axis=0, keepdims=True)
    spread = values.std(axis=0, keepdims=True)
    return (values - center) / np.maximum(spread, 1e-8)


def predict_larry_day2(root: str | Path = "data/real/larry", output: str | Path | None = None) -> dict[str, Any]:
    root = Path(root)
    selected_cells, heldout = _selected_day2_indices(root)
    truth = _ground_truth(root)
    genes = _read_gene_names(root / "stateFate_inVitro_gene_names.txt.gz")
    matrix, feature_names = _expression_features(root, selected_cells, genes)
    z = _safe_zscore(matrix)
    direction = np.array([1.0 if name in {"Spi1", "Irf8", "Cebpa", "Cebpe", "Mpo", "Csf1r", "Elane"} else -1.0 for name in feature_names])
    predictor = 1.0 / (1.0 + np.exp(-z @ direction / max(1.0, np.sqrt(len(direction)))))
    train = ~heldout
    test = heldout & np.isfinite(truth)
    if train.sum() < 20 or test.sum() < 20:
        raise ValueError("insufficient aligned LARRY train/test cells")
    report = evaluate_methods({"TW_expression_marker_proxy": predictor[test].tolist()}, truth[test].tolist(), benchmark="LARRY day-2 held-out expression proxy", holdout_fraction=float(test.mean()), bootstrap=500, permutations=500, provenance={"source": "https://github.com/scDiffEq/LARRY-dataset", "selected_cells": int(len(selected_cells)), "train_cells": int(train.sum()), "test_cells": int(test.sum()), "features": feature_names, "model": "pre-specified signed marker proxy", "velocity_used": False})
    payload = {"status": "heldout_predictive_benchmark", "train_cells": int(train.sum()), "test_cells": int(test.sum()), "features": feature_names, "prediction_summary": {"mean": float(predictor.mean()), "std": float(predictor.std())}, "championship": report.to_dict(), "claim_boundary": "This is a pre-specified expression-only predictive proxy. It is not Thermodynamic Waddington's thermodynamic estimate, not causal evidence, and not a wet-lab result."}
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        Path(output).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload
