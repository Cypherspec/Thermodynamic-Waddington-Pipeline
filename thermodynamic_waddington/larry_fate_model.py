from __future__ import annotations

"""Leakage-resistant real LARRY predictors with train-only feature selection.

The module compares a locked marker model with a training-only panel selected from
the official enriched-gene list. Held-out outcomes are never used for feature
selection, scaling, tuning, or intervention claims.
"""

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .benchmark_championship import evaluate_methods
from .larry_expression import _read_gene_names, stream_larry_matrix
from .larry_predictor import _ground_truth, _selected_day2_indices

DEFAULT_MARKERS = ("Spi1", "Irf8", "Cebpa", "Gata1", "Cebpe", "Mpo", "Csf1r", "Elane", "Gata2")

@dataclass(frozen=True)
class FateModelConfig:
    alpha: float = 1.0
    markers: tuple[str, ...] = DEFAULT_MARKERS
    feature_mode: str = "fixed"
    top_k_enriched: int = 64
    bootstrap_rounds: int = 500
    permutation_rounds: int = 500
    seed: int = 17

    def validate(self) -> None:
        if self.alpha <= 0 or self.top_k_enriched < 2:
            raise ValueError("alpha and top_k_enriched must be positive")
        if len(self.markers) < 2:
            raise ValueError("at least two fixed markers are required")
        if self.feature_mode not in {"fixed", "enriched_train_only"}:
            raise ValueError("feature_mode must be fixed or enriched_train_only")
        if self.bootstrap_rounds < 100 or self.permutation_rounds < 100:
            raise ValueError("resampling rounds must be at least 100")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _feature_panel(root: Path, config: FateModelConfig, train_cells: np.ndarray) -> list[str]:
    genes = set(_read_gene_names(root / "stateFate_inVitro_gene_names.txt.gz"))
    if config.feature_mode == "fixed":
        return [name for name in config.markers if name in genes]
    enriched = []
    path = root / "enriched_genes_IN_VITRO.txt"
    if not path.exists():
        raise FileNotFoundError("official enriched gene panel is required for enriched_train_only mode")
    seen = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1] in genes and parts[1] not in seen:
            seen.add(parts[1]); enriched.append(parts[1])
        if len(enriched) >= config.top_k_enriched * 4:
            break
    if len(enriched) < 2:
        raise ValueError("enriched panel contains too few genes")
    return enriched[:config.top_k_enriched]


def _expression_matrix(root: Path, selected_cells: np.ndarray, markers: Sequence[str]) -> tuple[np.ndarray, list[str]]:
    genes = _read_gene_names(root / "stateFate_inVitro_gene_names.txt.gz")
    wanted = [genes.index(name) for name in markers if name in genes]
    names = [genes[index] for index in wanted]
    index_map = {int(cell): row for row, cell in enumerate(selected_cells.tolist())}
    gene_map = {gene: column for column, gene in enumerate(wanted)}
    values = np.zeros((len(selected_cells), len(wanted)), dtype=float)
    for cell, gene, value in stream_larry_matrix(root, selected_cells=selected_cells.tolist(), selected_genes=wanted):
        row = index_map.get(int(cell)); column = gene_map.get(int(gene))
        if row is not None and column is not None: values[row, column] = float(value)
    if len(names) < 2: raise ValueError("fewer than two requested markers are present")
    return values, names


def _standardize(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    center = train.mean(axis=0); scale = np.maximum(train.std(axis=0), 1e-8)
    return (train - center) / scale, (test - center) / scale


def _ridge_fit(x: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    design = np.column_stack([np.ones(len(x)), x]); penalty = np.eye(design.shape[1]) * alpha; penalty[0, 0] = 0
    return np.linalg.solve(design.T @ design + penalty, design.T @ y)


def _predict(x: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.clip(np.column_stack([np.ones(len(x)), x]) @ weights, 0, 1)


def _run_model(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, test_y: np.ndarray, names: list[str], config: FateModelConfig) -> dict[str, Any]:
    weights = _ridge_fit(train_x, train_y, config.alpha); prediction = _predict(test_x, weights)
    report = evaluate_methods({"TW_fate_model": prediction.tolist()}, test_y.tolist(), benchmark="LARRY day-2 held-out", outcome_definition="published held-out fate probability", holdout_fraction=len(test_y)/(len(train_y)+len(test_y)), bootstrap=config.bootstrap_rounds, permutations=config.permutation_rounds, seed=config.seed, provenance={"feature_mode":config.feature_mode,"features":names,"fit":"ridge training partition only"})
    metric = report.metrics[0].to_dict()
    return {"metrics":metric,"championship":report.to_dict(),"weights":{"intercept":float(weights[0]),**{n:float(w) for n,w in zip(names,weights[1:])}},"prediction":prediction}


def run_real_larry_model(root: str | Path = "data/real/larry", output: str | Path | None = None, config: FateModelConfig | None = None) -> dict[str, Any]:
    config=config or FateModelConfig(); config.validate(); root=Path(root)
    cells, heldout = _selected_day2_indices(root); truth=_ground_truth(root)
    if len(cells)>25000: cells,heldout,truth=cells[:25000],heldout[:25000],truth[:25000]
    train=~heldout; test=heldout & np.isfinite(truth)
    names=_feature_panel(root,config,cells[train]); expression, names=_expression_matrix(root,cells,names)
    train_x,test_x=_standardize(expression[train],expression[test]); train_y,test_y=truth[train],truth[test]
    result=_run_model(train_x,train_y,test_x,test_y,names,config); result.pop('prediction')
    payload={"status":"real_public_heldout_prediction","model":{"name":"TW_fate_model","config":asdict(config),"features":names,"train_only_feature_selection":True,"train_only_standardization":True,"train_only_fit":True},"partitions":{"selected_cells":int(len(cells)),"training_cells":int(train.sum()),"heldout_cells":int(test.sum())},"metrics":result["metrics"],"weights":result["weights"],"championship":result["championship"],"provenance":{"source":"https://github.com/scDiffEq/LARRY-dataset","expression_sha256":_sha256(root/'stateFate_inVitro_normed_counts.mtx.gz'),"gene_list_sha256":_sha256(root/'stateFate_inVitro_gene_names.txt.gz'),"feature_panel_sha256":_sha256(root/'enriched_genes_IN_VITRO.txt') if (root/'enriched_genes_IN_VITRO.txt').exists() else None,"heldout_mask":str(root/'figure5/heldout_mask.npy'),"leakage_policy":"heldout outcomes never enter panel selection or fitting"},"claim_boundary":"Real public held-out predictive evidence only; not causal, wet-lab, or thermodynamic proof."}
    payload["fingerprint"]=hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    if output: Path(output).parent.mkdir(parents=True,exist_ok=True); Path(output).write_text(json.dumps(payload,indent=2,sort_keys=True),encoding='utf-8')
    return payload
