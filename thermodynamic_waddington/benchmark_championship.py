from __future__ import annotations

"""Leakage-resistant comparison of fate predictors on published held-out arrays.

This module deliberately treats published arrays as outcomes, not as training data.
It provides paired metrics, bootstrap intervals, permutation tests, calibration
summaries, and machine-readable claim boundaries. It does not make a predictor
causal merely because it wins a benchmark.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class PredictionMetrics:
    method: str
    n: int
    spearman: float
    pearson: float
    mae: float
    rmse: float
    log_loss: float
    brier: float
    calibration_gap: float
    top_decile_recall: float
    finite: bool
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class PairedComparison:
    left: str
    right: str
    n: int
    spearman_delta: float
    mae_delta: float
    bootstrap_low: float
    bootstrap_high: float
    permutation_p: float
    winner: str
    exchangeability_assumption: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ChampionshipReport:
    benchmark: str
    outcome_definition: str
    metrics: list[PredictionMetrics] = field(default_factory=list)
    comparisons: list[PairedComparison] = field(default_factory=list)
    holdout_fraction: float = 0.0
    status: str = "diagnostic_only"
    claim_boundary: str = "Predictive benchmark performance is not causal evidence or molecular free-energy measurement."
    provenance: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "benchmark": self.benchmark,
            "outcome_definition": self.outcome_definition,
            "metrics": [x.to_dict() for x in self.metrics],
            "comparisons": [x.to_dict() for x in self.comparisons],
            "holdout_fraction": self.holdout_fraction,
            "status": self.status,
            "claim_boundary": self.claim_boundary,
            "provenance": self.provenance,
        }


def _finite_pairs(prediction: Sequence[float], outcome: Sequence[float]) -> tuple[list[float], list[float]]:
    pairs = []
    for left, right in zip(prediction, outcome):
        try:
            x, y = float(left), float(right)
        except (TypeError, ValueError):
            continue
        if math.isfinite(x) and math.isfinite(y):
            pairs.append((x, y))
    return [x for x, _ in pairs], [y for _, y in pairs]


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def _rank(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: (values[i], i))
    ranks = [0.0] * len(values)
    pos = 0
    while pos < len(order):
        end = pos + 1
        while end < len(order) and values[order[end]] == values[order[pos]]:
            end += 1
        rank = 0.5 * (pos + end - 1) + 1.0
        for j in range(pos, end):
            ranks[order[j]] = rank
        pos = end
    return ranks


def _pearson(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) < 2:
        return float("nan")
    lm, rm = _mean(left), _mean(right)
    a = [x - lm for x in left]
    b = [x - rm for x in right]
    den = math.sqrt(sum(x * x for x in a) * sum(y * y for y in b))
    return sum(x * y for x, y in zip(a, b)) / den if den > 0 else 0.0


def _spearman(left: Sequence[float], right: Sequence[float]) -> float:
    return _pearson(_rank(left), _rank(right))


def _clip_probability(value: float, epsilon: float = 1e-8) -> float:
    return min(1.0 - epsilon, max(epsilon, value))


def _normalize_outcome(outcome: Sequence[float]) -> list[float]:
    if not outcome:
        return []
    lo, hi = min(outcome), max(outcome)
    if lo >= 0.0 and hi <= 1.0:
        return [float(x) for x in outcome]
    if hi == lo:
        return [0.5] * len(outcome)
    return [(float(x) - lo) / (hi - lo) for x in outcome]


def _calibration_gap(prediction: Sequence[float], outcome: Sequence[float], bins: int = 10) -> float:
    if not prediction:
        return float("nan")
    total = 0.0
    weight = 0
    for bucket in range(bins):
        low = bucket / bins
        high = (bucket + 1) / bins
        indices = [i for i, p in enumerate(prediction) if low <= p < high or (bucket == bins - 1 and p == high)]
        if indices:
            total += abs(_mean([prediction[i] for i in indices]) - _mean([outcome[i] for i in indices])) * len(indices)
            weight += len(indices)
    return total / weight if weight else float("nan")


def _top_recall(prediction: Sequence[float], outcome: Sequence[float], fraction: float = 0.1) -> float:
    if not prediction:
        return float("nan")
    n = max(1, int(math.ceil(len(prediction) * fraction)))
    predicted = set(sorted(range(len(prediction)), key=lambda i: prediction[i], reverse=True)[:n])
    actual = set(sorted(range(len(outcome)), key=lambda i: outcome[i], reverse=True)[:n])
    return len(predicted & actual) / len(actual)


def evaluate_prediction(method: str, prediction: Sequence[float], outcome: Sequence[float]) -> PredictionMetrics:
    pred, obs = _finite_pairs(prediction, outcome)
    warnings: list[str] = []
    if len(pred) < 20:
        warnings.append("fewer_than_20_finite_pairs")
    if len(pred) != len(prediction) or len(obs) != len(outcome):
        warnings.append("nonfinite_or_misaligned_values_dropped")
    probabilities = [_clip_probability(x if 0.0 <= x <= 1.0 else (x - min(pred)) / max(1e-12, max(pred) - min(pred))) for x in pred] if pred else []
    binary_obs = _normalize_outcome(obs)
    log_loss = _mean([-y * math.log(p) - (1.0 - y) * math.log(1.0 - p) for p, y in zip(probabilities, binary_obs)]) if pred else float("nan")
    errors = [x - y for x, y in zip(pred, obs)]
    return PredictionMetrics(
        method=method,
        n=len(pred),
        spearman=_spearman(pred, obs),
        pearson=_pearson(pred, obs),
        mae=_mean([abs(x) for x in errors]),
        rmse=math.sqrt(_mean([x * x for x in errors])),
        log_loss=log_loss,
        brier=_mean([(p - y) ** 2 for p, y in zip(probabilities, binary_obs)]) if pred else float("nan"),
        calibration_gap=_calibration_gap(probabilities, binary_obs),
        top_decile_recall=_top_recall(pred, obs),
        finite=bool(pred) and all(math.isfinite(x) for x in pred + obs),
        warnings=tuple(warnings),
    )


def _paired_delta(left: Sequence[float], right: Sequence[float], outcome: Sequence[float]) -> float:
    l, r, y = [], [], []
    for a, b, c in zip(left, right, outcome):
        if all(math.isfinite(float(x)) for x in (a, b, c)):
            l.append(float(a)); r.append(float(b)); y.append(float(c))
    return _spearman(l, y) - _spearman(r, y) if l else float("nan")


def paired_compare(left_name: str, left: Sequence[float], right_name: str, right: Sequence[float], outcome: Sequence[float], bootstrap: int = 1000, seed: int = 17, permutations: int = 1000) -> PairedComparison:
    triples = [(float(a), float(b), float(c)) for a, b, c in zip(left, right, outcome) if all(math.isfinite(float(x)) for x in (a, b, c))]
    if len(triples) < 2:
        raise ValueError("paired comparison requires at least two finite aligned observations")
    rng = random.Random(seed)
    observed = _paired_delta(left, right, outcome)
    bootstrap_values: list[float] = []
    for _ in range(max(100, bootstrap)):
        sample = [triples[rng.randrange(len(triples))] for _ in triples]
        bootstrap_values.append(_paired_delta([x[0] for x in sample], [x[1] for x in sample], [x[2] for x in sample]))
    bootstrap_values.sort()
    lo = bootstrap_values[max(0, int(0.025 * len(bootstrap_values)))]
    hi = bootstrap_values[min(len(bootstrap_values) - 1, int(0.975 * len(bootstrap_values)))]
    extreme = 0
    for _ in range(max(100, permutations)):
        swapped_left, swapped_right = [], []
        for a, b, _ in triples:
            if rng.random() < 0.5:
                swapped_left.append(a); swapped_right.append(b)
            else:
                swapped_left.append(b); swapped_right.append(a)
        delta = _paired_delta(swapped_left, swapped_right, [x[2] for x in triples])
        if abs(delta) >= abs(observed):
            extreme += 1
    p_value = (extreme + 1.0) / (max(100, permutations) + 1.0)
    winner = left_name if observed > 0 else right_name if observed < 0 else "tie"
    return PairedComparison(left_name, right_name, len(triples), observed, _mean([abs(a - c) - abs(b - c) for a, b, c in triples]), lo, hi, p_value, winner, "Within-cell paired exchangeability under the null")


def evaluate_methods(predictions: Mapping[str, Sequence[float]], outcome: Sequence[float], *, benchmark: str = "unnamed", outcome_definition: str = "published held-out outcome", holdout_fraction: float = 1.0, bootstrap: int = 1000, permutations: int = 1000, seed: int = 17, provenance: Mapping[str, Any] | None = None) -> ChampionshipReport:
    if not predictions:
        raise ValueError("at least one prediction method is required")
    metrics = [evaluate_prediction(name, values, outcome) for name, values in predictions.items()]
    comparisons: list[PairedComparison] = []
    names = list(predictions)
    for i, left in enumerate(names):
        for right in names[i + 1:]:
            comparisons.append(paired_compare(left, predictions[left], right, predictions[right], outcome, bootstrap=bootstrap, permutations=permutations, seed=seed + i))
    return ChampionshipReport(benchmark, outcome_definition, metrics, comparisons, float(holdout_fraction), "validated_predictive_comparison", "Predictive benchmark performance is not causal evidence or molecular free-energy measurement.", dict(provenance or {}))


def fingerprint_report(report: Mapping[str, Any]) -> str:
    payload = json.dumps(report, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_figure5_predictions(root: str | Path) -> tuple[dict[str, list[float]], list[float], dict[str, Any]]:
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("numpy is required to load the published benchmark") from exc
    directory = Path(root)
    names = {"PBA": "PBA_predictions.npy", "FateID": "FateID_predictions.npy", "WOT": "WOT_predictions.npy"}
    truth_path = directory / "smoothed_groundtruth_from_heldout.npy"
    mask_path = directory / "heldout_mask.npy"
    if not truth_path.exists() or not mask_path.exists():
        raise FileNotFoundError("Figure 5 held-out truth and mask are required")
    mask = np.asarray(np.load(mask_path), dtype=bool)
    truth = np.asarray(np.load(truth_path), dtype=float)[mask]
    predictions = {name: np.asarray(np.load(directory / filename), dtype=float)[mask].tolist() for name, filename in names.items()}
    provenance = {"root": str(directory), "heldout_cells": int(mask.sum()), "truth": truth_path.name, "mask": mask_path.name, "prediction_files": names}
    return predictions, truth.tolist(), provenance


def write_figure5_report(root: str | Path = "data/real/larry/figure5", output: str | Path = "experiments/larry_figure5_championship.json", *, bootstrap: int = 1000, permutations: int = 1000, seed: int = 17) -> dict[str, Any]:
    predictions, truth, provenance = load_figure5_predictions(root)
    report = evaluate_methods(predictions, truth, benchmark="Weinreb2020 LARRY Figure 5 held-out", outcome_definition="published smoothed ground-truth fate probability on the held-out mask", holdout_fraction=0.5, bootstrap=bootstrap, permutations=permutations, seed=seed, provenance=provenance)
    payload = report.to_dict()
    payload["fingerprint"] = fingerprint_report(payload)
    target = Path(output); target.parent.mkdir(parents=True, exist_ok=True); target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload
