from __future__ import annotations

"""Audited cross-study transportability for cell-fate predictors.

This module evaluates whether a fixed predictor behaves consistently across
independent datasets without converting predictive agreement into causal or
thermodynamic evidence. It requires explicit study provenance and keeps
cell-weighted and equal-study summaries separate.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class TransportabilityConfig:
    min_studies: int = 2
    min_cells_per_study: int = 20
    bootstrap_rounds: int = 500
    seed: int = 17
    max_i2: float = 0.75
    min_positive_sign_consistency: float = 0.75
    max_negative_control_abs_spearman: float = 0.20


@dataclass(frozen=True)
class StudyObservation:
    study_id: str
    predictions: tuple[float, ...]
    outcomes: tuple[float, ...]
    domain: str = "unspecified"
    provenance: Mapping[str, Any] = field(default_factory=dict)
    negative_control: bool = False

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "StudyObservation":
        required = ("study_id", "predictions", "outcomes")
        missing = [key for key in required if key not in value]
        if missing:
            raise ValueError(f"study is missing required fields: {', '.join(missing)}")
        return cls(
            study_id=str(value["study_id"]),
            predictions=tuple(float(item) for item in value["predictions"]),
            outcomes=tuple(float(item) for item in value["outcomes"]),
            domain=str(value.get("domain", "unspecified")),
            provenance=dict(value.get("provenance", {})),
            negative_control=bool(value.get("negative_control", False)),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "study_id": self.study_id,
            "predictions": list(self.predictions),
            "outcomes": list(self.outcomes),
            "domain": self.domain,
            "provenance": dict(self.provenance),
            "negative_control": self.negative_control,
        }


@dataclass(frozen=True)
class StudyMetric:
    study_id: str
    domain: str
    n: int
    spearman: float
    pearson: float
    top_decile_recall: float
    bootstrap_low: float
    bootstrap_high: float
    finite: bool
    negative_control: bool
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class TransportabilityReport:
    status: str
    metrics: tuple[StudyMetric, ...]
    pooled_equal_study: Mapping[str, float]
    pooled_cell_weighted: Mapping[str, float]
    heterogeneity: Mapping[str, float]
    negative_controls: Mapping[str, Any]
    blockers: tuple[str, ...]
    warnings: tuple[str, ...]
    claim_boundary: str
    provenance: Mapping[str, Any]
    fingerprint: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "metrics": [metric.to_dict() for metric in self.metrics],
            "pooled_equal_study": dict(self.pooled_equal_study),
            "pooled_cell_weighted": dict(self.pooled_cell_weighted),
            "heterogeneity": dict(self.heterogeneity),
            "negative_controls": dict(self.negative_controls),
            "blockers": list(self.blockers),
            "warnings": list(self.warnings),
            "claim_boundary": self.claim_boundary,
            "provenance": dict(self.provenance),
            "fingerprint": self.fingerprint,
        }


def _finite_pairs(predictions: Sequence[float], outcomes: Sequence[float]) -> tuple[list[float], list[float]]:
    pairs: list[tuple[float, float]] = []
    for prediction, outcome in zip(predictions, outcomes):
        try:
            left, right = float(prediction), float(outcome)
        except (TypeError, ValueError):
            continue
        if math.isfinite(left) and math.isfinite(right):
            pairs.append((left, right))
    return [left for left, _ in pairs], [right for _, right in pairs]


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def _rank(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    ranks = [0.0] * len(values)
    position = 0
    while position < len(order):
        end = position + 1
        while end < len(order) and values[order[end]] == values[order[position]]:
            end += 1
        average = 0.5 * (position + end - 1) + 1.0
        for cursor in range(position, end):
            ranks[order[cursor]] = average
        position = end
    return ranks


def _pearson(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) < 2:
        return float("nan")
    left_mean, right_mean = _mean(left), _mean(right)
    centered_left = [value - left_mean for value in left]
    centered_right = [value - right_mean for value in right]
    denominator = math.sqrt(sum(value * value for value in centered_left) * sum(value * value for value in centered_right))
    return sum(a * b for a, b in zip(centered_left, centered_right)) / denominator if denominator > 0.0 else 0.0


def _spearman(left: Sequence[float], right: Sequence[float]) -> float:
    return _pearson(_rank(left), _rank(right))


def _top_decile_recall(predictions: Sequence[float], outcomes: Sequence[float]) -> float:
    if not predictions:
        return float("nan")
    count = max(1, math.ceil(len(predictions) * 0.10))
    predicted = set(sorted(range(len(predictions)), key=lambda index: predictions[index], reverse=True)[:count])
    observed = set(sorted(range(len(outcomes)), key=lambda index: outcomes[index], reverse=True)[:count])
    return len(predicted & observed) / len(observed) if observed else float("nan")


def _bootstrap_interval(predictions: Sequence[float], outcomes: Sequence[float], rounds: int, seed: int) -> tuple[float, float]:
    if len(predictions) < 2:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    estimates: list[float] = []
    indices = list(range(len(predictions)))
    for _ in range(max(100, rounds)):
        sample = [indices[rng.randrange(len(indices))] for _ in indices]
        estimates.append(_spearman([predictions[index] for index in sample], [outcomes[index] for index in sample]))
    estimates.sort()
    low_index = max(0, min(len(estimates) - 1, int(0.025 * len(estimates))))
    high_index = max(0, min(len(estimates) - 1, int(0.975 * len(estimates))))
    return estimates[low_index], estimates[high_index]


def _study_metric(study: StudyObservation, config: TransportabilityConfig, seed_offset: int) -> StudyMetric:
    predictions, outcomes = _finite_pairs(study.predictions, study.outcomes)
    warnings: list[str] = []
    if len(study.predictions) != len(study.outcomes):
        warnings.append("prediction_outcome_lengths_differ")
    if len(predictions) != len(study.predictions):
        warnings.append("nonfinite_pairs_dropped")
    if len(predictions) < config.min_cells_per_study:
        warnings.append("below_minimum_cells_per_study")
    low, high = _bootstrap_interval(predictions, outcomes, config.bootstrap_rounds, config.seed + seed_offset)
    return StudyMetric(
        study_id=study.study_id,
        domain=study.domain,
        n=len(predictions),
        spearman=_spearman(predictions, outcomes),
        pearson=_pearson(predictions, outcomes),
        top_decile_recall=_top_decile_recall(predictions, outcomes),
        bootstrap_low=low,
        bootstrap_high=high,
        finite=bool(predictions) and all(math.isfinite(value) for value in predictions + outcomes),
        negative_control=study.negative_control,
        warnings=tuple(warnings),
    )


def _aggregate(metrics: Sequence[StudyMetric], equal_study: bool) -> dict[str, float]:
    eligible = [metric for metric in metrics if metric.n > 0 and math.isfinite(metric.spearman)]
    if not eligible:
        return {"studies": 0.0, "cells": 0.0, "spearman": float("nan"), "pearson": float("nan"), "top_decile_recall": float("nan")}
    if equal_study:
        weights = [1.0 for _ in eligible]
    else:
        weights = [float(metric.n) for metric in eligible]
    total = sum(weights)
    return {
        "studies": float(len(eligible)),
        "cells": float(sum(metric.n for metric in eligible)),
        "spearman": sum(metric.spearman * weight for metric, weight in zip(eligible, weights)) / total,
        "pearson": sum(metric.pearson * weight for metric, weight in zip(eligible, weights)) / total,
        "top_decile_recall": sum(metric.top_decile_recall * weight for metric, weight in zip(eligible, weights)) / total,
    }


def _heterogeneity(metrics: Sequence[StudyMetric]) -> dict[str, float]:
    eligible = [metric for metric in metrics if metric.n >= 4 and math.isfinite(metric.spearman)]
    if len(eligible) < 2:
        return {"studies": float(len(eligible)), "q": float("nan"), "i2": float("nan"), "between_study_variance": float("nan"), "sign_consistency": float("nan")}
    transformed: list[float] = []
    variances: list[float] = []
    for metric in eligible:
        clipped = min(0.999999, max(-0.999999, metric.spearman))
        transformed.append(math.atanh(clipped))
        variances.append(1.0 / max(1.0, metric.n - 3.0))
    weights = [1.0 / variance for variance in variances]
    weighted_mean = sum(value * weight for value, weight in zip(transformed, weights)) / sum(weights)
    q = sum(weight * (value - weighted_mean) ** 2 for value, weight in zip(transformed, weights))
    degrees = len(eligible) - 1
    i2 = max(0.0, (q - degrees) / q) if q > 0.0 else 0.0
    between = max(0.0, (q - degrees) / max(1e-12, sum(weights) - sum(weight * weight for weight in weights) / sum(weights)))
    signs = [metric.spearman > 0.0 for metric in eligible]
    return {"studies": float(len(eligible)), "q": q, "i2": i2, "between_study_variance": between, "sign_consistency": sum(signs) / len(signs)}


def _study_signature(study: StudyObservation) -> str:
    payload = json.dumps({"study_id": study.study_id, "domain": study.domain, "provenance": dict(study.provenance), "negative_control": study.negative_control}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def evaluate_transportability(studies: Sequence[StudyObservation | Mapping[str, Any]], config: TransportabilityConfig | None = None) -> TransportabilityReport:
    config = config or TransportabilityConfig()
    normalized = [study if isinstance(study, StudyObservation) else StudyObservation.from_mapping(study) for study in studies]
    blockers: list[str] = []
    warnings: list[str] = []
    if len(normalized) < config.min_studies:
        blockers.append("insufficient_independent_studies")
    identifiers = [study.study_id for study in normalized]
    signatures = [_study_signature(study) for study in normalized]
    if len(identifiers) != len(set(identifiers)):
        blockers.append("duplicate_study_id")
    if len(signatures) != len(set(signatures)):
        blockers.append("duplicate_study_signature")
    if any(not study.provenance for study in normalized):
        warnings.append("one_or_more_studies_missing_provenance")
    metrics = tuple(_study_metric(study, config, index * 1009) for index, study in enumerate(normalized))
    if any(metric.n < config.min_cells_per_study for metric in metrics):
        blockers.append("study_below_minimum_finite_cells")
    if any(not metric.finite for metric in metrics):
        blockers.append("nonfinite_study_metric")
    heterogeneity = _heterogeneity(metrics)
    if math.isfinite(heterogeneity["i2"]) and heterogeneity["i2"] > config.max_i2:
        blockers.append("excessive_between_study_heterogeneity")
    if math.isfinite(heterogeneity["sign_consistency"]) and heterogeneity["sign_consistency"] < config.min_positive_sign_consistency:
        blockers.append("inconsistent_effect_direction_across_studies")
    controls = [metric for metric in metrics if metric.negative_control]
    control_failures = [metric.study_id for metric in controls if math.isfinite(metric.spearman) and abs(metric.spearman) > config.max_negative_control_abs_spearman]
    if control_failures:
        blockers.append("negative_control_signal_detected")
    negative_controls = {
        "n": len(controls),
        "study_ids": [metric.study_id for metric in controls],
        "failures": control_failures,
        "threshold_abs_spearman": config.max_negative_control_abs_spearman,
        "status": "passed" if not control_failures else "failed",
    }
    status = "transportable_predictive_signal_under_audited_benchmark" if not blockers else "transportability_not_established"
    provenance = {
        "module": "thermodynamic_waddington.transportability",
        "study_ids": identifiers,
        "study_signatures": signatures,
        "domains": [study.domain for study in normalized],
        "config": asdict(config),
        "predictor_is_fixed": True,
        "analysis_is_causal": False,
    }
    payload = {
        "status": status,
        "metrics": [metric.to_dict() for metric in metrics],
        "pooled_equal_study": _aggregate(metrics, True),
        "pooled_cell_weighted": _aggregate(metrics, False),
        "heterogeneity": heterogeneity,
        "negative_controls": negative_controls,
        "blockers": blockers,
        "warnings": warnings,
        "claim_boundary": "Cross-study predictive consistency is not causal evidence, mechanistic evidence, or a molecular free-energy measurement.",
        "provenance": provenance,
    }
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=True).encode("utf-8")).hexdigest()
    return TransportabilityReport(
        status=status,
        metrics=metrics,
        pooled_equal_study=_aggregate(metrics, True),
        pooled_cell_weighted=_aggregate(metrics, False),
        heterogeneity=heterogeneity,
        negative_controls=negative_controls,
        blockers=tuple(blockers),
        warnings=tuple(warnings),
        claim_boundary=payload["claim_boundary"],
        provenance=provenance,
        fingerprint=fingerprint,
    )


def write_transportability_report(studies: Sequence[StudyObservation | Mapping[str, Any]], output: str | Path, config: TransportabilityConfig | None = None) -> dict[str, Any]:
    report = evaluate_transportability(studies, config)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = report.to_dict()
    target.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n", encoding="utf-8")
    return payload


def load_transportability_manifest(path: str | Path) -> tuple[list[StudyObservation], TransportabilityConfig]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, Mapping) or not isinstance(payload.get("studies"), list):
        raise ValueError("transportability manifest must contain a studies list")
    config_payload = payload.get("config", {})
    if not isinstance(config_payload, Mapping):
        raise ValueError("transportability manifest config must be an object")
    allowed = {field_name for field_name in TransportabilityConfig.__dataclass_fields__}
    config = TransportabilityConfig(**{key: value for key, value in config_payload.items() if key in allowed})
    return [StudyObservation.from_mapping(item) for item in payload["studies"]], config
