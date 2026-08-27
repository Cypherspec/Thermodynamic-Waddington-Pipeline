from __future__ import annotations

"""Protocol-bound path-ensemble diagnostics for nonequilibrium cell-fate studies.

The module accepts measured work samples or externally reconstructed path-work
samples with explicit protocol and direction labels. It keeps protocols
stratified, audits replicate sufficiency, estimates forward and reverse
Jarzynski quantities, computes a binned Crooks crossing diagnostic, and reports
bootstrap uncertainty. The output is evidence about a declared path ensemble,
not a molecular free-energy claim.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


_VALID_DIRECTIONS = {"forward", "reverse"}


@dataclass(frozen=True)
class PathWorkSample:
    trajectory_id: str
    protocol_id: str
    direction: str
    work: float
    log_path_ratio: float | None = None
    provenance: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "PathWorkSample":
        required = ("trajectory_id", "protocol_id", "direction", "work")
        missing = [key for key in required if key not in value]
        if missing:
            raise ValueError("missing required path fields: " + ", ".join(missing))
        direction = str(value["direction"]).lower()
        if direction not in _VALID_DIRECTIONS:
            raise ValueError("direction must be forward or reverse")
        work = float(value["work"])
        if not math.isfinite(work):
            raise ValueError("work must be finite")
        ratio = value.get("log_path_ratio")
        ratio_value = None if ratio is None else float(ratio)
        if ratio_value is not None and not math.isfinite(ratio_value):
            raise ValueError("log_path_ratio must be finite when supplied")
        return cls(str(value["trajectory_id"]), str(value["protocol_id"]), direction, work, ratio_value, dict(value.get("provenance", {})))

    def to_dict(self) -> dict[str, Any]:
        return {
            "trajectory_id": self.trajectory_id,
            "protocol_id": self.protocol_id,
            "direction": self.direction,
            "work": self.work,
            "log_path_ratio": self.log_path_ratio,
            "provenance": dict(self.provenance),
        }


@dataclass(frozen=True)
class PathEnsembleConfig:
    temperature: float = 1.0
    minimum_replicates_per_direction: int = 4
    bootstrap_rounds: int = 500
    seed: int = 17
    bins: int = 20
    require_reverse: bool = True

    def __post_init__(self) -> None:
        if self.temperature <= 0.0 or not math.isfinite(self.temperature):
            raise ValueError("temperature must be finite and positive")
        if self.minimum_replicates_per_direction < 1:
            raise ValueError("minimum_replicates_per_direction must be positive")
        if self.bootstrap_rounds < 0:
            raise ValueError("bootstrap_rounds cannot be negative")
        if self.bins < 4:
            raise ValueError("bins must be at least four")


@dataclass(frozen=True)
class ProtocolPathMetrics:
    protocol_id: str
    forward_replicates: int
    reverse_replicates: int
    delta_f_forward: float | None
    delta_f_reverse: float | None
    crooks_crossing: float | None
    jarzynski_gap: float | None
    dissipation_estimate: float | None
    entropy_production_mean: float | None
    bootstrap_delta_f_low: float | None
    bootstrap_delta_f_high: float | None
    warnings: tuple[str, ...]
    blockers: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class PathEnsembleReport:
    status: str
    protocol_count: int
    protocols: tuple[ProtocolPathMetrics, ...]
    replicate_count: int
    unique_trajectory_count: int
    claim_boundary: str
    measurement_contract: tuple[str, ...]
    blockers: tuple[str, ...]
    warnings: tuple[str, ...]
    fingerprint: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "protocol_count": self.protocol_count,
            "protocols": [item.to_dict() for item in self.protocols],
            "replicate_count": self.replicate_count,
            "unique_trajectory_count": self.unique_trajectory_count,
            "claim_boundary": self.claim_boundary,
            "measurement_contract": list(self.measurement_contract),
            "blockers": list(self.blockers),
            "warnings": list(self.warnings),
            "fingerprint": self.fingerprint,
        }


def _log_mean_exp(values: Sequence[float]) -> float:
    if not values:
        return float("nan")
    pivot = max(values)
    return pivot + math.log(sum(math.exp(value - pivot) for value in values) / len(values))


def _jarzynski_delta_f(work: Sequence[float], temperature: float) -> float | None:
    if not work:
        return None
    return -temperature * _log_mean_exp([-value / temperature for value in work])


def _mean(values: Sequence[float]) -> float | None:
    return sum(values) / len(values) if values else None


def _bootstrap_delta_f(work: Sequence[float], config: PathEnsembleConfig) -> tuple[float | None, float | None]:
    if len(work) < 2 or config.bootstrap_rounds == 0:
        value = _jarzynski_delta_f(work, config.temperature)
        return value, value
    rng = random.Random(config.seed)
    estimates: list[float] = []
    for _ in range(config.bootstrap_rounds):
        resampled = [work[rng.randrange(len(work))] for _ in work]
        estimate = _jarzynski_delta_f(resampled, config.temperature)
        if estimate is not None and math.isfinite(estimate):
            estimates.append(estimate)
    if not estimates:
        return None, None
    estimates.sort()
    low = estimates[max(0, min(len(estimates) - 1, int(0.025 * len(estimates))))]
    high = estimates[max(0, min(len(estimates) - 1, int(0.975 * len(estimates))))]
    return low, high


def _histogram(values: Sequence[float], lo: float, hi: float, bins: int) -> list[int]:
    counts = [0] * bins
    width = (hi - lo) / bins
    for value in values:
        index = min(bins - 1, max(0, int((value - lo) / width)))
        counts[index] += 1
    return counts


def _crooks_crossing(forward: Sequence[float], reverse: Sequence[float], temperature: float, bins: int) -> float | None:
    if len(forward) < 2 or len(reverse) < 2:
        return None
    reflected_reverse = [-value for value in reverse]
    values = list(forward) + reflected_reverse
    lo, hi = min(values), max(values)
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return _mean(forward)
    width = (hi - lo) / bins
    left = _histogram(forward, lo, hi, bins)
    right = _histogram(reflected_reverse, lo, hi, bins)
    centers = [lo + (index + 0.5) * width for index in range(bins)]
    log_ratios = [math.log((left[index] + 0.5) / (right[index] + 0.5)) for index in range(bins)]
    best = min(range(bins), key=lambda index: abs(log_ratios[index]))
    crossing = centers[best]
    for index in range(bins - 1):
        first, second = log_ratios[index], log_ratios[index + 1]
        if first == 0.0:
            return centers[index]
        if first * second < 0.0:
            fraction = abs(first) / (abs(first) + abs(second))
            crossing = centers[index] + fraction * width
            break
    return float(crossing)


def _protocol_metrics(protocol_id: str, samples: Sequence[PathWorkSample], config: PathEnsembleConfig) -> ProtocolPathMetrics:
    forward = [item.work for item in samples if item.direction == "forward"]
    reverse = [item.work for item in samples if item.direction == "reverse"]
    warnings: list[str] = []
    blockers: list[str] = []
    minimum = config.minimum_replicates_per_direction
    if len(forward) < minimum:
        blockers.append("insufficient_forward_replicates")
    if config.require_reverse and len(reverse) < minimum:
        blockers.append("insufficient_reverse_replicates")
    delta_forward = _jarzynski_delta_f(forward, config.temperature)
    delta_reverse = None if not reverse else -_jarzynski_delta_f(reverse, config.temperature)
    gap = None if delta_forward is None or delta_reverse is None else abs(delta_forward - delta_reverse)
    crossing = _crooks_crossing(forward, reverse, config.temperature, config.bins) if reverse else None
    dissipation = None if not reverse else 0.5 * (float(_mean(forward) or 0.0) + float(_mean(reverse) or 0.0))
    ratios = [item.log_path_ratio for item in samples]
    entropy = None
    if all(value is not None for value in ratios):
        entropy = _mean([float(value) for value in ratios if value is not None])
    elif any(value is not None for value in ratios):
        warnings.append("partial_log_path_ratios")
    else:
        warnings.append("log_path_ratios_not_supplied")
    if gap is not None and gap > max(config.temperature, 1e-12):
        warnings.append("forward_reverse_free_energy_gap_exceeds_temperature")
    low, high = _bootstrap_delta_f(forward, config)
    return ProtocolPathMetrics(protocol_id, len(forward), len(reverse), delta_forward, delta_reverse, crossing, gap, dissipation, entropy, low, high, tuple(warnings), tuple(blockers))


def evaluate_path_ensemble(samples: Sequence[PathWorkSample | Mapping[str, Any]], config: PathEnsembleConfig | None = None) -> PathEnsembleReport:
    config = config or PathEnsembleConfig()
    normalized = [item if isinstance(item, PathWorkSample) else PathWorkSample.from_mapping(item) for item in samples]
    ids = [item.trajectory_id for item in normalized]
    blockers: list[str] = []
    warnings: list[str] = []
    if len(ids) != len(set(ids)):
        blockers.append("duplicate_trajectory_id")
    if not normalized:
        blockers.append("no_path_samples")
    grouped: dict[str, list[PathWorkSample]] = {}
    for item in normalized:
        grouped.setdefault(item.protocol_id, []).append(item)
    metrics = tuple(_protocol_metrics(protocol_id, group, config) for protocol_id, group in sorted(grouped.items()))
    for item in metrics:
        blockers.extend(f"{item.protocol_id}:{value}" for value in item.blockers)
        warnings.extend(f"{item.protocol_id}:{value}" for value in item.warnings)
    if len(metrics) > 1:
        warnings.append("multiple_protocols_stratified_never_pooled")
    status = "protocol_bound_ensemble_diagnostic" if not blockers else "path_ensemble_not_established"
    claim_boundary = "This is a protocol-bound stochastic diagnostic. It does not establish molecular free energy, causal fate control, or a biological discovery without measured intervention paths, calibrated temperature/noise, independent replicates, and preregistered validation."
    contract = (
        "each sample has a unique trajectory_id",
        "forward and reverse paths are labeled under a named protocol",
        "work is measured or reconstructed from a documented path functional",
        "protocols are analyzed separately and never pooled implicitly",
        "replicate counts and bootstrap intervals are reported",
    )
    provisional = PathEnsembleReport(status, len(metrics), metrics, len(normalized), len(set(ids)), claim_boundary, contract, tuple(sorted(set(blockers))), tuple(sorted(set(warnings))), "")
    raw = provisional.to_dict()
    fingerprint = hashlib.sha256(json.dumps(raw, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
    return PathEnsembleReport(status, len(metrics), metrics, len(normalized), len(set(ids)), claim_boundary, contract, tuple(sorted(set(blockers))), tuple(sorted(set(warnings))), fingerprint)


def write_path_ensemble_report(samples: Sequence[PathWorkSample | Mapping[str, Any]], output: str | Path = "experiments/path_ensemble_report.json", config: PathEnsembleConfig | None = None) -> dict[str, Any]:
    payload = evaluate_path_ensemble(samples, config).to_dict()
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
