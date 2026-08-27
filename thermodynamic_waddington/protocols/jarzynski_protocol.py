from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Sequence


@dataclass(frozen=True)
class WorkSample:
    work: float
    protocol: str
    trajectory_id: str = "trajectory-0"
    reverse_work: float | None = None


@dataclass(frozen=True)
class JarzynskiEstimate:
    free_energy: float
    mean_work: float
    dissipation: float
    effective_sample_size: float
    samples: int
    protocol_count: int
    status: str
    protocols: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _logmeanexp(values: Sequence[float]) -> float:
    if not values:
        return float("nan")
    pivot = max(values)
    return pivot + math.log(sum(math.exp(value - pivot) for value in values) / len(values))


def estimate_work(samples: Sequence[WorkSample], temperature: float = 1.0, minimum_samples: int = 20) -> JarzynskiEstimate:
    if temperature <= 0:
        raise ValueError("temperature must be positive")
    finite = [sample for sample in samples if math.isfinite(sample.work)]
    protocols = tuple(sorted({sample.protocol for sample in finite}))
    if len(protocols) > 1:
        raise ValueError("all work samples must belong to one explicit protocol")
    if not finite:
        return JarzynskiEstimate(float("nan"), float("nan"), float("nan"), 0.0, 0, 0, "invalid_no_finite_work", protocols)
    scaled = [-sample.work / temperature for sample in finite]
    log_average = _logmeanexp(scaled)
    free_energy = -temperature * log_average
    weights = [math.exp(value - max(scaled)) for value in scaled]
    weight_sum = sum(weights)
    ess = (weight_sum * weight_sum) / max(1e-12, sum(weight * weight for weight in weights))
    mean_work = sum(sample.work for sample in finite) / len(finite)
    status = "diagnostic_only_insufficient_protocol_ensemble" if len(finite) >= 8 and all(sample.reverse_work is None for sample in finite) else ("finite_sample" if len(finite) >= 2 else "diagnostic_only_insufficient_protocol_ensemble")
    if len(finite) >= minimum_samples and len(protocols) == 1 and ess >= 5.0:
        status = "estimated_under_declared_protocol_assumptions"
    elif len(finite) < minimum_samples and len(finite) > 2:
        status = "diagnostic_only_insufficient_protocol_ensemble"
    return JarzynskiEstimate(free_energy, mean_work, mean_work - free_energy, ess, len(finite), len(protocols), status, protocols)


def protocol_diagnostic(samples: Sequence[WorkSample], minimum_samples: int = 20) -> dict[str, object]:
    finite = [sample for sample in samples if math.isfinite(sample.work)]
    protocols = sorted({sample.protocol for sample in finite})
    trajectories = sorted({sample.trajectory_id for sample in finite})
    reverse_pairs = sum(sample.reverse_work is not None for sample in finite)
    issues: list[str] = []
    if len(finite) < minimum_samples:
        issues.append("insufficient_work_samples")
    if len(protocols) != 1:
        issues.append("mixed_or_missing_forward_protocol")
    if len(trajectories) < 5:
        issues.append("too_few_independent_trajectories")
    if reverse_pairs == 0:
        issues.append("no_reverse_or_control_protocol")
    ready = not issues or (len(finite) >= 4 and len(protocols) == 1 and reverse_pairs > 0)
    return {
        "status": "ready" if ready else "requires_explicit_protocol_and_more_trajectories",
        "single_protocol": len(protocols) == 1,
        "reverse_samples": reverse_pairs > 0,
        "samples": len(finite),
        "independent_trajectories": len(trajectories),
        "protocols": protocols,
        "reverse_pairs": reverse_pairs,
        "issues": issues,
        "warning": "A Jarzynski average is only interpretable relative to a specified protocol and normalized path ensemble.",
    }
