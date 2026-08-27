from __future__ import annotations

"""Adversarial, deterministic falsification checks for the effective landscape.

These checks are deliberately not biological validation. They probe whether the
estimator reacts in ways that would reveal leakage, coordinate dependence,
velocity sign errors, graph sparsity, or unsupported thermodynamic language.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

from .config import FitConfig
from .model import fit_landscape
from .synthetic import make_synthetic_dataset


@dataclass(frozen=True)
class FalsificationCase:
    name: str
    hypothesis: str
    passed: bool
    statistic: float
    threshold: float
    interpretation: str
    claim_boundary: str = "Software stress test only; not biological or thermodynamic evidence."

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _fit(bundle: dict[str, Any], seed: int) -> Any:
    config = FitConfig(
        neighbors=10,
        dimensions=min(4, len(bundle["expression"][0])),
        bootstrap_replicates=3,
        max_paths=30,
        seed=seed,
    )
    return fit_landscape(
        bundle["expression"],
        bundle.get("velocity"),
        embedding=bundle.get("embedding"),
        labels=bundle.get("labels"),
        cell_ids=bundle.get("cell_ids"),
        config=config,
        metadata={"falsification": True},
    )


def _energy_vector(fit: Any) -> list[float]:
    values = getattr(fit, "energies", None)
    if values is None:
        values = getattr(fit, "energy", None)
    return [float(x) for x in values or []]


def _rank_correlation(a: Sequence[float], b: Sequence[float]) -> float:
    n = min(len(a), len(b))
    if n < 2:
        return float("nan")
    def ranks(values: Sequence[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: (values[i], i))
        out = [0.0] * len(values)
        start = 0
        while start < len(order):
            end = start + 1
            while end < len(order) and values[order[end]] == values[order[start]]:
                end += 1
            rank = (start + end - 1) / 2.0 + 1.0
            for index in order[start:end]:
                out[index] = rank
            start = end
        return out
    ar, br = ranks(list(a)[:n]), ranks(list(b)[:n])
    am, bm = sum(ar) / n, sum(br) / n
    da = sum((x - am) ** 2 for x in ar)
    db = sum((y - bm) ** 2 for y in br)
    return sum((x - am) * (y - bm) for x, y in zip(ar, br)) / math.sqrt(da * db) if da and db else 0.0


def _clone_bundle(bundle: dict[str, Any]) -> dict[str, Any]:
    return json.loads(json.dumps(bundle))


def _base_bundle(seed: int) -> dict[str, Any]:
    result = make_synthetic_dataset(cells=96, genes=10, seed=seed)
    if hasattr(result, "to_dict"):
        result = result.to_dict()
    return result


def _case_seed_determinism(bundle: dict[str, Any]) -> FalsificationCase:
    first, second = _fit(bundle, 41), _fit(bundle, 41)
    a, b = _energy_vector(first), _energy_vector(second)
    statistic = max((abs(x - y) for x, y in zip(a, b)), default=float("inf"))
    return FalsificationCase("seed_determinism", "same inputs and seed produce the same energy vector", statistic <= 1e-12, statistic, 1e-12, "A failure indicates nondeterministic numerics or unstable serialization.")


def _case_velocity_sign(bundle: dict[str, Any]) -> FalsificationCase:
    normal = _fit(bundle, 43)
    reversed_bundle = _clone_bundle(bundle)
    reversed_bundle["velocity"] = [[-float(v) for v in row] for row in reversed_bundle.get("velocity", [])]
    reversed_fit = _fit(reversed_bundle, 43)
    a, b = _energy_vector(normal), _energy_vector(reversed_fit)
    statistic = _rank_correlation(a, b)
    passed = _finite(statistic) and statistic < 0.999999
    return FalsificationCase("velocity_sign_sensitivity", "reversing velocity must not be silently identical", passed, statistic if _finite(statistic) else 0.0, 0.999999, "A perfect identity can reveal ignored velocity inputs or a zero-velocity placeholder.")


def _case_permutation(bundle: dict[str, Any]) -> FalsificationCase:
    permuted = _clone_bundle(bundle)
    rng = random.Random(53)
    rng.shuffle(permuted["velocity"])
    first, second = _fit(bundle, 47), _fit(permuted, 47)
    statistic = _rank_correlation(_energy_vector(first), _energy_vector(second))
    passed = _finite(statistic) and statistic < 0.999999
    return FalsificationCase("velocity_permutation_sensitivity", "cell-velocity pairing carries information", passed, statistic if _finite(statistic) else 0.0, 0.999999, "A failure suggests the fit does not use cell-specific velocity pairing.")


def _case_coordinate_invariance(bundle: dict[str, Any]) -> FalsificationCase:
    transformed = _clone_bundle(bundle)
    transformed["embedding"] = [[3.0 * float(row[0]) + 7.0, -2.0 * float(row[1]) - 4.0] for row in transformed["embedding"]]
    first, second = _fit(bundle, 59), _fit(transformed, 59)
    statistic = _rank_correlation(_energy_vector(first), _energy_vector(second))
    passed = _finite(statistic) and statistic > 0.95
    return FalsificationCase("affine_coordinate_stability", "affine embedding units should preserve energy ordering", passed, statistic if _finite(statistic) else 0.0, 0.95, "A failure indicates scale-sensitive geometry or an undocumented coordinate convention.")


def _case_finite_output(bundle: dict[str, Any]) -> FalsificationCase:
    fit = _fit(bundle, 61)
    values = _energy_vector(fit)
    bad = sum(1 for value in values if not _finite(value))
    statistic = float(bad)
    return FalsificationCase("finite_output", "ordinary finite inputs produce finite energies", bad == 0, statistic, 0.0, "Non-finite output is a hard numerical failure, not a biological result.")


def _case_null_velocity(bundle: dict[str, Any]) -> FalsificationCase:
    null_bundle = _clone_bundle(bundle)
    null_bundle["velocity"] = [[0.0 for _ in row] for row in null_bundle["velocity"]]
    fit = _fit(null_bundle, 67)
    audit = getattr(fit, "diagnostics", {}) or {}
    serialized = json.dumps(audit, sort_keys=True).lower()
    exposed = any(token in serialized for token in ("unobserved", "placeholder", "zero velocity", "velocity_observed"))
    return FalsificationCase("null_velocity_disclosure", "zero velocity input is disclosed as a limitation", exposed, 1.0 if exposed else 0.0, 1.0, "The model must not turn a missing or null velocity field into evidence of thermodynamic equilibrium.")


def run_falsification_suite(seed: int = 17) -> dict[str, Any]:
    bundle = _base_bundle(seed)
    cases = [
        _case_seed_determinism(bundle),
        _case_velocity_sign(bundle),
        _case_permutation(bundle),
        _case_coordinate_invariance(bundle),
        _case_finite_output(bundle),
        _case_null_velocity(bundle),
    ]
    payload: dict[str, Any] = {
        "suite": "TW-adversarial-falsification-v1",
        "seed": seed,
        "cases": [case.to_dict() for case in cases],
        "passed": sum(case.passed for case in cases),
        "total": len(cases),
        "status": "passed" if all(case.passed for case in cases) else "failed",
        "claim_boundary": "Passing these checks validates implementation behavior under synthetic perturbations only. It does not establish biological causality, molecular free energy, or wet-lab truth.",
    }
    payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return payload


def write_falsification_report(output: str | Path = "experiments/falsification_report.json", seed: int = 17) -> dict[str, Any]:
    payload = run_falsification_suite(seed)
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
