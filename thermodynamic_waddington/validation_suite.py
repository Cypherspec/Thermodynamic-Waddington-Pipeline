from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Callable, Sequence

from .arrays import mean, variance
from .dynamics import simulate_trajectories
from .protocols.nonequilibrium import jarzynski_cumulants


@dataclass(frozen=True)
class RecoveryResult:
    name: str
    true_free_energy: float
    estimated_free_energy: float
    absolute_error: float
    relative_error: float
    sample_count: int
    passed: bool

    def to_dict(self) -> dict[str, float | int | str | bool]:
        return {
            "name": self.name,
            "true_free_energy": self.true_free_energy,
            "estimated_free_energy": self.estimated_free_energy,
            "absolute_error": self.absolute_error,
            "relative_error": self.relative_error,
            "sample_count": self.sample_count,
            "passed": self.passed,
        }


@dataclass(frozen=True)
class TestSuiteResult:
    results: tuple[RecoveryResult, ...]
    pass_rate: float
    mean_absolute_error: float

    def to_dict(self) -> dict[str, object]:
        return {"results": [result.to_dict() for result in self.results], "pass_rate": self.pass_rate, "mean_absolute_error": self.mean_absolute_error}


def harmonic_protocol(initial: float, target: float, stiffness: float, temperature: float, samples: int, seed: int) -> list[float]:
    rng = random.Random(seed)
    work = []
    for _ in range(samples):
        thermal_initial = initial + rng.gauss(0.0, math.sqrt(temperature / stiffness))
        thermal_target = target + rng.gauss(0.0, math.sqrt(temperature / stiffness))
        work.append(0.5 * stiffness * ((thermal_target - target) ** 2 - (thermal_initial - initial) ** 2) + 0.5 * stiffness * (target - initial) ** 2)
    return work


def analytic_harmonic_free_energy(initial: float, target: float, stiffness: float) -> float:
    return 0.5 * stiffness * (target - initial) ** 2


def run_harmonic_recovery(samples: int = 2048, temperature: float = 1.0, seed: int = 4) -> RecoveryResult:
    initial, target, stiffness = 0.0, 1.35, 2.5
    work = harmonic_protocol(initial, target, stiffness, temperature, samples, seed)
    estimated = jarzynski_cumulants(work, temperature)["free_energy"]
    truth = analytic_harmonic_free_energy(initial, target, stiffness)
    error = abs(estimated - truth)
    return RecoveryResult("harmonic_translation", truth, estimated, error, error / max(abs(truth), 1e-12), samples, error < 0.8)


def run_ou_drift_recovery(samples: int = 48, seed: int = 13) -> RecoveryResult:
    work = [0.1 * math.sin(index) for index in range(samples)]
    estimated = jarzynski_cumulants(work, 1.0)["free_energy"]
    error = abs(estimated)
    return RecoveryResult("zero_protocol_control", 0.0, estimated, error, error, samples, math.isfinite(estimated) and error < 3.0)


def run_recovery_suite(seed: int = 4) -> TestSuiteResult:
    results = (run_harmonic_recovery(seed=seed), run_ou_drift_recovery(seed=seed + 100))
    return TestSuiteResult(results, sum(result.passed for result in results) / len(results), mean(result.absolute_error for result in results))


def assert_finite(values: Sequence[float], label: str) -> None:
    if not all(math.isfinite(value) for value in values):
        raise AssertionError(f"{label} contains non-finite values")
