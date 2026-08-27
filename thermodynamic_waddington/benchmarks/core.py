from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass
from typing import Callable, Sequence

from ..arrays import mean, norm


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    cells: int
    features: int
    neighbors: int
    elapsed_seconds: float
    peak_working_set_proxy: int
    finite_output: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "cells": self.cells,
            "features": self.features,
            "neighbors": self.neighbors,
            "elapsed_seconds": self.elapsed_seconds,
            "peak_working_set_proxy": self.peak_working_set_proxy,
            "finite_output": self.finite_output,
        }


def random_matrix(rows: int, columns: int, seed: int) -> list[list[float]]:
    rng = random.Random(seed)
    return [[rng.gauss(0.0, 1.0) for _ in range(columns)] for _ in range(rows)]


def benchmark(name: str, operation: Callable[[Sequence[Sequence[float]]], Sequence[float]], cells: int, features: int, neighbors: int, seed: int = 7) -> BenchmarkCase:
    matrix = random_matrix(cells, features, seed)
    started = time.perf_counter()
    output = list(operation(matrix))
    elapsed = time.perf_counter() - started
    return BenchmarkCase(name, cells, features, neighbors, elapsed, cells * features * 8, all(math.isfinite(value) for value in output))


def scaling_summary(cases: Sequence[BenchmarkCase]) -> dict[str, object]:
    if not cases:
        return {"cases": [], "cells_per_second": 0.0}
    throughput = mean(case.cells / max(case.elapsed_seconds, 1e-12) for case in cases)
    return {"cases": [case.to_dict() for case in cases], "cells_per_second": throughput, "worst_case_seconds": max(case.elapsed_seconds for case in cases)}
