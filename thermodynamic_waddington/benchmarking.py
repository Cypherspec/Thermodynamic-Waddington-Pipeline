from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .arrays import mean


@dataclass(frozen=True)
class BenchmarkResult:
    method: str
    n: int
    mean_absolute_error: float
    rank_correlation: float
    basin_recall: float
    notes: str

    def to_dict(self) -> dict[str, float | int | str]:
        return {"method": self.method, "n": self.n, "mean_absolute_error": self.mean_absolute_error, "rank_correlation": self.rank_correlation, "basin_recall": self.basin_recall, "notes": self.notes}


def rank_correlation(first: Sequence[float], second: Sequence[float]) -> float:
    if len(first) != len(second) or not first:
        return 0.0
    left = sorted(range(len(first)), key=lambda i: first[i])
    right = sorted(range(len(second)), key=lambda i: second[i])
    position = {node: index for index, node in enumerate(right)}
    d2 = sum((index - position[node]) ** 2 for index, node in enumerate(left))
    n = len(first)
    return 1.0 - 6.0 * d2 / max(1, n * (n * n - 1))


def compare_energy_methods(reference: Sequence[float], candidates: dict[str, Sequence[float]], reference_basins: set[int] | None = None) -> list[BenchmarkResult]:
    output: list[BenchmarkResult] = []
    reference_basins = reference_basins or set()
    for method, values in candidates.items():
        error = mean([abs(a - b) for a, b in zip(reference, values)]) if len(reference) == len(values) else float("inf")
        predicted = set(sorted(range(len(values)), key=lambda i: values[i])[: max(1, len(reference_basins))])
        recall = len(predicted & reference_basins) / max(1, len(reference_basins))
        output.append(BenchmarkResult(method, len(values), error, rank_correlation(reference, values), recall, "diagnostic comparison; not a claim of biological superiority"))
    return sorted(output, key=lambda item: item.mean_absolute_error)
