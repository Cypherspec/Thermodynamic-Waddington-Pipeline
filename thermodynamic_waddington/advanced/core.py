from __future__ import annotations

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from typing import Any, Iterable, Mapping, Sequence

Number = int | float
Row = Sequence[Number]
Matrix = Sequence[Row]


def finite(value: Number, default: float = 0.0) -> float:
    value = float(value)
    return value if math.isfinite(value) else default


def vector(values: Any) -> list[float]:
    if values is None:
        return []
    if isinstance(values, (int, float)):
        return [finite(values)]
    return [finite(value) for value in values]


def matrix(values: Any) -> list[list[float]]:
    if values is None:
        return []
    if isinstance(values, (int, float)):
        return [[finite(values)]]
    rows = list(values)
    if not rows:
        return []
    if isinstance(rows[0], (int, float)):
        return [vector(rows)]
    return [vector(row) for row in rows]


def flatten(values: Any) -> list[float]:
    return [value for row in matrix(values) for value in row]


def clamp(value: Number, low: Number, high: Number) -> float:
    return max(float(low), min(float(high), finite(value)))


def mean(values: Iterable[Number]) -> float:
    items = [finite(value) for value in values]
    return sum(items) / len(items) if items else 0.0


def variance(values: Iterable[Number]) -> float:
    items = [finite(value) for value in values]
    if len(items) < 2:
        return 0.0
    center = mean(items)
    return sum((value - center) ** 2 for value in items) / (len(items) - 1)


def quantile(values: Iterable[Number], probability: Number) -> float:
    items = sorted(finite(value) for value in values)
    if not items:
        return 0.0
    position = clamp(probability, 0.0, 1.0) * (len(items) - 1)
    low = int(math.floor(position))
    high = int(math.ceil(position))
    return items[low] + (items[high] - items[low]) * (position - low)


def dot(left: Row, right: Row) -> float:
    return sum(finite(a) * finite(b) for a, b in zip(left, right))


def norm(row: Row) -> float:
    return math.sqrt(max(0.0, dot(row, row)))


def distance(left: Row, right: Row) -> float:
    width = min(len(left), len(right))
    return math.sqrt(sum((finite(left[i]) - finite(right[i])) ** 2 for i in range(width)))


def normalize(row: Row) -> list[float]:
    scale = norm(row)
    return [finite(value) / scale for value in row] if scale > 1e-12 else [0.0 for _ in row]


def softmax(values: Row, temperature: Number = 1.0) -> list[float]:
    items = vector(values)
    if not items:
        return []
    temperature = max(1e-12, finite(temperature, 1.0))
    pivot = max(items)
    weights = [math.exp((value - pivot) / temperature) for value in items]
    total = sum(weights)
    return [weight / total for weight in weights] if total else [1.0 / len(items)] * len(items)


def logsumexp(values: Row) -> float:
    items = vector(values)
    if not items:
        return float('-inf')
    pivot = max(items)
    return pivot + math.log(sum(math.exp(value - pivot) for value in items))


def entropy(values: Row) -> float:
    return -sum(probability * math.log(max(probability, 1e-300)) for probability in softmax(values))


def covariance(left: Matrix, right: Matrix | None = None) -> list[list[float]]:
    a = matrix(left)
    b = matrix(right if right is not None else left)
    width = min(max((len(row) for row in a), default=0), max((len(row) for row in b), default=0))
    if not a or not b or not width:
        return []
    means_a = [mean(row[j] for row in a if j < len(row)) for j in range(width)]
    means_b = [mean(row[j] for row in b if j < len(row)) for j in range(width)]
    n = max(1, min(len(a), len(b)) - 1)
    return [[sum((a[i][j] - means_a[j]) * (b[i][k] - means_b[k]) for i in range(min(len(a), len(b))) if j < len(a[i]) and k < len(b[i])) / n for k in range(width)] for j in range(width)]


def nearest(values: Matrix, query: Row, limit: int = 8) -> list[tuple[int, float]]:
    scored = [(index, distance(row, query)) for index, row in enumerate(matrix(values))]
    return sorted(scored, key=lambda item: (item[1], item[0]))[:max(0, int(limit))]


def hash_payload(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(',', ':'), default=str)
    return hashlib.sha256(payload.encode('utf-8')).hexdigest()


def deterministic_rng(seed: int = 17) -> random.Random:
    return random.Random(int(seed))


@dataclass(frozen=True)
class ObservationContract:
    cells: int
    features: int
    velocity_cells: int
    velocity_features: int
    finite_fraction: float
    aligned: bool
    warnings: tuple[str, ...] = ()

    @classmethod
    def inspect(cls, expression: Any, velocity: Any) -> 'ObservationContract':
        x = matrix(expression)
        v = matrix(velocity)
        x_count = sum(len(row) for row in x)
        v_count = sum(len(row) for row in v)
        finite_count = sum(math.isfinite(value) for value in flatten(x)) + sum(math.isfinite(value) for value in flatten(v))
        total = x_count + v_count
        warnings: list[str] = []
        aligned = bool(x) and len(x) == len(v) and max((len(row) for row in x), default=0) == max((len(row) for row in v), default=0)
        if not aligned:
            warnings.append('expression_velocity_shape_mismatch')
        if not x:
            warnings.append('empty_expression')
        if not v:
            warnings.append('empty_velocity')
        return cls(len(x), max((len(row) for row in x), default=0), len(v), max((len(row) for row in v), default=0), finite_count / max(1, total), aligned, tuple(warnings))

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class EvidenceLedger:
    seed: int = 17
    records: list[dict[str, Any]] = field(default_factory=list)

    def add(self, stage: str, value: Any, status: str = 'diagnostic', **metadata: Any) -> dict[str, Any]:
        record = {'stage': stage, 'status': status, 'digest': hash_payload(value), 'metadata': metadata}
        self.records.append(record)
        return record

    def summary(self) -> dict[str, Any]:
        return {'seed': self.seed, 'stages': len(self.records), 'digest': hash_payload(self.records), 'records': list(self.records)}


@dataclass
class OnlineMoments:
    count: int = 0
    mean_value: float = 0.0
    m2: float = 0.0

    def update(self, value: Number) -> None:
        value = finite(value)
        self.count += 1
        delta = value - self.mean_value
        self.mean_value += delta / self.count
        self.m2 += delta * (value - self.mean_value)

    @property
    def variance_value(self) -> float:
        return self.m2 / (self.count - 1) if self.count > 1 else 0.0

    def to_dict(self) -> dict[str, Any]:
        return {'count': self.count, 'mean': self.mean_value, 'variance': self.variance_value}


class ReproducibleEngine:
    def __init__(self, seed: int = 17):
        self.seed = int(seed)
        self.rng = deterministic_rng(seed)
        self.ledger = EvidenceLedger(seed)

    def evaluate(self, rows: Any, scale: Number = 1.0) -> dict[str, Any]:
        values = matrix(rows)
        scores = [0.5 * norm(row) ** 2 / max(1e-12, finite(scale, 1.0)) for row in values]
        return {'scores': scores, 'mean': mean(scores), 'uncertainty': math.sqrt(variance(scores)), 'count': len(scores)}

    def checkpoint(self, stage: str, value: Any, **metadata: Any) -> dict[str, Any]:
        return self.ledger.add(stage, value, **metadata)

    def snapshot(self) -> dict[str, Any]:
        return {'seed': self.seed, 'ledger': self.ledger.summary()}
