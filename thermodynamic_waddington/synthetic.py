from __future__ import annotations

import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from .arrays import norm


@dataclass
class SyntheticDataset:
    expression: list[list[float]]
    velocity: list[list[float]]
    embedding: list[list[float]]
    labels: list[str]
    lineage_outcomes: list[str]
    cell_ids: list[str]
    metadata: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), separators=(",", ":")))

    @classmethod
    def load(cls, path: str | Path) -> "SyntheticDataset":
        return cls(**json.loads(Path(path).read_text()))


def _potential(x: float, y: float) -> float:
    wells = [(-2.2, -1.3, 0.35), (-1.0, 1.8, 0.30), (1.7, 0.9, 0.22), (2.0, -1.5, 0.18)]
    baseline = 0.08 * (x * x + y * y) + 0.03 * x * y
    return baseline - sum(depth * math.exp(-((x - cx) ** 2 + (y - cy) ** 2) / 0.75) for cx, cy, depth in wells)


def _gradient(x: float, y: float) -> tuple[float, float]:
    epsilon = 1e-4
    dx = (_potential(x + epsilon, y) - _potential(x - epsilon, y)) / (2 * epsilon)
    dy = (_potential(x, y + epsilon) - _potential(x, y - epsilon)) / (2 * epsilon)
    return dx, dy


def _gene_response(x: float, y: float, gene: int, rng: random.Random) -> float:
    phase = gene * 0.71
    nonlinear = math.sin((gene + 1) * x * 0.55 + phase) + math.cos((gene + 2) * y * 0.45 - phase)
    interaction = 0.16 * x * y if gene % 3 == 0 else 0.05 * (x * x - y * y)
    return nonlinear + interaction + rng.gauss(0, 0.035)


def make_synthetic_dataset(cells: int = 800, genes: int = 16, seed: int = 17) -> SyntheticDataset:
    if cells < 20 or genes < 4:
        raise ValueError("synthetic benchmark needs at least 20 cells and 4 genes")
    rng = random.Random(seed)
    embedding: list[list[float]] = []
    expression: list[list[float]] = []
    velocity: list[list[float]] = []
    labels: list[str] = []
    outcomes: list[str] = []
    wells = [("HSC", -2.2, -1.3), ("erythroid", -1.0, 1.8), ("myeloid", 1.7, 0.9), ("megakaryocyte", 2.0, -1.5)]
    for index in range(cells):
        well_index = 0 if index < cells // 4 else rng.randrange(len(wells))
        label, cx, cy = wells[well_index]
        radius = abs(rng.gauss(0.0, 0.42 if label != "HSC" else 0.58))
        angle = rng.random() * 2 * math.pi
        x = cx + radius * math.cos(angle)
        y = cy + radius * math.sin(angle)
        gx, gy = _gradient(x, y)
        noise = 0.12 + 0.03 * (abs(x) + abs(y))
        embedding.append([x, y])
        expression.append([_gene_response(x, y, gene, rng) for gene in range(genes)])
        velocity.append([-gx + rng.gauss(0, noise), -gy + rng.gauss(0, noise)] + [rng.gauss(0, noise * 0.25) for _ in range(max(0, genes - 2))])
        labels.append(label)
        outcome = label if label != "HSC" else rng.choices(["erythroid", "myeloid", "megakaryocyte"], weights=[0.42, 0.38, 0.20])[0]
        outcomes.append(outcome)
    gene_names = [f"synthetic_program_{index:03d}" for index in range(genes)]
    return SyntheticDataset(expression, velocity, embedding, labels, outcomes, [f"synthetic_{i:05d}" for i in range(cells)], {"seed": seed, "potential": "four-well nonlinear benchmark", "temperature": 1.0, "generator": "thermodynamic_waddington.synthetic", "gene_names": gene_names})
