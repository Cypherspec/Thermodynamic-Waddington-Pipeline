from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

from .arrays import mean, variance
from .config import FitConfig
from .model import LandscapeFit, fit_landscape


@dataclass
class OnlineLandscape:
    """Bounded-memory mini-batch accumulator for time-resolved experiments."""

    config: FitConfig = field(default_factory=FitConfig)
    expression: list[list[float]] = field(default_factory=list)
    velocity: list[list[float]] = field(default_factory=list)
    embedding: list[list[float]] | None = None
    labels: list[str] | None = None
    cell_ids: list[str] | None = None
    batches_seen: int = 0
    observations_seen: int = 0

    def update(self, expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], embedding: Sequence[Sequence[float]] | None = None, labels: Sequence[str] | None = None, cell_ids: Sequence[str] | None = None) -> None:
        if len(expression) != len(velocity):
            raise ValueError("expression and velocity batch sizes must match")
        if not expression:
            return
        self.expression.extend([list(row) for row in expression])
        self.velocity.extend([list(row) for row in velocity])
        if embedding is not None:
            if self.embedding is None:
                self.embedding = []
            self.embedding.extend([list(row) for row in embedding])
        if labels is not None:
            if self.labels is None:
                self.labels = []
            self.labels.extend(labels)
        if cell_ids is not None:
            if self.cell_ids is None:
                self.cell_ids = []
            self.cell_ids.extend(cell_ids)
        self.batches_seen += 1
        self.observations_seen += len(expression)

    def snapshot(self) -> LandscapeFit:
        return fit_landscape(self.expression, self.velocity, self.embedding, self.labels, self.cell_ids, self.config, {"mode": "online", "batches_seen": self.batches_seen, "observations_seen": self.observations_seen})

    def diagnostics(self) -> dict[str, float | int]:
        if not self.velocity:
            return {"batches_seen": self.batches_seen, "observations_seen": self.observations_seen, "mean_velocity": 0.0, "velocity_variance": 0.0}
        speeds = [sum(value * value for value in row) ** 0.5 for row in self.velocity]
        return {"batches_seen": self.batches_seen, "observations_seen": self.observations_seen, "mean_velocity": mean(speeds), "velocity_variance": variance(speeds)}

    def clear(self) -> None:
        self.expression.clear()
        self.velocity.clear()
        if self.embedding is not None:
            self.embedding.clear()
        if self.labels is not None:
            self.labels.clear()
        if self.cell_ids is not None:
            self.cell_ids.clear()
        self.batches_seen = 0
        self.observations_seen = 0
