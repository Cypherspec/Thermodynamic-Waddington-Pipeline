from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, variance
from .config import FitConfig
from .model import LandscapeFit, fit_landscape


@dataclass(frozen=True)
class EnsembleSummary:
    members: int
    mean_energy: tuple[float, ...]
    energy_spread: tuple[float, ...]
    mean_diffusion: tuple[float, ...]
    attractor_consensus: tuple[int, ...]

    def to_dict(self) -> dict[str, object]:
        return {"members": self.members, "mean_energy": list(self.mean_energy), "energy_spread": list(self.energy_spread), "mean_diffusion": list(self.mean_diffusion), "attractor_consensus": list(self.attractor_consensus)}


def fit_ensemble(expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], config: FitConfig, members: int = 8, seed: int = 7) -> tuple[LandscapeFit, ...]:
    if members < 1:
        raise ValueError("members must be positive")
    rng = random.Random(seed)
    fits = []
    for index in range(members):
        indices = [row for row in range(len(expression)) if rng.random() > 0.08]
        if len(indices) < config.neighbors + 2:
            indices = list(range(len(expression)))
        expr = [list(expression[row]) for row in indices]
        vel = [list(velocity[row]) for row in indices]
        member_config = FitConfig(**{**config.__dict__, "seed": seed + index})
        fits.append(fit_landscape(expr, vel, config=member_config, metadata={"ensemble_member": index, "resampled_cells": len(indices)}))
    return tuple(fits)


def summarize_ensemble(fits: Sequence[LandscapeFit], consensus_threshold: float = 0.5) -> EnsembleSummary:
    if not fits:
        return EnsembleSummary(0, (), (), (), ())
    n = len(fits[0].energies)
    energies = [[fit.energies[index] for fit in fits if index < len(fit.energies)] for index in range(n)]
    diffusions = [[fit.diffusions[index] for fit in fits if index < len(fit.diffusions)] for index in range(n)]
    counts: dict[int, int] = {}
    for fit in fits:
        for node in fit.attractors:
            counts[node] = counts.get(node, 0) + 1
    consensus = tuple(sorted(node for node, count in counts.items() if count / len(fits) >= consensus_threshold))
    return EnsembleSummary(len(fits), tuple(mean(values) for values in energies), tuple(variance(values) ** 0.5 for values in energies), tuple(mean(values) for values in diffusions), consensus)
