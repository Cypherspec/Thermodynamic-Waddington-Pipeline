from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, variance
from .graph import Edge, adjacency


@dataclass(frozen=True)
class ControlPolicy:
    actions: list[int]
    expected_energy_drop: float
    robustness: float
    intervention_cost: float

    def to_dict(self) -> dict[str, object]:
        return {"actions": self.actions, "expected_energy_drop": self.expected_energy_drop, "robustness": self.robustness, "intervention_cost": self.intervention_cost}


def optimize_policy(edges: Sequence[Edge], energies: Sequence[float], target: int, horizon: int = 6, budget: float = 3.0) -> ControlPolicy:
    outgoing = adjacency(edges, len(energies))
    current = target
    actions: list[int] = []
    drop = 0.0
    cost = 0.0
    for _ in range(horizon):
        candidates = [edge for edge in outgoing[current] if energies[edge.target] < energies[current]]
        if not candidates:
            break
        edge = min(candidates, key=lambda item: energies[item.target])
        local_cost = max(edge.distance, 1e-9)
        if cost + local_cost > budget:
            break
        actions.append(edge.target)
        drop += energies[current] - energies[edge.target]
        cost += local_cost
        current = edge.target
    robustness = 1.0 / (1.0 + variance([energies[node] for node in actions]) if actions else 1.0)
    return ControlPolicy(actions, drop, robustness, cost)


def policy_rollout(policy: ControlPolicy, energies: Sequence[float], perturbation: float = 0.05) -> dict[str, float]:
    if not policy.actions:
        return {"nominal_drop": 0.0, "perturbed_drop": 0.0, "retention": 1.0}
    nominal = sum(max(0.0, energies[policy.actions[i - 1]] - energies[node]) for i, node in enumerate(policy.actions) if i)
    perturbed = max(0.0, nominal - perturbation * math.sqrt(len(policy.actions)))
    return {"nominal_drop": nominal, "perturbed_drop": perturbed, "retention": perturbed / max(nominal, 1e-9)}
