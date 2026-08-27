from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from .arrays import mean, variance
from .graph import Edge, adjacency


@dataclass(frozen=True)
class ThermodynamicLedger:
    potential_work: float
    housekeeping_heat: float
    excess_heat: float
    entropy_production: float
    current_strength: float
    detailed_balance_violation: float

    def as_dict(self) -> dict:
        return {"potential_work": self.potential_work, "housekeeping_heat": self.housekeeping_heat, "excess_heat": self.excess_heat, "entropy_production": self.entropy_production, "current_strength": self.current_strength, "detailed_balance_violation": self.detailed_balance_violation}


def edge_current(edge: Edge, energies: Sequence[float], temperature: float) -> float:
    forward = math.exp(-max(-50.0, min(50.0, (energies[edge.target] - energies[edge.source]) / max(temperature, 1e-9))))
    reverse = math.exp(max(-50.0, min(50.0, (energies[edge.target] - energies[edge.source]) / max(temperature, 1e-9))))
    return (forward - reverse) * edge.alignment


def ledger(edges: Sequence[Edge], energies: Sequence[float], temperature: float = 1.0) -> ThermodynamicLedger:
    if not edges:
        return ThermodynamicLedger(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    currents = [edge_current(edge, energies, temperature) for edge in edges]
    potential = sum(energies[edge.target] - energies[edge.source] for edge in edges) / len(edges)
    housekeeping = sum(abs(current) * max(0.0, energies[edge.target] - energies[edge.source]) for current, edge in zip(currents, edges)) / len(edges)
    excess = sum(abs(current) for current in currents) / len(currents)
    entropy = sum(current * current for current in currents) / len(currents)
    violation = variance(currents) ** 0.5
    return ThermodynamicLedger(potential, housekeeping, excess, entropy, mean([abs(value) for value in currents]), violation)


def local_dissipation(edges: Sequence[Edge], energies: Sequence[float], n: int, temperature: float = 1.0) -> list[float]:
    result = [0.0 for _ in range(n)]
    counts = [0 for _ in range(n)]
    for edge in edges:
        value = edge_current(edge, energies, temperature) ** 2
        result[edge.source] += value
        counts[edge.source] += 1
    return [value / max(1, count) for value, count in zip(result, counts)]


def effective_thermodynamic_summary(energies: Sequence[float], diffusions: Sequence[float], edges: Sequence[Edge], temperature: float = 1.0) -> dict[str, object]:
    record = ledger(edges, energies, temperature)
    local = local_dissipation(edges, energies, len(energies), temperature)
    return {"ledger": record.as_dict(), "local_dissipation_mean": mean(local) if local else 0.0, "local_dissipation_max": max(local) if local else 0.0, "diffusion_mean": mean(diffusions) if diffusions else 0.0, "interpretation": "effective stochastic diagnostic; not a calibrated molecular thermodynamic measurement"}
