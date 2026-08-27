from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class Basin:
    representative: int
    members: tuple[int, ...]
    minimum_energy: float
    mass: float

    def to_dict(self) -> dict[str, object]:
        return {
            "representative": self.representative,
            "members": list(self.members),
            "minimum_energy": self.minimum_energy,
            "mass": self.mass,
        }


def watershed_basins(energies: Sequence[float], edges: Sequence[tuple[int, int]], tolerance: float = 0.0) -> list[Basin]:
    outgoing: list[list[int]] = [[] for _ in energies]
    for source, target in edges:
        outgoing[source].append(target)
    roots = []
    for node, energy in enumerate(energies):
        if all(energies[target] >= energy - tolerance for target in outgoing[node]):
            roots.append(node)
    if not roots:
        return []
    members = {root: [] for root in roots}
    for node in range(len(energies)):
        nearest = min(roots, key=lambda root: abs(energies[root] - energies[node]))
        members[nearest].append(node)
    return [Basin(root, tuple(members[root]), float(energies[root]), len(members[root]) / max(1, len(energies))) for root in roots]

def barrier_heights(energies: Sequence[float], edges: Sequence[tuple[int, int]]) -> list[float]:
    """Return directed uphill barriers for a compact topology representation."""
    return [max(0.0, float(energies[target]) - float(energies[source])) for source, target in edges]

