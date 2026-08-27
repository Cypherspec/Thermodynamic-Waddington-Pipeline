from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class PersistencePair:
    minimum: int
    saddle: int
    birth: float
    death: float
    persistence: float

    def to_dict(self) -> dict[str, float | int]:
        return {"minimum": self.minimum, "saddle": self.saddle, "birth": self.birth, "death": self.death, "persistence": self.persistence}


def lower_star_pairs(energies: Sequence[float], edges: Sequence[tuple[int, int]]) -> list[PersistencePair]:
    order = sorted(range(len(energies)), key=lambda index: energies[index])
    active: set[int] = set()
    components: dict[int, int] = {}
    pairs: list[PersistencePair] = []

    def root(component: int) -> int:
        while components.get(component, component) != component:
            component = components[component]
        return component

    for node in order:
        active.add(node)
        components[node] = node
        adjacent = [other for left, right in edges for other in ((right,) if left == node else (left,) if right == node else ()) if other in active]
        for other in adjacent:
            left = root(node)
            right = root(other)
            if left == right:
                continue
            if energies[left] <= energies[right]:
                components[right] = left
                pairs.append(PersistencePair(right, node, energies[right], energies[node], energies[node] - energies[right]))
            else:
                components[left] = right
                pairs.append(PersistencePair(left, node, energies[left], energies[node], energies[node] - energies[left]))
    return sorted(pairs, key=lambda pair: pair.persistence, reverse=True)
