from __future__ import annotations

from collections import defaultdict
from typing import Sequence

from ..graph import Edge


def strongly_connected_components(edges: Sequence[Edge], n: int) -> list[list[int]]:
    outgoing: dict[int, list[int]] = defaultdict(list)
    for edge in edges:
        if isinstance(edge, dict):
            outgoing[int(edge["source"])].append(int(edge["target"]))
        else:
            outgoing[edge.source].append(edge.target)
    index = 0
    stack: list[int] = []
    on_stack: set[int] = set()
    indices: dict[int, int] = {}
    lowlink: dict[int, int] = {}
    components: list[list[int]] = []

    def visit(node: int) -> None:
        nonlocal index
        indices[node] = index
        lowlink[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)
        for target in outgoing[node]:
            if target not in indices:
                visit(target)
                lowlink[node] = min(lowlink[node], lowlink[target])
            elif target in on_stack:
                lowlink[node] = min(lowlink[node], indices[target])
        if lowlink[node] == indices[node]:
            component = []
            while True:
                target = stack.pop()
                on_stack.remove(target)
                component.append(target)
                if target == node:
                    break
            components.append(sorted(component))

    for node in range(n):
        if node not in indices:
            visit(node)
    return sorted(components, key=lambda component: (len(component), component))


def cycle_summary(fit) -> dict[str, object]:
    components = strongly_connected_components(fit.edges, len(fit.energies))
    cyclic = [component for component in components if len(component) > 1]
    return {"components": len(components), "cyclic_components": len(cyclic), "largest_cycle": max((len(c) for c in cyclic), default=0)}
