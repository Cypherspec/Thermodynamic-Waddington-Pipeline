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

    # iterative Tarjan so deep graphs do not overflow the recursion limit
    for root in range(n):
        if root in indices:
            continue
        work: list[tuple[int, int]] = [(root, 0)]
        while work:
            node, ptr = work[-1]
            if ptr == 0:
                indices[node] = index
                lowlink[node] = index
                index += 1
                stack.append(node)
                on_stack.add(node)
            neighbors = outgoing[node]
            descended = False
            i = ptr
            while i < len(neighbors):
                target = neighbors[i]
                if target not in indices:
                    work[-1] = (node, i + 1)
                    work.append((target, 0))
                    descended = True
                    break
                if target in on_stack:
                    lowlink[node] = min(lowlink[node], indices[target])
                i += 1
            if descended:
                continue
            if lowlink[node] == indices[node]:
                component = []
                while True:
                    popped = stack.pop()
                    on_stack.discard(popped)
                    component.append(popped)
                    if popped == node:
                        break
                components.append(sorted(component))
            work.pop()
            if work:
                parent = work[-1][0]
                lowlink[parent] = min(lowlink[parent], lowlink[node])
    return sorted(components, key=lambda component: (len(component), component))


def cycle_summary(fit) -> dict[str, object]:
    components = strongly_connected_components(fit.edges, len(fit.energies))
    cyclic = [component for component in components if len(component) > 1]
    return {"components": len(components), "cyclic_components": len(cyclic), "largest_cycle": max((len(c) for c in cyclic), default=0)}
