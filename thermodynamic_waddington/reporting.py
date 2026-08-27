from __future__ import annotations

import json
from dataclasses import asdict
from html import escape
from pathlib import Path

from .model import LandscapeFit


def fit_manifest(fit: LandscapeFit) -> dict:
    return {
        "schema": "twaddington.manifest.v1",
        "metadata": fit.metadata,
        "config": fit.config,
        "n_cells": len(fit.energies),
        "n_edges": len(fit.edges),
        "attractor_count": len(fit.attractors),
        "diagnostics": fit.diagnostics,
    }


def write_manifest(fit: LandscapeFit, path: str | Path) -> None:
    Path(path).write_text(json.dumps(fit_manifest(fit), indent=2, sort_keys=True))


def markdown_summary(fit: LandscapeFit) -> str:
    lines = ["# Landscape fit", "", f"- Cells: {len(fit.energies)}", f"- Directed edges: {len(fit.edges)}", f"- Candidate attractors: {len(fit.attractors)}", "", "## Diagnostics", ""]
    lines.extend(f"- **{key}**: {value}" for key, value in fit.diagnostics.items())
    lines.extend(["", "## Scientific status", "", "This is an effective, model-dependent landscape estimate. It is not an equilibrium free energy and should not be interpreted as a direct thermodynamic state function."])
    return "\n".join(lines) + "\n"


def write_markdown_summary(fit: LandscapeFit, path: str | Path) -> None:
    Path(path).write_text(markdown_summary(fit))


def graphviz_dot(fit: LandscapeFit, max_edges: int = 600) -> str:
    lines = ["digraph landscape {", "  graph [rankdir=LR, bgcolor=transparent];", "  node [shape=circle, style=filled, fontname=Helvetica];"]
    attractors = set(fit.attractors)
    for index, energy in enumerate(fit.energies):
        color = "#ffb000" if index in attractors else "#4f87ff"
        lines.append(f'  n{index} [label="{index}", fillcolor="{color}", tooltip="energy={energy:.4f}"];')
    for edge in fit.edges[:max_edges]:
        lines.append(f'  n{edge["source"]} -> n{edge["target"]} [label="{edge["work"]:.3f}", color="#7d8ca3"];')
    lines.append("}")
    return "\n".join(lines) + "\n"


def write_graphviz(fit: LandscapeFit, path: str | Path, max_edges: int = 600) -> None:
    Path(path).write_text(graphviz_dot(fit, max_edges))
