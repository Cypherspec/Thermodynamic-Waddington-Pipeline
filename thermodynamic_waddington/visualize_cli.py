from __future__ import annotations

import argparse
import json
from pathlib import Path

from .model import LandscapeFit
from .reporting import write_graphviz, write_manifest, write_markdown_summary
from .visual_diagnostics import diagnostics_payload
from .visual_queue import render_queue_html, queue_from_fit_stages, queue_terminal_summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tw-visualize", description="Generate visual, graph, manifest, diagnostic, and queue artifacts")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--title", default="Thermodynamic Waddington / Visual Queue")
    parser.add_argument("--max-edges", type=int, default=600)
    args = parser.parse_args(argv)
    fit = LandscapeFit.load(args.input)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    visual_diagnostics = diagnostics_payload(fit)
    (args.out_dir / "visual_diagnostics.json").write_text(json.dumps(visual_diagnostics, indent=2, sort_keys=True))
    stages = [
        ("load", "Loaded serialized LandscapeFit", {"cells": len(fit.energies), "edges": len(fit.edges)}),
        ("diagnostics", "Collected model diagnostics", {"keys": sorted(fit.diagnostics), "overall": visual_diagnostics["overall"]}),
        ("geometry", "Prepared embedding and energy coordinates", {"dimensions": len(fit.embedding[0]) if fit.embedding else 0}),
        ("topology", "Prepared attractor and basin overlays", {"attractors": len(fit.attractors)}),
        ("export", "Wrote visual artifacts", {"output_dir": str(args.out_dir)}),
    ]
    queue = queue_from_fit_stages(stages)
    queue.save(args.out_dir / "visual_queue.json")
    render_queue_html(queue, args.out_dir / "visual_queue.html", title=args.title)
    write_graphviz(fit, args.out_dir / "landscape.dot", max_edges=args.max_edges)
    write_manifest(fit, args.out_dir / "fit_manifest.json")
    write_markdown_summary(fit, args.out_dir / "fit_summary.md")
    print(queue_terminal_summary(queue))
    print(f"wrote visual artifacts: {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
