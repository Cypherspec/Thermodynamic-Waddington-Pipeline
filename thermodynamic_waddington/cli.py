from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .config import FitConfig
from .io import read_bundle
from .model import LandscapeFit, fit_landscape
from .render import render_report
from .synthetic import SyntheticDataset, make_synthetic_dataset


def _fit_from_bundle(bundle: dict[str, object], args: argparse.Namespace) -> LandscapeFit:
    config = FitConfig(neighbors=args.neighbors, dimensions=args.dimensions, bootstrap_replicates=args.bootstrap, seed=args.seed, reference_index=args.reference, temperature=args.temperature)
    return fit_landscape(bundle["expression"], bundle["velocity"], config=config, embedding=bundle.get("embedding"), labels=bundle.get("labels"), cell_ids=bundle.get("cell_ids"), lineage_outcomes=bundle.get("lineage_outcomes"), metadata={**(bundle.get("metadata") or {}), "gene_names": bundle.get("gene_names", (bundle.get("metadata") or {}).get("gene_names", []))})


def command_simulate(args: argparse.Namespace) -> int:
    dataset = make_synthetic_dataset(args.cells, args.genes, args.seed)
    if args.gene_names:
        names = [line.strip() for line in args.gene_names.read_text().splitlines() if line.strip()]
        if len(names) != args.genes:
            raise ValueError("--gene-names must contain exactly --genes names")
        payload = dataset.to_dict()
        payload["gene_names"] = names
        payload["metadata"]["gene_names"] = names
        Path(args.out).write_text(json.dumps(payload, separators=(",", ":")))
    else:
        dataset.save(args.out)
    print(f"wrote synthetic dataset: {args.out}")
    return 0


def command_fit(args: argparse.Namespace) -> int:
    fit = _fit_from_bundle(read_bundle(args.input), args)
    fit.save(args.out)
    print(json.dumps({"out": str(args.out), "cells": len(fit.energies), "edges": len(fit.edges), "attractors": len(fit.attractors), "coverage": fit.diagnostics["coverage"]}, indent=2))
    return 0


def command_report(args: argparse.Namespace) -> int:
    fit = LandscapeFit.load(args.input)
    render_report(fit, args.out)
    print(f"wrote report: {args.out}")
    return 0


def command_audit(args: argparse.Namespace) -> int:
    fit = LandscapeFit.load(args.input)
    audit = fit.diagnostics.get("scientific_audit", {})
    print(json.dumps(audit, indent=2))
    return 0

def command_validate(args: argparse.Namespace) -> int:
    fit = LandscapeFit.load(args.input)
    checks = {"coverage": fit.diagnostics.get("coverage"), "mean_velocity_alignment": fit.diagnostics.get("mean_velocity_alignment"), "reversibility_gap": fit.diagnostics.get("reversibility_gap"), "attractors": fit.attractors, "status": "diagnostic-only; inspect assumptions before biological interpretation"}
    print(json.dumps(checks, indent=2))
    return 0


def parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="twaddington", description="Estimate and audit an effective cell-fate free-energy landscape from velocity fields")
    sub = parser.add_subparsers(dest="command", required=True)
    simulate = sub.add_parser("simulate", help="generate a nonlinear synthetic benchmark")
    simulate.add_argument("--cells", type=int, default=800)
    simulate.add_argument("--genes", type=int, default=16)
    simulate.add_argument("--seed", type=int, default=17)
    simulate.add_argument("--out", type=Path, required=True)
    simulate.add_argument("--gene-names", type=Path, help="optional newline-delimited real gene names")
    simulate.set_defaults(func=command_simulate)
    fit = sub.add_parser("fit", help="fit a landscape from a JSON bundle")
    fit.add_argument("--input", type=Path, required=True)
    fit.add_argument("--out", type=Path, required=True)
    fit.add_argument("--neighbors", type=int, default=24)
    fit.add_argument("--dimensions", type=int, default=6)
    fit.add_argument("--bootstrap", type=int, default=24)
    fit.add_argument("--seed", type=int, default=17)
    fit.add_argument("--reference", type=int, default=0)
    fit.add_argument("--temperature", type=float, default=1.0)
    fit.set_defaults(func=command_fit)
    report = sub.add_parser("report", help="render an interactive-free HTML/SVG report")
    report.add_argument("--input", type=Path, required=True)
    report.add_argument("--out", type=Path, required=True)
    report.set_defaults(func=command_report)
    audit = sub.add_parser("audit", help="print scientific identifiability and claim-boundary audit")
    audit.add_argument("--input", type=Path, required=True)
    audit.set_defaults(func=command_audit)
    validate = sub.add_parser("validate", help="print diagnostic checks")
    validate.add_argument("--input", type=Path, required=True)
    validate.set_defaults(func=command_validate)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, OSError, KeyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
