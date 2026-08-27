from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import read_bundle, write_bundle
from .model import fit_landscape
from .real_data import anndata_to_bundle, catalog_payload, download_dataset, inspect_anndata, paul15_to_bundle
from .benchmark_protocol import run_external_validation_report
from .causal_validation import StudyConfig, analyze_study, create_preregistration, load_cells
from .wetlab_protocol import build_hematopoietic_fate_protocol, write_protocol


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="tw-real", description="Opt-in real-data acquisition and AnnData conversion")
    sub = parser.add_subparsers(dest="command", required=True)
    catalog = sub.add_parser("catalog", help="print the curated dataset catalog")
    catalog.add_argument("--out", type=Path)
    catalog.set_defaults(handler=cmd_catalog)
    download = sub.add_parser("download", help="download one explicitly selected public dataset")
    download.add_argument("name", choices=["paul15", "pancreas", "larry_in_vitro"])
    download.add_argument("--out", type=Path, required=True)
    download.set_defaults(handler=cmd_download)
    paul = sub.add_parser("paul15-bundle", help="convert the downloaded Paul15 HDF5 expression matrix into a project bundle")
    paul.add_argument("path", type=Path)
    paul.add_argument("--out", type=Path, required=True)
    paul.add_argument("--max-cells", type=int)
    paul.add_argument("--max-genes", type=int)
    paul.set_defaults(handler=cmd_paul15_bundle)
    inspect = sub.add_parser("inspect", help="inspect an AnnData/H5AD file without loading its matrix")
    inspect.add_argument("path", type=Path)
    inspect.set_defaults(handler=cmd_inspect)
    convert = sub.add_parser("convert", help="convert AnnData into the project JSON bundle")
    convert.add_argument("path", type=Path)
    convert.add_argument("--out", type=Path, required=True)
    convert.add_argument("--expression-key", default="X")
    convert.add_argument("--velocity-key", default="velocity")
    convert.add_argument("--label-key")
    convert.add_argument("--max-cells", type=int)
    convert.add_argument("--max-genes", type=int)
    convert.set_defaults(handler=cmd_convert)
    fit = sub.add_parser("fit", help="convert AnnData and fit a landscape in one explicit command")
    fit.add_argument("path", type=Path)
    fit.add_argument("--out", type=Path, required=True)
    fit.add_argument("--expression-key", default="X")
    fit.add_argument("--velocity-key", default="velocity")
    fit.add_argument("--label-key")
    fit.add_argument("--max-cells", type=int)
    fit.add_argument("--max-genes", type=int)
    fit.add_argument("--neighbors", type=int, default=24)
    fit.add_argument("--dimensions", type=int, default=6)
    fit.add_argument("--bootstrap", type=int, default=24)
    fit.add_argument("--seed", type=int, default=17)
    fit.add_argument("--reference", type=int, default=0)
    fit.add_argument("--temperature", type=float, default=1.0)
    fit.set_defaults(handler=cmd_fit)
    validate = sub.add_parser("validate-external", help="run the reproducible public pancreas and LARRY validation protocol")
    validate.add_argument("--out", type=Path, default=Path("experiments/external_validation.json"))
    validate.add_argument("--pancreas", type=Path, default=Path("data/real/endocrinogenesis_day15.h5ad"))
    validate.add_argument("--larry-root", type=Path, default=Path("data/real/larry"))
    validate.set_defaults(handler=cmd_validate_external)
    prereg = sub.add_parser("preregister", help="write a locked wet-lab protocol and causal-analysis preregistration")
    prereg.add_argument("--target-fate", required=True)
    prereg.add_argument("--control", default="vehicle")
    prereg.add_argument("--intervention", action="append", required=True)
    prereg.add_argument("--protocol-id", default="TW-HEM-FATE-001")
    prereg.add_argument("--out", type=Path, required=True)
    prereg.set_defaults(handler=cmd_preregister)
    analyze = sub.add_parser("analyze-study", help="analyze measured donor-paired perturbation outcomes without upgrading unsupported claims")
    analyze.add_argument("path", type=Path)
    analyze.add_argument("--target-fate", required=True)
    analyze.add_argument("--control", default="vehicle")
    analyze.add_argument("--out", type=Path, required=True)
    analyze.add_argument("--minimum-donors", type=int, default=3)
    analyze.add_argument("--minimum-cells-per-group", type=int, default=10)
    analyze.add_argument("--minimum-effect", type=float, default=0.10)
    analyze.add_argument("--bootstrap", type=int, default=2000)
    analyze.add_argument("--seed", type=int, default=17)
    analyze.set_defaults(handler=cmd_analyze_study)
    return parser


def cmd_catalog(args: argparse.Namespace) -> int:
    payload = catalog_payload()
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(json.dumps(payload, indent=2))
    return 0


def cmd_download(args: argparse.Namespace) -> int:
    print(json.dumps(download_dataset(args.name, args.out).to_dict(), indent=2))
    return 0


def cmd_paul15_bundle(args: argparse.Namespace) -> int:
    payload = paul15_to_bundle(args.path, args.max_cells, args.max_genes)
    write_bundle(args.out, payload)
    print(json.dumps({"out": str(args.out), "cells": len(payload["expression"]), "genes": len(payload["gene_names"]), "velocity_status": payload["metadata"]["velocity_status"]}, indent=2))
    return 0


def cmd_inspect(args: argparse.Namespace) -> int:
    print(json.dumps(inspect_anndata(args.path), indent=2))
    return 0


def cmd_convert(args: argparse.Namespace) -> int:
    payload = anndata_to_bundle(args.path, args.expression_key, args.velocity_key, args.label_key, args.max_cells, args.max_genes)
    write_bundle(args.out, payload)
    print(f"wrote bundle: {args.out}")
    return 0


def cmd_fit(args: argparse.Namespace) -> int:
    payload = anndata_to_bundle(args.path, args.expression_key, args.velocity_key, args.label_key, args.max_cells, args.max_genes)
    from .config import FitConfig
    fit = fit_landscape(payload["expression"], payload["velocity"], FitConfig(neighbors=args.neighbors, dimensions=args.dimensions, bootstrap_replicates=args.bootstrap, seed=args.seed, reference_index=args.reference, temperature=args.temperature), labels=payload.get("labels"), cell_ids=payload.get("cell_ids"), metadata={**(payload.get("metadata") or {}), "gene_names": payload.get("gene_names", (payload.get("metadata") or {}).get("gene_names", []))})
    fit.save(args.out)
    print(json.dumps({"out": str(args.out), "cells": len(fit.energies), "edges": len(fit.edges), "coverage": fit.diagnostics.get("coverage")}, indent=2))
    return 0


def cmd_validate_external(args: argparse.Namespace) -> int:
    report = run_external_validation_report(args.out, args.pancreas, args.larry_root)
    print(json.dumps({"out": str(args.out), "status": report["status"], "benchmarks": [item["name"] for item in report["results"]]}, indent=2))
    return 0

def cmd_preregister(args: argparse.Namespace) -> int:
    config = StudyConfig(args.target_fate, args.control, bootstrap_rounds=max(100, args.bootstrap if hasattr(args, "bootstrap") else 2000))
    protocol = build_hematopoietic_fate_protocol(args.target_fate, args.intervention)
    payload = {"protocol": protocol.to_dict(), "analysis": create_preregistration(config, args.intervention, args.protocol_id)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(json.dumps({"out": str(args.out), "protocol_id": args.protocol_id, "status": payload["analysis"]["status"]}, indent=2))
    return 0


def cmd_analyze_study(args: argparse.Namespace) -> int:
    cells = load_cells(args.path)
    report = analyze_study(cells, StudyConfig(args.target_fate, args.control, minimum_donors=args.minimum_donors, minimum_cells_per_group=args.minimum_cells_per_group, minimum_effect=args.minimum_effect, bootstrap_rounds=args.bootstrap, seed=args.seed), args.path.stem)
    report.save(args.out)
    print(json.dumps({"out": str(args.out), "status": report.status, "interventions": len(report.interventions), "blockers": len(report.blockers)}, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.handler(args)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"error: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
