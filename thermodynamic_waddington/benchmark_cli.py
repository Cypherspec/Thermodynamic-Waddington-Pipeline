from __future__ import annotations

import argparse
import json
from pathlib import Path

from .benchmarking import compare_energy_methods
from .config import FitConfig
from .io import read_bundle
from .model import fit_landscape
from .synthetic import make_synthetic_dataset


def _fit(dataset, args):
    config = FitConfig(neighbors=args.neighbors, dimensions=args.dimensions, bootstrap_replicates=args.bootstrap, seed=args.seed)
    return fit_landscape(dataset.expression, dataset.velocity, config=config, embedding=dataset.embedding, labels=dataset.labels, cell_ids=dataset.cell_ids, lineage_outcomes=dataset.lineage_outcomes)


def _density_energy(expression: list[list[float]]) -> list[float]:
    row_sums = [sum(max(0.0, value) for value in row) for row in expression]
    total = sum(row_sums) or 1.0
    return [-value / total for value in row_sums]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tw-benchmark", description="Compare landscape estimators on synthetic or bundle data")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cells", type=int, default=128)
    parser.add_argument("--genes", type=int, default=12)
    parser.add_argument("--neighbors", type=int, default=12)
    parser.add_argument("--dimensions", type=int, default=5)
    parser.add_argument("--bootstrap", type=int, default=4)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args(argv)
    if args.input:
        payload = read_bundle(args.input)
        class Dataset:
            expression = payload["expression"]
            velocity = payload["velocity"]
            embedding = payload.get("embedding")
            labels = payload.get("labels")
            cell_ids = payload.get("cell_ids")
            lineage_outcomes = payload.get("lineage_outcomes")
        dataset = Dataset()
    else:
        dataset = make_synthetic_dataset(args.cells, args.genes, args.seed)
    fit = _fit(dataset, args)
    reference = fit.energies
    candidates = {"density_only": _density_energy(dataset.expression), "velocity_landscape": fit.energies}
    results = compare_energy_methods(reference, candidates, set(fit.attractors))
    payload = {"results": [item.to_dict() for item in results], "reference": "velocity_landscape", "scientific_note": "Benchmark rankings are dataset- and preprocessing-dependent; they are not evidence of biological truth."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2))
    print(f"wrote benchmark report: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
