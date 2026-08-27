from __future__ import annotations

import argparse
import json
from pathlib import Path

from .validation_suite import run_recovery_suite
from .protocols.nonequilibrium import jarzynski_cumulants
from .benchmarks import benchmark, random_matrix, scaling_summary
from .linalg import pca


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tw-experiment", description="Synthetic recovery and scaling experiments")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--cells", type=int, default=256)
    args = parser.parse_args(argv)
    recovery = run_recovery_suite(args.seed)
    cases = []
    for index, cells in enumerate((max(32, args.cells // 4), max(64, args.cells // 2), args.cells)):
        cases.append(benchmark("pca", lambda matrix: [value for row in pca(matrix, min(5, len(matrix[0])), args.seed)[0] for value in row], cells, 12, 0, args.seed + index))
    payload = {"recovery_suite": recovery.to_dict(), "scaling": scaling_summary(cases), "nonequilibrium_estimator_smoke": jarzynski_cumulants([0.0, 0.2, 0.4, 0.8], 1.0), "scientific_note": "These are algorithmic recovery and scaling checks, not evidence that an inferred landscape is a biological free-energy state function."}
    args.out.write_text(json.dumps(payload, indent=2))
    print(f"wrote experiment report: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
