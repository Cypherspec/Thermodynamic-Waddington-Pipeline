"""Reproduce every benchmark and figure in one command.

Runs each benchmark in sequence. The data-dependent benchmarks are skipped with
a clear message if the pancreas h5ad is not present (see the README Data
section). Synthetic benchmarks always run.

    python benchmarks/run_all.py
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

DATA = Path("data/real/endocrinogenesis_day15.h5ad")

SYNTHETIC = [
    ["benchmarks/runtime_scaling.py"],
    ["benchmarks/plot_scaling.py"],
    ["benchmarks/scaling.py"],
]
REAL_DATA = [
    ["benchmarks/pancreas_entropy_validation.py"],
    ["benchmarks/plot_pancreas_validation.py"],
    ["benchmarks/free_energy_calibration.py"],
    ["benchmarks/irreversibility_detection.py"],
    ["benchmarks/predictive_ordering.py"],
    ["benchmarks/statistical_validation.py"],
]


def run(cmd) -> bool:
    print("::", " ".join(cmd), flush=True)
    start = time.perf_counter()
    result = subprocess.run([sys.executable, *cmd])
    ok = result.returncode == 0
    print(f"   {'ok' if ok else 'FAIL'} in {time.perf_counter() - start:.1f}s", flush=True)
    return ok


def main() -> None:
    ok = True
    for cmd in SYNTHETIC:
        ok = run(cmd) and ok
    if DATA.exists():
        for cmd in REAL_DATA:
            ok = run(cmd) and ok
    else:
        print(f"skip real-data benchmarks: {DATA} not found (see README Data section)")
    print("ALL OK" if ok else "SOME FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
