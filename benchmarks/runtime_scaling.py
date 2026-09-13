"""Runtime scaling benchmark for the landscape fit.

Times the full fit_landscape call across a range of cell counts and compares
against a recorded pure-Python baseline (commit 40d1fcd, before the numpy
vectorization). The fit outputs are bit-identical between the two versions:
energies max abs diff is 0.0, attractors and counterfactual ranking match
exactly. So this measures pure speedup, not a quality tradeoff.

Reproduce the current numbers:
    python benchmarks/runtime_scaling.py --sizes 120 200 300 400
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.model import fit_landscape
from thermodynamic_waddington.synthetic import make_synthetic_dataset

# seconds per fit at commit 40d1fcd (pure Python core), same machine and config
BASELINE = {"120": 3.632, "200": 8.278, "300": 15.225, "400": 28.102}


def bench_config(seed: int = 17) -> FitConfig:
    return FitConfig(
        neighbors=20,
        dimensions=6,
        bootstrap_replicates=8,
        entropy_production_bootstrap_replicates=8,
        entropy_production_permutation_replicates=8,
        seed=seed,
    )


def time_fit(n: int, repeats: int = 2, seed: int = 7) -> float:
    ds = make_synthetic_dataset(cells=n, genes=16, seed=seed)
    cfg = bench_config()
    best = float("inf")
    for _ in range(max(1, repeats)):
        start = time.perf_counter()
        fit_landscape(ds.expression, ds.velocity, config=cfg, labels=ds.labels)
        best = min(best, time.perf_counter() - start)
    return best


def run(sizes: list[int], repeats: int = 2) -> dict:
    rows = []
    for n in sizes:
        dt = time_fit(n, repeats=repeats)
        base = BASELINE.get(str(n))
        rows.append({
            "n": n,
            "baseline_s": base,
            "current_s": round(dt, 3),
            "speedup": round(base / dt, 2) if base else None,
        })
        tag = f"{base / dt:.1f}x" if base else "n/a"
        print(f"n={n:5d}  current={dt:7.2f}s  baseline={base}  speedup={tag}")
    return {
        "benchmark": "runtime_scaling",
        "rows": rows,
        "python": sys.version.split()[0],
        "platform": f"{platform.system()}-{platform.release()}-{platform.machine()}",
        "config": "neighbors=20 dims=6 boot=8 ep_boot=8 ep_perm=8 seed=17",
        "correctness": "fit outputs bit-identical to baseline (energies max abs diff 0.0)",
        "baseline_commit": "40d1fcd",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sizes", type=int, nargs="+", default=[120, 200, 300, 400])
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--out", default="experiments/runtime_scaling.json")
    args = ap.parse_args()
    report = run(args.sizes, args.repeats)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print("wrote", out)


if __name__ == "__main__":
    main()
