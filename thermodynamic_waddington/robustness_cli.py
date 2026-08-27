from __future__ import annotations
import argparse, json
from .robustness import write_robustness_report

def main() -> int:
    parser = argparse.ArgumentParser(description="Run specification-multiverse robustness checks")
    parser.add_argument("--out", default="experiments/robustness_report.json")
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--cells", type=int, default=64)
    parser.add_argument("--genes", type=int, default=8)
    args = parser.parse_args()
    report = write_robustness_report(args.out, seed=args.seed, cells=args.cells, genes=args.genes)
    print(json.dumps({"out": args.out, "status": report["status"], "fingerprint": report["fingerprint"]}, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
