from __future__ import annotations

import argparse
import json

from .falsification import write_falsification_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Run pre-registered adversarial landscape diagnostics")
    parser.add_argument("--out", default="experiments/falsification_report.json")
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    payload = write_falsification_report(args.out, seed=args.seed)
    print(json.dumps({"out": args.out, "status": "completed", "fingerprint": payload["fingerprint"], "cases": len(payload["cases"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
