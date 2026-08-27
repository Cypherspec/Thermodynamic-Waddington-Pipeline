from __future__ import annotations

import argparse
import json
from .reproduce import run_reproduction


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the evidence-bounded reproducibility workflow")
    parser.add_argument("--root", default=".")
    parser.add_argument("--out", default="experiments/reproduction_run.json")
    args = parser.parse_args()
    payload = run_reproduction(args.root, args.out)
    print(json.dumps({"out": args.out, "status": payload["status"], "fingerprint": payload["fingerprint"]}, indent=2, sort_keys=True))
    return 0 if payload["status"] == "reproducible_review_ready" else 1


if __name__ == "__main__":
    raise SystemExit(main())
