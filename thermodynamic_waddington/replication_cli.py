from __future__ import annotations

import argparse
import json

from .replication import load_studies, write_replication_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate independent-study replication without pooling cells")
    parser.add_argument("input", help="JSON file containing a studies list")
    parser.add_argument("--out", default="experiments/replication_report.json")
    parser.add_argument("--bootstrap", type=int, default=256)
    args = parser.parse_args()
    result = write_replication_report(load_studies(args.input), args.out, args.bootstrap)
    print(json.dumps({"out": args.out, "status": result["status"], "independent_studies": result["independent_studies"], "eligible_studies": result["eligible_studies"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
