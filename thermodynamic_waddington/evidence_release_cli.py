from __future__ import annotations

import argparse
import json

from .evidence_release import write_release_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Build an integrity-checked, evidence-bounded release manifest")
    parser.add_argument("--root", default=".")
    parser.add_argument("--out", default="experiments/evidence_release.json")
    args = parser.parse_args()
    report = write_release_report(args.out, args.root)
    print(json.dumps({"out": args.out, "status": report["status"], "passed_checks": report["passed_checks"], "total_checks": report["total_checks"], "fingerprint": report["fingerprint"]}, indent=2, sort_keys=True))
    return 0 if report["status"] == "release_ready_evidence_bounded" else 1


if __name__ == "__main__":
    raise SystemExit(main())
