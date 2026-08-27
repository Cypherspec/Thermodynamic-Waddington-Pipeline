from __future__ import annotations

import argparse
import json

from .intervention_trial import TrialConfig, load_observations, write_preregistration, write_trial_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Analyze preregistered measured intervention outcomes")
    sub = parser.add_subparsers(dest="command", required=True)
    pre = sub.add_parser("preregister")
    pre.add_argument("--out", default="experiments/intervention_preregistration.json")
    run = sub.add_parser("analyze")
    run.add_argument("input")
    run.add_argument("--out", default="experiments/intervention_trial_report.json")
    run.add_argument("--preregistration-hash", required=True)
    args = parser.parse_args()
    if args.command == "preregister":
        result = write_preregistration(args.out)
        print(json.dumps({"out": args.out, "preregistration_hash": result["preregistration_hash"]}, indent=2, sort_keys=True))
    else:
        result = write_trial_report(load_observations(args.input), args.out, TrialConfig(), args.preregistration_hash)
        print(json.dumps({"out": args.out, "status": result["status"], "blockers": result["blockers"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
