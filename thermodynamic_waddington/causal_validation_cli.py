from __future__ import annotations

import argparse
import json
from pathlib import Path

from .causal_validation import StudyConfig, analyze_study, create_preregistration, load_cells


def main() -> int:
    parser = argparse.ArgumentParser(description="Locked donor-paired validation gate for measured intervention studies")
    sub = parser.add_subparsers(dest="command", required=True)
    pre = sub.add_parser("preregister")
    pre.add_argument("--target-fate", required=True)
    pre.add_argument("--control", required=True)
    pre.add_argument("--interventions", nargs="+", required=True)
    pre.add_argument("--out", required=True)
    ana = sub.add_parser("analyze")
    ana.add_argument("--input", required=True)
    ana.add_argument("--target-fate", required=True)
    ana.add_argument("--control", required=True)
    ana.add_argument("--study-id", default="registered-study")
    ana.add_argument("--out", required=True)
    args = parser.parse_args()
    if args.command == "preregister":
        config = StudyConfig(target_fate=args.target_fate, control_intervention=args.control)
        payload = create_preregistration(config, args.interventions)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({"out": args.out, "status": payload["status"], "preregistration_hash": payload["preregistration_hash"]}, indent=2))
        return 0
    config = StudyConfig(target_fate=args.target_fate, control_intervention=args.control)
    report = analyze_study(load_cells(args.input), config, args.study_id)
    report.save(args.out)
    print(json.dumps({"out": args.out, "status": report.status, "claim_boundary": report.claim_boundary, "blockers": len(report.blockers)}, indent=2))
    return 0 if report.status != "invalid_study_schema" else 2


if __name__ == "__main__":
    raise SystemExit(main())
