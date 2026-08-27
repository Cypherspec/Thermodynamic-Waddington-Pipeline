from __future__ import annotations

import argparse
import json

from .decision_engine import DecisionConfig, decide_from_files


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the evidence-gated Thermodynamic Waddington decision engine")
    parser.add_argument("--trial")
    parser.add_argument("--transportability")
    parser.add_argument("--robustness")
    parser.add_argument("--path-ensemble")
    parser.add_argument("--out", default="experiments/decision_report.json")
    parser.add_argument("--minimum-effect", type=float, default=0.10)
    args = parser.parse_args()
    payload = decide_from_files(
        trial_path=args.trial,
        transportability_path=args.transportability,
        robustness_path=args.robustness,
        path_ensemble_path=args.path_ensemble,
        config=DecisionConfig(minimum_effect=args.minimum_effect),
    )
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"out": args.out, "decision": payload["decision"], "attained_level": payload["attained_level"], "blockers": len(payload["blockers"]), "fingerprint": payload["fingerprint"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
