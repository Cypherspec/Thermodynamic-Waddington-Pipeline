from __future__ import annotations

import argparse
import json

from .translation_gate import write_translation_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate whether evidence supports a biological claim")
    parser.add_argument("--root", default=".")
    parser.add_argument("--out", default="experiments/translation_gate.json")
    args = parser.parse_args()
    payload = write_translation_report(args.out, args.root)
    decision = payload["decision"]
    print(json.dumps({"out": args.out, "status": decision["status"], "attained_level": decision["attained_level"], "target_level": decision["target_level"], "blockers": len(decision["blockers"]), "fingerprint": payload["fingerprint"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
