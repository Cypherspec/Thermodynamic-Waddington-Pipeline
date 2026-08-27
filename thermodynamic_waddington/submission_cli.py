from __future__ import annotations

import argparse
import json
from .review_submission import write_submission_dossier


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate an evidence-bounded review submission dossier")
    parser.add_argument("--out", default="experiments/submission_dossier.json")
    args = parser.parse_args()
    payload = write_submission_dossier(args.out)
    print(json.dumps({"out": args.out, "status": payload["submission_status"], "fingerprint": payload["fingerprint"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
