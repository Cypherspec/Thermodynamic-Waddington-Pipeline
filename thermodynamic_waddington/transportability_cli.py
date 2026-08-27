from __future__ import annotations

import argparse
import json

from .transportability import load_transportability_manifest, write_transportability_report


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate fixed-predictor transportability across independent studies")
    parser.add_argument("manifest")
    parser.add_argument("--out", default="experiments/transportability.json")
    args = parser.parse_args()
    studies, config = load_transportability_manifest(args.manifest)
    payload = write_transportability_report(studies, args.out, config)
    print(json.dumps({"out": args.out, "status": payload["status"], "fingerprint": payload["fingerprint"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
