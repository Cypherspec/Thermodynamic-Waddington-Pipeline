from __future__ import annotations

import argparse
import json
from pathlib import Path

from .discovery_engine import DiscoveryConfig, build_discovery_ledger, validate_public_references


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tw-discover", description="Create conservative, preregisterable fate hypotheses")
    parser.add_argument("--out", type=Path, default=Path("experiments/discovery_ledger.json"))
    parser.add_argument("--larry-root", type=Path, default=Path("data/real/larry"))
    parser.add_argument("--permutations", type=int, default=499)
    args = parser.parse_args(argv)
    ledger = build_discovery_ledger(config=DiscoveryConfig(permutation_count=args.permutations))
    payload = ledger.to_dict()
    payload["public_reference_validation"] = validate_public_references(args.larry_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"out": str(args.out), "status": payload["claim_status"], "hypotheses": len(payload["hypotheses"]), "blockers": len(payload["blockers"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
