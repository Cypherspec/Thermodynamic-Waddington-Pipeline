from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from .award_capsule import build_capsule, fingerprint


def verify(path: str | Path) -> dict[str, object]:
    target = Path(path)
    payload = json.loads(target.read_text(encoding="utf-8"))
    stored = payload.pop("fingerprint", None)
    computed = fingerprint(payload)
    records = payload.get("artifacts", [])
    missing = [item["path"] for item in records if not item.get("exists")]
    return {"valid_fingerprint": stored == computed, "stored_fingerprint": stored, "computed_fingerprint": computed, "missing_artifacts": missing, "valid": stored == computed and not missing}


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify a Thermodynamic Waddington research capsule")
    parser.add_argument("--verify", required=True)
    args = parser.parse_args()
    report = verify(args.verify)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
