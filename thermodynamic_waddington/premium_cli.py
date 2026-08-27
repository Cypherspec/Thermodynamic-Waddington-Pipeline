from __future__ import annotations

import argparse
import json
from pathlib import Path

from .premium_platform import run_premium_platform, write_premium_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="tw-premium", description="Run the evidence-gated premium research platform")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--out", type=Path, default=Path("experiments/premium_platform_run.json"))
    parser.add_argument("--no-refresh", action="store_true", help="reuse existing computational artifacts")
    args = parser.parse_args(argv)
    report = run_premium_platform(args.root, refresh=not args.no_refresh)
    if args.out != args.root / "experiments/premium_platform_run.json":
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(args.out), "status": report["status"], "claim_level": report["claim_level"], "blockers": len(report["blockers"]), "fingerprint": report["fingerprint"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
