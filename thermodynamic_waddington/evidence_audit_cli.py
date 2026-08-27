from __future__ import annotations
import argparse, json
from .evidence_audit import write_evidence_audit

def main() -> int:
    parser = argparse.ArgumentParser(description="Audit evidence and claim boundaries across a Thermodynamic Waddington release")
    parser.add_argument("--out", default="experiments/evidence_audit.json")
    args = parser.parse_args()
    payload = write_evidence_audit(args.out)
    print(json.dumps({"out": args.out, "status": payload["status"], "passed": payload["passed"], "total": payload["total"], "fingerprint": payload["fingerprint"]}, indent=2, sort_keys=True))
    return 0 if payload["status"] != "audit_incomplete" else 1

if __name__ == "__main__":
    raise SystemExit(main())
