from __future__ import annotations

import argparse
import json

from .evidence_graph import write_evidence_graph


def main() -> int:
    parser = argparse.ArgumentParser(description="Build an auditable evidence graph for a Thermodynamic Waddington release")
    parser.add_argument("--root", default=".")
    parser.add_argument("--out", default="experiments/evidence_graph.json")
    args = parser.parse_args()
    payload = write_evidence_graph(args.out, args.root)
    print(json.dumps({"out": args.out, "status": payload["status"], "nodes": payload["counts"]["nodes"], "edges": payload["counts"]["edges"], "blockers": len(payload["blockers"]), "fingerprint": payload["fingerprint"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
