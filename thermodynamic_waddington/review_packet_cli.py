from __future__ import annotations
import argparse, json
from .review_packet import write_review_packet

def main() -> int:
    p = argparse.ArgumentParser(description="Generate an evidence-bounded reviewer packet")
    p.add_argument("--out", default="experiments/review_packet.json")
    a = p.parse_args()
    report = write_review_packet(a.out)
    print(json.dumps({"out": a.out, "status": report["status"], "fingerprint": report["fingerprint"]}, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
