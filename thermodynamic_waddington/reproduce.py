from __future__ import annotations

"""One-command, evidence-bounded reproducibility runner."""

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .award_capsule import fingerprint, write_capsule
from .review_packet import write_review_packet
from .falsification import write_falsification_report
from .evidence_audit import write_evidence_audit
from .review_submission import write_submission_dossier


def _run(command: list[str], root: Path) -> dict[str, Any]:
    completed = subprocess.run(command, cwd=root, text=True, capture_output=True)
    return {
        "command": " ".join(command),
        "returncode": completed.returncode,
        "passed": completed.returncode == 0,
        "stdout_tail": completed.stdout[-2000:],
        "stderr_tail": completed.stderr[-2000:],
    }


def run_reproduction(root: str | Path = ".", output: str | Path = "experiments/reproduction_run.json") -> dict[str, Any]:
    root = Path(root).resolve()
    steps = [
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-q"],
        [sys.executable, "-m", "compileall", "-q", "thermodynamic_waddington"],
    ]
    results = [_run(command, root) for command in steps]
    capsule = write_capsule(root / "experiments/research_capsule.json", root)
    packet = write_review_packet(root / "experiments/review_packet.json", root)
    falsification = write_falsification_report(root / "experiments/falsification_report.json")
    evidence_audit = write_evidence_audit(root / "experiments/evidence_audit.json", root)
    submission = write_submission_dossier(root / "experiments/submission_dossier.json", root)
    payload: dict[str, Any] = {
        "run_id": "TW-REPRODUCTION-v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "runtime": {"python": sys.version.split()[0], "platform": platform.platform()},
        "steps": results,
        "capsule_fingerprint": capsule["fingerprint"],
        "packet_fingerprint": packet["fingerprint"],
        "falsification_fingerprint": falsification["fingerprint"],
        "evidence_audit_fingerprint": evidence_audit["fingerprint"],
        "submission_fingerprint": submission["fingerprint"],
        "status": "reproducible_review_ready" if all(step["passed"] for step in results) else "reproduction_failed",
        "claim_boundary": "This verifies software and public computational evidence; it does not create wet-lab measurements or establish a biological discovery.",
    }
    payload["fingerprint"] = fingerprint(payload)
    target = root / output
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
