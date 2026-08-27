from __future__ import annotations

import hashlib
import json
import platform
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ReleaseCheck:
    name: str
    passed: bool
    severity: str
    detail: str
    evidence: list[str]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(root: Path, relative: str) -> Any | None:
    path = root / relative
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _check_artifact(root: Path, relative: str, required: bool = True) -> ReleaseCheck:
    path = root / relative
    exists = path.is_file() and path.stat().st_size > 0
    return ReleaseCheck(
        name=f"artifact:{relative}",
        passed=exists or not required,
        severity="blocker" if required else "warning",
        detail="present and non-empty" if exists else ("optional artifact absent" if not required else "required artifact is absent or empty"),
        evidence=[relative] if exists else [],
    )


def build_release_report(root: str | Path = ".") -> dict[str, Any]:
    base = Path(root).resolve()
    checks: list[ReleaseCheck] = []
    for relative in (
        "README.md",
        "experiments/reproduction_run.json",
        "experiments/evidence_audit.json",
        "experiments/research_capsule.json",
        "experiments/review_packet.json",
        "experiments/falsification_report.json",
        "experiments/wetlab_validation_protocol.json",
    ):
        checks.append(_check_artifact(base, relative))

    capsule = _load_json(base, "experiments/research_capsule.json")
    if isinstance(capsule, dict):
        status = capsule.get("status")
        checks.append(ReleaseCheck(
            name="claim_status_is_bounded",
            passed=status in {"review_ready_not_discovery", "computationally_reproducible_evidence_bounded", "not_established"},
            severity="blocker",
            detail=f"capsule status={status!r}; discovery claims remain gated",
            evidence=["experiments/research_capsule.json"],
        ))
        checks.append(ReleaseCheck(
            name="capsule_fingerprint_present",
            passed=bool(capsule.get("fingerprint")),
            severity="blocker",
            detail="release capsule carries an integrity fingerprint" if capsule.get("fingerprint") else "missing capsule fingerprint",
            evidence=["experiments/research_capsule.json"],
        ))
    else:
        checks.append(ReleaseCheck("claim_status_is_bounded", False, "blocker", "capsule is unavailable or malformed", []))

    audit = _load_json(base, "experiments/evidence_audit.json")
    audit_passed = isinstance(audit, dict) and audit.get("status") in {"computationally_reproducible_evidence_bounded", "review_ready_not_discovery"}
    checks.append(ReleaseCheck(
        name="evidence_audit_passes",
        passed=audit_passed,
        severity="blocker",
        detail=f"evidence audit status={audit.get('status')!r}" if isinstance(audit, dict) else "evidence audit unavailable",
        evidence=["experiments/evidence_audit.json"] if isinstance(audit, dict) else [],
    ))

    code_files = sorted((base / "thermodynamic_waddington").rglob("*.py"))
    code_hash = hashlib.sha256()
    for path in code_files:
        code_hash.update(str(path.relative_to(base)).encode())
        code_hash.update(_sha256(path).encode())
    passed = sum(check.passed for check in checks)
    blockers = [check.detail for check in checks if not check.passed and check.severity == "blocker"]
    payload = {
        "schema": "thermodynamic-waddington/evidence-release-v1",
        "status": "release_ready_evidence_bounded" if not blockers else "blocked",
        "passed_checks": passed,
        "total_checks": len(checks),
        "blockers": blockers,
        "checks": [check.to_dict() for check in checks],
        "code": {"python_files": len(code_files), "sha256": code_hash.hexdigest()},
        "runtime": {"python": sys.version.split()[0], "platform": platform.platform()},
        "scientific_boundary": "The software can release reproducible computational evidence; it cannot convert public or simulated data into a causal biological discovery without measured intervention replication.",
    }
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    payload["fingerprint"] = fingerprint
    return payload


def write_release_report(path: str | Path, root: str | Path = ".") -> dict[str, Any]:
    report = build_release_report(root)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report
