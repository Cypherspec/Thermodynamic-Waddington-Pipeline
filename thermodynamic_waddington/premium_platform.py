from __future__ import annotations

"""Premium, evidence-gated orchestration for Thermodynamic Waddington.

This module turns the research package into a reproducible platform surface:
one command refreshes the computational evidence bundle, checks provenance,
and emits a machine-readable release decision. It never upgrades an
observational result into a causal biological claim.
"""

import hashlib
import json
import platform
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

from .award_capsule import build_capsule
from .evidence_audit import build_evidence_audit
from .evidence_graph import build_evidence_graph
from .evidence_release import build_release_report
from .falsification import run_falsification_suite
from .reproduce import run_reproduction
from .review_packet import build_review_packet
from .robustness import run_synthetic_robustness


@dataclass(frozen=True)
class PlatformRun:
    schema: str
    status: str
    release_status: str
    claim_level: str
    generated_artifacts: list[str]
    checks: dict[str, Any]
    blockers: list[str]
    warnings: list[str]
    fingerprint: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _write_json(root: Path, relative: str, payload: Any) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(payload, "to_dict"):
        payload = payload.to_dict()
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def _status(release: dict[str, Any], graph: dict[str, Any]) -> tuple[str, str, list[str]]:
    blockers = list(release.get("blockers", [])) + list(graph.get("blockers", []))
    blockers = list(dict.fromkeys(blockers))
    if blockers:
        return "evidence_bounded_blocked_for_causal_claims", "computationally_reproducible", blockers
    return "premium_release_ready_evidence_bounded", "computationally_reproducible", []


def run_premium_platform(root: str | Path = ".", refresh: bool = True) -> dict[str, Any]:
    base = Path(root).resolve()
    generated: list[str] = []
    if refresh:
        reproduction = run_reproduction(base, base / "experiments/reproduction_run.json")
        generated.append("experiments/reproduction_run.json")
        _write_json(base, "experiments/falsification_report.json", run_falsification_suite())
        generated.append("experiments/falsification_report.json")
        _write_json(base, "experiments/robustness_report.json", run_synthetic_robustness())
        generated.append("experiments/robustness_report.json")
    else:
        reproduction = None

    audit = build_evidence_audit(base)
    _write_json(base, "experiments/evidence_audit.json", audit)
    generated.append("experiments/evidence_audit.json")
    capsule = build_capsule(base).to_dict()
    _write_json(base, "experiments/research_capsule.json", capsule)
    generated.append("experiments/research_capsule.json")
    packet = build_review_packet(base)
    _write_json(base, "experiments/review_packet.json", packet)
    generated.append("experiments/review_packet.json")
    graph = build_evidence_graph(base).to_dict()
    _write_json(base, "experiments/evidence_graph.json", graph)
    generated.append("experiments/evidence_graph.json")
    release = build_release_report(base)
    release_status, claim_level, blockers = _status(release, graph)
    warnings = list(graph.get("warnings", []))
    checks = {
        "release": release,
        "audit": audit,
        "capsule": {"status": capsule.get("status"), "fingerprint": capsule.get("fingerprint")},
        "review_packet": {"status": packet.get("status"), "fingerprint": packet.get("fingerprint")},
        "evidence_graph": {"status": graph.get("status"), "nodes": graph.get("counts", {}).get("nodes"), "edges": graph.get("counts", {}).get("edges")},
        "reproduction": {"status": reproduction.get("status")} if isinstance(reproduction, dict) else {"status": "not_refreshed"},
    }
    payload_without_fp = {
        "schema": "thermodynamic-waddington/premium-platform-v1",
        "status": release_status,
        "release_status": release.get("status"),
        "claim_level": claim_level,
        "generated_artifacts": list(dict.fromkeys(generated)),
        "checks": checks,
        "blockers": blockers,
        "warnings": warnings,
        "runtime": {"python": sys.version.split()[0], "platform": platform.platform()},
        "scientific_boundary": "Premium orchestration improves reproducibility, provenance, and reviewability; it does not establish causality without measured randomized intervention data, independent donors, lineage/barcode data, viability, orthogonal assays, and preregistered endpoints.",
    }
    fingerprint = hashlib.sha256(json.dumps(payload_without_fp, sort_keys=True, default=str).encode()).hexdigest()
    result = {**payload_without_fp, "fingerprint": fingerprint}
    _write_json(base, "experiments/premium_platform_run.json", result)
    return result


def write_premium_report(path: str | Path = "experiments/premium_platform_run.json", root: str | Path = ".") -> dict[str, Any]:
    result = run_premium_platform(root)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return result
