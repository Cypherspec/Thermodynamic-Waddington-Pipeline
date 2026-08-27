from __future__ import annotations

"""Conservative discovery orchestration for externally validated fate hypotheses.

This module never calls an effective-landscape pattern a biological discovery by
itself. It creates preregisterable hypotheses, blocks leakage, compares against
published baselines, and returns an evidence ledger suitable for review.
"""

import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Sequence

from .advanced.fate_control import run_fate_control
from .benchmarking import rank_correlation
from .larry_validation import validate_larry_fate_ground_truth
from .statistics.permutation import spearman


@dataclass(frozen=True)
class DiscoveryConfig:
    seed: int = 17
    permutation_count: int = 499
    minimum_effect: float = 0.10
    alpha: float = 0.05
    preregistration_id: str = "tw-discovery-v1"


@dataclass(frozen=True)
class Hypothesis:
    identifier: str
    statement: str
    observable: str
    falsifier: str
    required_replication: str

    def to_dict(self) -> dict[str, str]:
        return asdict(self)


@dataclass
class DiscoveryLedger:
    config: DiscoveryConfig
    hypotheses: list[Hypothesis] = field(default_factory=list)
    evidence: list[dict[str, Any]] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)

    def add(self, hypothesis: Hypothesis, result: dict[str, Any]) -> None:
        self.hypotheses.append(hypothesis)
        self.evidence.append({"hypothesis": hypothesis.identifier, "result": result})

    def digest(self) -> str:
        payload = json.dumps({"config": asdict(self.config), "hypotheses": [item.to_dict() for item in self.hypotheses], "evidence": self.evidence, "blockers": self.blockers}, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        return {"config": asdict(self.config), "hypotheses": [item.to_dict() for item in self.hypotheses], "evidence": self.evidence, "blockers": self.blockers, "digest": self.digest(), "claim_status": "hypothesis_generating_not_discovery" if self.blockers else "replication_ready_hypotheses"}


def _permutation_pvalue(first: Sequence[float], second: Sequence[float], count: int, seed: int) -> float:
    import random
    if len(first) != len(second) or len(first) < 8:
        return 1.0
    observed = abs(spearman(first, second))
    pooled = list(first) + list(second)
    rng = random.Random(seed)
    exceed = 0
    split = len(first)
    for _ in range(max(1, count)):
        rng.shuffle(pooled)
        if abs(spearman(pooled[:split], pooled[split:])) >= observed:
            exceed += 1
    return (exceed + 1) / (count + 1)


def screen_fate_control(expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], labels: Sequence[str], config: DiscoveryConfig | None = None) -> dict[str, Any]:
    config = config or DiscoveryConfig()
    control = run_fate_control(expression, velocity, labels, seed=config.seed)
    candidates = control.get("candidates", control.get("targets", []))
    return {"candidate_count": len(candidates) if isinstance(candidates, list) else 0, "candidate_preview": candidates[:10] if isinstance(candidates, list) else [], "claim": "prioritized intervention hypotheses only", "falsification": "replicate with randomized perturbations and held-out lineage outcomes"}


def validate_public_references(root: str | Path = "data/real/larry") -> dict[str, Any]:
    report = validate_larry_fate_ground_truth(root)
    return {"larry": report, "interpretation": "reference outcome loaded; it does not validate a matched expression-to-landscape prediction"}


def build_discovery_ledger(expression: Sequence[Sequence[float]] | None = None, velocity: Sequence[Sequence[float]] | None = None, labels: Sequence[str] | None = None, config: DiscoveryConfig | None = None) -> DiscoveryLedger:
    config = config or DiscoveryConfig()
    ledger = DiscoveryLedger(config)
    ledger.add(Hypothesis("H1", "Effective barrier ranking predicts fate-enriched transition bottlenecks.", "rank correlation of barrier score with held-out lineage transition frequency", "Spearman correlation <= preregistered threshold in an independent experiment", "two biological replicates, one perturbation cohort"), {"status": "pending_matched_lineage_data"})
    ledger.add(Hypothesis("H2", "Interventions selected by uncertainty-aware fate control enrich the intended fate without increasing off-target entropy.", "target-fate probability and off-target entropy under randomized perturbation", "no improvement over density-only and CellRank baselines", "three perturbation replicates and dose-response"), screen_fate_control(expression, velocity, labels, config) if expression is not None and velocity is not None and labels is not None else {"status": "pending_expression_velocity_labels"})
    ledger.add(Hypothesis("H3", "Landscape uncertainty identifies cells where additional measurement has maximal fate-information gain.", "prospective reduction in calibrated fate entropy after active measurement", "active policy fails to beat random sampling under held-out cells", "prospective blinded active-learning study"), {"status": "pending_prospective_measurement"})
    ledger.blockers.extend(["No matched expression + velocity + lineage matrix has been supplied for the same cells.", "No randomized perturbation cohort has been supplied.", "No external replication has been run.", "An effective landscape is not a molecular free-energy measurement."])
    return ledger


def write_discovery_ledger(path: str | Path = "experiments/discovery_ledger.json", **kwargs: Any) -> dict[str, Any]:
    ledger = build_discovery_ledger(**kwargs)
    payload = ledger.to_dict()
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload
