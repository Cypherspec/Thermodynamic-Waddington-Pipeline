from __future__ import annotations

"""Formal claim ladder for reproducible biology and discovery review.

No model output can skip levels. The gate is intentionally stricter than a
benchmark leaderboard because publication, causal, and translational claims
need different evidence.
"""

import hashlib
import json
from dataclasses import asdict, dataclass, field
from typing import Any, Mapping, Sequence

LEVELS = (
    "computational_pattern",
    "heldout_predictive_result",
    "replicated_predictive_result",
    "experimentally_supported_candidate",
    "independently_replicated_effect",
    "mechanistically_supported_effect",
    "clinical_or_therapeutic_relevance_not_established",
)


@dataclass(frozen=True)
class EvidenceItem:
    name: str
    present: bool
    provenance: str
    notes: str = ""


@dataclass
class ClaimDecision:
    attained_level: str
    next_level: str | None
    accepted: bool
    blockers: list[str] = field(default_factory=list)
    evidence: list[dict[str, Any]] = field(default_factory=list)
    claim_text: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self) | {"fingerprint": fingerprint(asdict(self))}


def fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def decide(evidence: Sequence[EvidenceItem], requested: str = "experimentally_supported_candidate") -> ClaimDecision:
    if requested not in LEVELS:
        raise ValueError(f"unknown claim level: {requested}")
    present = {item.name for item in evidence if item.present}
    blockers: list[str] = []
    rules = {
        "heldout_predictive_result": ("heldout_outcome", "no held-out outcome benchmark"),
        "replicated_predictive_result": ("independent_dataset", "no independent dataset replication"),
        "experimentally_supported_candidate": ("randomized_intervention", "no randomized measured intervention"),
        "independently_replicated_effect": ("independent_biological_replicate", "no independent biological replication"),
        "mechanistically_supported_effect": ("orthogonal_mechanism", "no orthogonal mechanistic assay"),
    }
    attained = "computational_pattern"
    for level in LEVELS[1:]:
        key, message = rules[level]
        if key in present:
            attained = level
        else:
            blockers.append(message)
            break
    if requested == "clinical_or_therapeutic_relevance_not_established":
        blockers.append("clinical or therapeutic relevance is never established by this software alone")
    accepted = LEVELS.index(attained) >= LEVELS.index(requested)
    next_level = None if attained == LEVELS[-1] else LEVELS[LEVELS.index(attained) + 1]
    return ClaimDecision(attained, next_level, accepted, blockers, [item.to_dict() if hasattr(item, "to_dict") else asdict(item) for item in evidence], f"Evidence supports {attained}; it does not support any stronger claim unless all gates are cleared.")
