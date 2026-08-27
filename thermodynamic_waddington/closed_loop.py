from __future__ import annotations

"""Closed-loop experimental design for cell-fate hypotheses.

The engine is deliberately asymmetric about evidence: computational predictions
can propose experiments, but only measured perturbation outcomes can update a
claim from ``hypothesis`` to ``experimentally_supported``.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class InterventionCandidate:
    name: str
    effect: tuple[float, ...]
    cost: float = 1.0
    feasibility: float = 1.0
    modality: str = "unmeasured_counterfactual"
    negative_control: bool = False
    safety_notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MeasuredOutcome:
    intervention: str
    replicate: str
    donor: str
    cell_count: int
    target_probability: float
    control_probability: float
    observed: bool = True
    batch: str = "unspecified"
    timepoint: str = "unspecified"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ClosedLoopConfig:
    target_fate: str
    temperature: float = 1.0
    exploration_weight: float = 0.5
    uncertainty_weight: float = 0.75
    information_weight: float = 0.5
    cost_weight: float = 0.25
    minimum_biological_replicates: int = 3
    maximum_interventions_per_round: int = 4
    seed: int = 17

    def validate(self) -> None:
        if not self.target_fate:
            raise ValueError("target_fate is required")
        if self.temperature <= 0:
            raise ValueError("temperature must be positive")
        if self.minimum_biological_replicates < 2:
            raise ValueError("at least two biological replicates are required")
        if self.maximum_interventions_per_round < 1:
            raise ValueError("maximum_interventions_per_round must be positive")
        for value in (self.exploration_weight, self.uncertainty_weight, self.information_weight, self.cost_weight):
            if value < 0:
                raise ValueError("utility weights must be non-negative")


@dataclass(frozen=True)
class PolicyScore:
    intervention: str
    baseline_target_mass: float
    predicted_target_gain: float
    epistemic_uncertainty: float
    expected_information_gain: float
    cost: float
    feasibility: float
    utility: float
    status: str
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ClosedLoopState:
    round_index: int
    config: ClosedLoopConfig
    scores: list[PolicyScore] = field(default_factory=list)
    selected: list[str] = field(default_factory=list)
    observations: list[MeasuredOutcome] = field(default_factory=list)
    posterior: dict[str, dict[str, float]] = field(default_factory=dict)
    claim_status: str = "hypothesis_only"
    blockers: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "round_index": self.round_index,
            "config": asdict(self.config),
            "scores": [score.to_dict() for score in self.scores],
            "selected": list(self.selected),
            "observations": [observation.to_dict() for observation in self.observations],
            "posterior": self.posterior,
            "claim_status": self.claim_status,
            "blockers": list(self.blockers),
            "fingerprint": fingerprint(self),
        }


class ClosedLoopEngine:
    def __init__(self, config: ClosedLoopConfig):
        config.validate()
        self.config = config
        self.rng = random.Random(config.seed)
        self.state = ClosedLoopState(0, config)

    def score_candidates(
        self,
        expression: Sequence[Sequence[float]],
        velocity: Sequence[Sequence[float]],
        labels: Sequence[str],
        energies: Sequence[float],
        candidates: Sequence[InterventionCandidate],
    ) -> list[PolicyScore]:
        if not expression or len(expression) != len(velocity) or len(expression) != len(labels):
            raise ValueError("expression, velocity, and labels must be nonempty and cell-aligned")
        target_cells = [i for i, label in enumerate(labels) if str(label) == self.config.target_fate]
        baseline = len(target_cells) / max(1, len(labels))
        scores: list[PolicyScore] = []
        for candidate in candidates:
            if candidate.cost <= 0 or candidate.feasibility <= 0:
                continue
            gains: list[float] = []
            for row, flow, energy in zip(expression, velocity, energies):
                width = min(len(row), len(candidate.effect), len(flow))
                perturbation = sum(float(candidate.effect[j]) * (1.0 + float(flow[j])) for j in range(width))
                fate_alignment = 1.0 if str(labels[len(gains)]) == self.config.target_fate else -0.35
                gain = math.tanh(perturbation / max(self.config.temperature, 1e-12)) * fate_alignment
                gain *= 1.0 / (1.0 + abs(float(energy)))
                gains.append(gain)
            mean_gain = sum(gains) / max(1, len(gains))
            variance = sum((value - mean_gain) ** 2 for value in gains) / max(1, len(gains) - 1)
            uncertainty = math.sqrt(variance)
            information = uncertainty * math.log1p(len(expression))
            warnings = list(candidate.safety_notes)
            if candidate.modality != "measured_intervention":
                warnings.append("counterfactual_only_not_measured")
            if candidate.negative_control:
                warnings.append("negative_control")
            utility = (
                mean_gain
                + self.config.exploration_weight * information
                + self.config.uncertainty_weight * uncertainty
                + self.config.information_weight * information
                - self.config.cost_weight * candidate.cost
            ) * candidate.feasibility
            scores.append(PolicyScore(candidate.name, baseline, mean_gain, uncertainty, information, candidate.cost, candidate.feasibility, utility, "candidate", tuple(warnings)))
        return sorted(scores, key=lambda item: (-item.utility, item.intervention))

    def select(self, scores: Sequence[PolicyScore]) -> list[str]:
        selected: list[str] = []
        for score in scores:
            if score.intervention in selected or score.utility <= 0:
                continue
            selected.append(score.intervention)
            if len(selected) >= self.config.maximum_interventions_per_round:
                break
        self.state.selected = selected
        self.state.scores = list(scores)
        return selected

    def update(self, observations: Sequence[MeasuredOutcome]) -> ClosedLoopState:
        for observation in observations:
            if not observation.observed:
                continue
            if observation.cell_count < 1:
                raise ValueError("measured outcomes require a positive cell_count")
            if not 0.0 <= observation.target_probability <= 1.0 or not 0.0 <= observation.control_probability <= 1.0:
                raise ValueError("measured probabilities must lie in [0, 1]")
            self.state.observations.append(observation)
            posterior = self.state.posterior.setdefault(observation.intervention, {"alpha": 1.0, "beta": 1.0, "replicates": 0.0})
            observed_target = float(observation.target_probability)
            observed_control = float(observation.control_probability)
            posterior["alpha"] += observed_target
            posterior["beta"] += max(0.0, 1.0 - observed_target)
            posterior["replicates"] += 1.0
            posterior["control_mean"] = ((posterior.get("control_mean", 0.0) * (posterior["replicates"] - 1.0)) + observed_control) / posterior["replicates"]
            posterior["effect_mean"] = posterior["alpha"] / max(1e-12, posterior["alpha"] + posterior["beta"]) - posterior["control_mean"]
        measured_by_intervention: dict[str, list[MeasuredOutcome]] = {}
        for observation in self.state.observations:
            measured_by_intervention.setdefault(observation.intervention, []).append(observation)
        replicated = [items for items in measured_by_intervention.values() if len({item.donor for item in items}) >= self.config.minimum_biological_replicates]
        self.state.claim_status = "experimentally_supported_candidate" if replicated else "hypothesis_only"
        self.state.blockers = [] if replicated else ["insufficient independent biological replicates", "no causal claim without randomized perturbation and matched controls"]
        self.state.round_index += 1
        return self.state

    def design_round(self, expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], labels: Sequence[str], energies: Sequence[float], candidates: Sequence[InterventionCandidate]) -> ClosedLoopState:
        scores = self.score_candidates(expression, velocity, labels, energies, candidates)
        self.select(scores)
        self.state.claim_status = "ready_for_wetlab_validation"
        self.state.blockers = ["no measured intervention outcomes have been supplied"]
        return self.state


def fingerprint(value: Any) -> str:
    if isinstance(value, ClosedLoopState):
        payload = {"round_index": value.round_index, "selected": value.selected, "scores": [score.to_dict() for score in value.scores], "observations": [item.to_dict() for item in value.observations], "posterior": value.posterior, "claim_status": value.claim_status}
    else:
        payload = value
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def write_round_plan(state: ClosedLoopState, path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(state.to_dict(), indent=2, sort_keys=True), encoding="utf-8")
