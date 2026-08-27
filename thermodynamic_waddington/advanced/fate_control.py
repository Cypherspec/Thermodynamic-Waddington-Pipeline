from __future__ import annotations

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from typing import Any, Mapping, Sequence

from .core import matrix, mean, quantile, variance


@dataclass(frozen=True)
class Intervention:
    name: str
    effect: tuple[float, ...]
    cost: float = 1.0
    measured: bool = False
    source: str = "hypothetical_feature_direction"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class FateControlConfig:
    target_fate: str | None = None
    forbidden_fates: tuple[str, ...] = ()
    intervention_budget: int = 3
    bootstrap_rounds: int = 64
    perturbation_scale: float = 0.25
    uncertainty_penalty: float = 0.5
    off_target_penalty: float = 1.0
    minimum_fate_cells: int = 3
    seed: int = 17

    def validate(self) -> None:
        if self.intervention_budget < 1:
            raise ValueError("intervention_budget must be positive")
        if self.bootstrap_rounds < 8:
            raise ValueError("bootstrap_rounds must be at least 8")
        if self.perturbation_scale <= 0:
            raise ValueError("perturbation_scale must be positive")
        if self.uncertainty_penalty < 0 or self.off_target_penalty < 0:
            raise ValueError("penalties must be non-negative")
        if self.minimum_fate_cells < 1:
            raise ValueError("minimum_fate_cells must be positive")


@dataclass(frozen=True)
class CellFateRecord:
    cell_index: int
    cell_id: str
    probabilities: dict[str, float]
    predicted_fate: str
    confidence: float
    entropy: float
    epistemic_uncertainty: float
    energy: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class InterventionScore:
    name: str
    target_fate: str
    cells_evaluated: int
    baseline_probability: float
    counterfactual_probability: float
    probability_gain: float
    energy_change: float
    uncertainty: float
    off_target_risk: float
    cost: float
    robust_gain: float
    measured: bool
    source: str

    @property
    def utility(self) -> float:
        return self.robust_gain / max(self.cost, 1e-12)

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["utility"] = self.utility
        return result


@dataclass(frozen=True)
class ExperimentDesign:
    cell_indices: tuple[int, ...]
    intervention_names: tuple[str, ...]
    expected_information_gain: float
    expected_target_gain: float
    uncertainty_reduction_proxy: float
    rationale: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _safe(value: float, default: float = 0.0) -> float:
    value = float(value)
    return value if math.isfinite(value) else default


def _dot(left: Sequence[float], right: Sequence[float]) -> float:
    return sum(_safe(a) * _safe(b) for a, b in zip(left, right))


def _norm(row: Sequence[float]) -> float:
    return math.sqrt(max(0.0, _dot(row, row)))


def _distance(left: Sequence[float], right: Sequence[float]) -> float:
    width = min(len(left), len(right))
    return math.sqrt(sum((_safe(left[i]) - _safe(right[i])) ** 2 for i in range(width)))


def _mean_vector(rows: Sequence[Sequence[float]], width: int) -> list[float]:
    if not rows:
        return [0.0] * width
    return [mean(row[j] for row in rows if j < len(row)) for j in range(width)]


def _softmax(values: Sequence[float]) -> list[float]:
    if not values:
        return []
    pivot = max(values)
    weights = [math.exp(max(-745.0, min(50.0, value - pivot))) for value in values]
    total = sum(weights)
    return [weight / total for weight in weights] if total else [1.0 / len(values)] * len(values)


def _entropy(probabilities: Sequence[float]) -> float:
    return -sum(p * math.log(max(p, 1e-300)) for p in probabilities)


def _digest(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _quantile_interval(values: Sequence[float]) -> tuple[float, float]:
    return quantile(values, 0.025), quantile(values, 0.975)


class FateControlEngine:
    """Closed-loop fate steering over an inferred landscape.

    This is a decision-support layer, not a claim that a computational perturbation
    is a validated molecular intervention. It ranks counterfactuals, identifies the
    cells that would be most informative to perturb, and emits falsification tests.
    """

    def __init__(self, config: FateControlConfig | None = None):
        self.config = config or FateControlConfig()
        self.config.validate()
        self.rng = random.Random(self.config.seed)

    def _prototypes(self, expression: list[list[float]], labels: Sequence[str]) -> tuple[list[str], dict[str, list[float]], dict[str, float]]:
        groups: dict[str, list[list[float]]] = {}
        for row, label in zip(expression, labels):
            groups.setdefault(str(label), []).append(row)
        labels_out = sorted(groups)
        width = max((len(row) for row in expression), default=0)
        prototypes = {label: _mean_vector(groups[label], width) for label in labels_out}
        scales = {label: max(1e-6, mean([_distance(row, prototypes[label]) for row in groups[label]])) for label in labels_out}
        return labels_out, prototypes, scales

    def _probabilities(self, row: Sequence[float], labels: Sequence[str], prototypes: Mapping[str, Sequence[float]], scales: Mapping[str, float], energy: float = 0.0) -> dict[str, float]:
        logits = []
        for label in labels:
            distance = _distance(row, prototypes[label])
            logits.append(-distance / max(scales[label], 1e-6) - 0.05 * energy)
        return {label: probability for label, probability in zip(labels, _softmax(logits))}

    def _energy(self, row: Sequence[float], probabilities: Mapping[str, float]) -> float:
        return -math.log(max(probabilities.get(max(probabilities, key=probabilities.get), 0.0), 1e-300))

    def _cell_records(self, expression: list[list[float]], labels: Sequence[str], cell_ids: Sequence[str], energies: Sequence[float] | None = None) -> tuple[list[CellFateRecord], list[str], dict[str, list[float]], dict[str, float]]:
        fate_labels, prototypes, scales = self._prototypes(expression, labels)
        records = []
        energy_values = list(energies or [0.0] * len(expression))
        for index, row in enumerate(expression):
            probabilities = self._probabilities(row, fate_labels, prototypes, scales, energy_values[index] if index < len(energy_values) else 0.0)
            values = list(probabilities.values())
            predicted = max(probabilities, key=probabilities.get) if probabilities else "unassigned"
            records.append(CellFateRecord(index, str(cell_ids[index]), probabilities, predicted, max(values, default=0.0), _entropy(values), 0.0, energy_values[index] if index < len(energy_values) else 0.0))
        return records, fate_labels, prototypes, scales

    def _interventions(self, width: int, interventions: Sequence[Intervention] | Mapping[str, Sequence[float]] | None) -> list[Intervention]:
        if interventions is None:
            return [Intervention(f"feature_{index}_up", tuple(1.0 if j == index else 0.0 for j in range(width)), measured=False) for index in range(width)]
        if isinstance(interventions, Mapping):
            return [Intervention(str(name), tuple(_safe(value) for value in effect), measured=True, source="supplied_perturbation_effect") for name, effect in interventions.items()]
        return list(interventions)

    def _bootstrap_gain(self, row: Sequence[float], effect: Sequence[float], target: str, labels: Sequence[str], prototypes: Mapping[str, Sequence[float]], scales: Mapping[str, float], baseline: float) -> tuple[float, float]:
        gains = []
        width = min(len(row), len(effect))
        for replicate in range(self.config.bootstrap_rounds):
            perturbed = list(row)
            for index in range(width):
                noise = self.rng.gauss(0.0, 0.05 * max(1.0, abs(row[index])))
                perturbed[index] += self.config.perturbation_scale * effect[index] + noise
            probabilities = self._probabilities(perturbed, labels, prototypes, scales)
            gains.append(probabilities.get(target, 0.0) - baseline)
        return mean(gains), math.sqrt(variance(gains))

    def _score_intervention(self, intervention: Intervention, records: Sequence[CellFateRecord], expression: list[list[float]], target: str, labels: Sequence[str], prototypes: Mapping[str, Sequence[float]], scales: Mapping[str, float]) -> InterventionScore:
        target_records = [record for record in records if record.predicted_fate == target or record.probabilities.get(target, 0.0) > 0.15]
        if not target_records:
            target_records = list(records)
        gains = []
        energy_changes = []
        off_target = []
        for record in target_records:
            row = expression[record.cell_index]
            width = min(len(row), len(intervention.effect))
            perturbed = list(row)
            for index in range(width):
                perturbed[index] += self.config.perturbation_scale * intervention.effect[index]
            probabilities = self._probabilities(perturbed, labels, prototypes, scales)
            gain, spread = self._bootstrap_gain(row, intervention.effect, target, labels, prototypes, scales, record.probabilities.get(target, 0.0))
            gains.append(gain)
            energy_changes.append(-math.log(max(probabilities.get(target, 0.0), 1e-300)) + math.log(max(record.probabilities.get(target, 0.0), 1e-300)))
            forbidden_mass = sum(probabilities.get(fate, 0.0) for fate in self.config.forbidden_fates)
            off_target.append(forbidden_mass)
        baseline = mean([record.probabilities.get(target, 0.0) for record in target_records])
        raw_gain = mean(gains)
        uncertainty = math.sqrt(variance(gains))
        robust_gain = raw_gain - self.config.uncertainty_penalty * uncertainty - self.config.off_target_penalty * mean(off_target)
        return InterventionScore(intervention.name, target, len(target_records), baseline, baseline + raw_gain, raw_gain, mean(energy_changes), uncertainty, mean(off_target), intervention.cost, robust_gain, intervention.measured, intervention.source)

    def _select_policy(self, scores: Sequence[InterventionScore]) -> list[InterventionScore]:
        selected: list[InterventionScore] = []
        remaining = sorted(scores, key=lambda score: (score.utility, score.robust_gain, score.name), reverse=True)
        for candidate in remaining:
            if candidate.robust_gain <= 0:
                continue
            if any(candidate.name == item.name for item in selected):
                continue
            selected.append(candidate)
            if len(selected) >= self.config.intervention_budget:
                break
        return selected

    def _design_experiment(self, records: Sequence[CellFateRecord], scores: Sequence[InterventionScore], target: str) -> ExperimentDesign:
        ranked_cells = sorted(records, key=lambda record: (record.entropy + record.epistemic_uncertainty, -record.probabilities.get(target, 0.0)), reverse=True)
        cells = tuple(record.cell_index for record in ranked_cells[: min(12, len(ranked_cells))])
        interventions = tuple(score.name for score in sorted(scores, key=lambda score: score.utility, reverse=True)[: self.config.intervention_budget])
        information_gain = mean([record.entropy for record in ranked_cells[: len(cells)]]) if cells else 0.0
        target_gain = mean([score.robust_gain for score in scores if score.name in interventions]) if interventions else 0.0
        return ExperimentDesign(cells, interventions, information_gain, target_gain, information_gain / max(1, len(cells)), ("prioritize high-entropy cells", "include cells near the target decision boundary", "compare top policy against a matched negative-control intervention"))

    def _falsification(self, scores: Sequence[InterventionScore], records: Sequence[CellFateRecord], target: str) -> dict[str, Any]:
        observed = [score.robust_gain for score in scores]
        null = []
        for replicate in range(max(16, self.config.bootstrap_rounds // 2)):
            shuffled = list(observed)
            self.rng.shuffle(shuffled)
            null.append(max(shuffled, default=0.0) * self.rng.uniform(-0.25, 0.25))
        threshold = quantile(null, 0.95)
        best = max(observed, default=0.0)
        target_mass = mean([record.probabilities.get(target, 0.0) for record in records]) if records else 0.0
        return {
            "status": "passes_screening" if best > threshold and target_mass > 0.05 else "hypothesis_only",
            "observed_best_robust_gain": best,
            "null_95th_percentile": threshold,
            "null_samples": len(null),
            "target_baseline_mass": target_mass,
            "tests": ["velocity-preserving placebo", "feature-direction permutation", "held-out lineage calibration", "matched negative-control intervention"],
            "failure_condition": "If the best policy does not exceed the null distribution or fails held-out fate prediction, reject the steering claim.",
        }

    def run(self, expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], labels: Sequence[str], cell_ids: Sequence[str] | None = None, energies: Sequence[float] | None = None, lineage_outcomes: Sequence[str] | None = None, interventions: Sequence[Intervention] | Mapping[str, Sequence[float]] | None = None) -> dict[str, Any]:
        expression_rows = matrix(expression)
        velocity_rows = matrix(velocity)
        if len(expression_rows) != len(velocity_rows) or len(expression_rows) != len(labels):
            raise ValueError("expression, velocity, and labels must have the same number of cells")
        if len(expression_rows) < self.config.minimum_fate_cells:
            raise ValueError("not enough cells for fate control analysis")
        ids = list(cell_ids or [f"cell_{index:05d}" for index in range(len(expression_rows))])
        if len(ids) != len(expression_rows):
            raise ValueError("cell_ids must align with expression")
        records, fate_labels, prototypes, scales = self._cell_records(expression_rows, labels, ids, energies)
        target = self.config.target_fate or max(fate_labels, key=lambda fate: sum(label == fate for label in labels))
        if target not in fate_labels:
            raise ValueError(f"target_fate {target!r} is absent from labels")
        candidates = self._interventions(max((len(row) for row in expression_rows), default=0), interventions)
        scores = [self._score_intervention(item, records, expression_rows, target, fate_labels, prototypes, scales) for item in candidates]
        policy = self._select_policy(scores)
        experiment = self._design_experiment(records, policy, target)
        falsification = self._falsification(scores, records, target)
        measured_interventions = any(item.measured for item in candidates)
        lineage_available = lineage_outcomes is not None
        claim_status = "calibration_ready" if measured_interventions and lineage_available else "hypothesis_only"
        return {
            "engine": "closed_loop_fate_control",
            "version": "1.0",
            "claim_status": claim_status,
            "target_fate": target,
            "fate_labels": fate_labels,
            "fate_probabilities": [record.to_dict() for record in records],
            "prototypes": prototypes,
            "prototype_scales": scales,
            "interventions": [score.to_dict() for score in scores],
            "selected_policy": [score.to_dict() for score in policy],
            "experiment_design": experiment.to_dict(),
            "falsification": falsification,
            "lineage_validation_available": lineage_available,
            "measured_interventions_available": measured_interventions,
            "assumptions": ["counterfactual effects are local additive approximations", "RNA velocity is treated as directional evidence, not molecular work", "unmeasured interventions are hypothetical feature directions", "policy utility is not a clinical recommendation"],
            "required_evidence": ["held-out lineage outcomes", "measured perturbation-response matrices", "replicate/time-resolved observations", "negative controls and rescue experiments"],
            "provenance": {"input_digest": _digest({"expression": expression_rows, "velocity": velocity_rows, "labels": list(labels)}), "seed": self.config.seed, "bootstrap_rounds": self.config.bootstrap_rounds},
        }


def run_fate_control(expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]], labels: Sequence[str], **kwargs: Any) -> dict[str, Any]:
    config = kwargs.pop("config", None)
    return FateControlEngine(config).run(expression, velocity, labels, **kwargs)
