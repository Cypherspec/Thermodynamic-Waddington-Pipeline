from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any, Sequence

from .arrays import mean, variance
from .cell_atlas import CellInterpretation, ProgramScore, interpret_cell


@dataclass(frozen=True)
class Prediction:
    key: str
    title: str
    state: str
    probability: float
    confidence: str
    evidence: tuple[str, ...]
    caveats: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CellPrediction:
    cell_index: int
    cell_id: str
    predicted_identity: str
    identity_probability: float
    identity_confidence: str
    developmental_stage: str
    stage_probability: float
    damage_state: str
    damage_probability: float
    viability_state: str
    viability_probability: float
    cell_cycle_state: str
    dominant_function: str
    predictions: tuple[Prediction, ...]
    evidence: tuple[str, ...]
    limitations: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "cell_index": self.cell_index,
            "cell_id": self.cell_id,
            "predicted_identity": self.predicted_identity,
            "identity_probability": self.identity_probability,
            "identity_confidence": self.identity_confidence,
            "developmental_stage": self.developmental_stage,
            "stage_probability": self.stage_probability,
            "damage_state": self.damage_state,
            "damage_probability": self.damage_probability,
            "viability_state": self.viability_state,
            "viability_probability": self.viability_probability,
            "cell_cycle_state": self.cell_cycle_state,
            "dominant_function": self.dominant_function,
            "predictions": [item.to_dict() for item in self.predictions],
            "evidence": list(self.evidence),
            "limitations": list(self.limitations),
        }


IDENTITY_PROGRAMS = {
    "erythroid": ("hematopoietic_erythroid", "Erythroid-like"),
    "myeloid": ("hematopoietic_myeloid", "Myeloid-like"),
    "megakaryocyte": ("hematopoietic_megakaryocyte", "Megakaryocyte-like"),
}


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(max(-60.0, min(60.0, -value))))


def _confidence(probability: float) -> str:
    return "high" if probability >= 0.8 else "moderate" if probability >= 0.6 else "low"


def _score_map(interpretation: CellInterpretation) -> dict[str, ProgramScore]:
    return {item.key: item for item in interpretation.scores if item.z_score is not None}


def _top(score_map: dict[str, ProgramScore], keys: Sequence[str]) -> tuple[str, float]:
    candidates = [(key, float(score_map[key].z_score or 0.0)) for key in keys if key in score_map]
    return max(candidates, key=lambda item: item[1]) if candidates else ("unavailable", 0.0)


def predict_cell(interpretation: CellInterpretation, energy: float | None = None, uncertainty: float | None = None) -> CellPrediction:
    scores = _score_map(interpretation)
    identity_key, identity_z = _top(scores, [item[0] for item in IDENTITY_PROGRAMS.values()])
    identity_name = IDENTITY_PROGRAMS.get(identity_key, (identity_key, "Unknown / unclassified"))[1]
    identity_probability = _sigmoid(identity_z - 0.55)
    plasticity = float(scores.get("transcriptional_plasticity", ProgramScore("", "", "", 0, 0, "", (), (), "", "")).z_score or 0.0)
    cycle = float(scores.get("cell_cycle", ProgramScore("", "", "", 0, 0, "", (), (), "", "")).z_score or 0.0)
    stress_keys = ["mitochondrial_stress", "endoplasmic_reticulum", "interferon_inflammation", "apoptosis"]
    stress_values = [float(scores[key].z_score or 0.0) for key in stress_keys if key in scores]
    stress = mean(stress_values)
    damage_probability = _sigmoid(stress - 0.75)
    viability_signal = float(scores.get("apoptosis", ProgramScore("", "", "", 0, 0, "", (), (), "", "")).z_score or 0.0)
    survival_signal = float(scores.get("mitochondrial_energy", ProgramScore("", "", "", 0, 0, "", (), (), "", "")).z_score or 0.0)
    viability_probability = _sigmoid(survival_signal - viability_signal)
    damage_state = "high stress signature" if damage_probability >= 0.7 else "moderate stress signature" if damage_probability >= 0.45 else "low stress signature"
    viability_state = "viability-associated programs" if viability_probability >= 0.6 else "mixed viability signals"
    cell_cycle_state = "cycling-like" if cycle >= 1.0 else "quiescent-like" if cycle <= -0.5 else "indeterminate cycle state"
    stage_signal = max(0.0, min(1.0, 0.5 + 0.2 * plasticity + 0.15 * identity_z))
    stage = "early / plastic" if stage_signal < 0.42 else "committed / transitional" if stage_signal < 0.68 else "late / fate-associated"
    function_key, function_z = _top(scores, ["mitochondrial_energy", "ribosome_translation", "endoplasmic_reticulum", "golgi_trafficking", "lysosome_autophagy", "cytoskeleton_motility", "antigen_presentation"])
    function = scores.get(function_key).name if function_key in scores else "No dominant function available"
    uncertainty_penalty = min(0.35, (uncertainty or 0.0) / 10.0)
    identity_probability = max(0.0, identity_probability - uncertainty_penalty)
    evidence = ("expression-derived marker programs", "relative score against available genes")
    limitations = ("Damage is a stress-program prediction, not a diagnosis of physical damage.", "Identity and stage require a reference atlas and batch-aware calibration for reliable biological use.", "A single RNA snapshot cannot establish causal mechanism or future fate.")
    predictions = (
        Prediction("identity", "Cell identity", identity_name, identity_probability, _confidence(identity_probability), (f"top identity program: {identity_key}",), limitations[:1]),
        Prediction("stage", "Developmental stage", stage, stage_signal, _confidence(stage_signal), (f"plasticity signal={plasticity:.2f}",), limitations[1:]),
        Prediction("damage", "Damage / stress", damage_state, damage_probability, _confidence(damage_probability), (f"aggregate stress z-score={stress:.2f}",), limitations[:1]),
        Prediction("viability", "Viability", viability_state, viability_probability, _confidence(viability_probability), (f"energy z-score={survival_signal:.2f}",), limitations[:1]),
        Prediction("cell_cycle", "Cell cycle", cell_cycle_state, _sigmoid(cycle), _confidence(_sigmoid(cycle)), (f"cell-cycle z-score={cycle:.2f}",), limitations[:1]),
        Prediction("function", "Dominant function", function, _sigmoid(function_z), _confidence(_sigmoid(function_z)), (f"program z-score={function_z:.2f}",), limitations[:1]),
    )
    return CellPrediction(interpretation.cell_index, interpretation.cell_id, identity_name, identity_probability, _confidence(identity_probability), stage, stage_signal, damage_state, damage_probability, viability_state, viability_probability, cell_cycle_state, function, predictions, evidence, limitations)


def predict_population(expression: Sequence[Sequence[float]], gene_names: Sequence[str], cell_ids: Sequence[str], labels: Sequence[str], energies: Sequence[float] | None = None, uncertainties: Sequence[float] | None = None) -> list[CellPrediction]:
    interpretations = [interpret_cell(index, values, gene_names, cell_ids[index], labels[index], energies[index] if energies else None, uncertainties[index] if uncertainties else None) for index, values in enumerate(expression)]
    return [predict_cell(item, energies[item.cell_index] if energies else None, uncertainties[item.cell_index] if uncertainties else None) for item in interpretations]


def prediction_summary(predictions: Sequence[CellPrediction]) -> dict[str, Any]:
    identities: dict[str, int] = {}
    damages: dict[str, int] = {}
    stages: dict[str, int] = {}
    for prediction in predictions:
        identities[prediction.predicted_identity] = identities.get(prediction.predicted_identity, 0) + 1
        damages[prediction.damage_state] = damages.get(prediction.damage_state, 0) + 1
        stages[prediction.developmental_stage] = stages.get(prediction.developmental_stage, 0) + 1
    return {"cells": len(predictions), "identity_counts": identities, "damage_counts": damages, "stage_counts": stages, "limitations": ["Predictions are calibrated heuristics over expression programs unless trained against a reference atlas."]}
