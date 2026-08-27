from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any, Sequence

from .cell_atlas import PROGRAMS, CellInterpretation, interpret_cell
from .arrays import mean


@dataclass(frozen=True)
class Evidence:
    domain: str
    statement: str
    strength: float
    basis: tuple[str, ...]
    boundary: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CellPhenotype:
    cell_index: int
    cell_id: str
    observed_label: str
    identity: str
    identity_confidence: float
    developmental_stage: str
    stage_confidence: float
    damage_state: str
    damage_score: float
    viability_state: str
    viability_score: float
    cycle_state: str
    dominant_outputs: tuple[str, ...]
    compartment_activity: dict[str, float]
    evidence: tuple[Evidence, ...]
    caveats: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["evidence"] = [item.to_dict() for item in self.evidence]
        payload["dominant_outputs"] = list(self.dominant_outputs)
        payload["caveats"] = list(self.caveats)
        return payload


def _z(scores: dict[str, Any], key: str) -> float:
    item = scores.get(key)
    return float(item.z_score or 0.0) if item is not None and item.z_score is not None else 0.0


def _sigmoid(value: float) -> float:
    return 1.0 / (1.0 + math.exp(max(-50.0, min(50.0, -value))))


def _confidence(score: float, uncertainty: float | None) -> float:
    penalty = min(0.35, max(0.0, float(uncertainty or 0.0)) / 12.0)
    return max(0.0, min(1.0, _sigmoid(score) - penalty))


def interpret_phenotype(interpretation: CellInterpretation, energy: float | None = None, uncertainty: float | None = None) -> CellPhenotype:
    scores = {item.key: item for item in interpretation.scores if item.z_score is not None}
    lineage = [(key, _z(scores, key)) for key in ("hematopoietic_erythroid", "hematopoietic_myeloid", "hematopoietic_megakaryocyte") if key in scores]
    identity_key, identity_z = max(lineage, key=lambda item: item[1]) if lineage else ("unknown", 0.0)
    identity_names = {"hematopoietic_erythroid": "erythroid-like", "hematopoietic_myeloid": "myeloid-like", "hematopoietic_megakaryocyte": "megakaryocyte-like", "unknown": "unclassified"}
    identity = identity_names[identity_key]
    plasticity = _z(scores, "transcriptional_plasticity")
    cycle = _z(scores, "cell_cycle")
    stress = mean([_z(scores, key) for key in ("mitochondrial_stress", "endoplasmic_reticulum", "interferon_inflammation", "apoptosis")])
    energy_signal = _z(scores, "mitochondrial_energy")
    damage_score = _sigmoid(stress - 0.65)
    viability_score = _sigmoid(energy_signal - _z(scores, "apoptosis"))
    damage_state = "high stress / damage-associated signature" if damage_score >= 0.7 else "moderate stress signature" if damage_score >= 0.45 else "low stress signature"
    viability_state = "viability-associated expression" if viability_score >= 0.65 else "mixed viability evidence" if viability_score >= 0.4 else "low viability-associated expression"
    cycle_state = "cycling-like" if cycle >= 1.0 else "quiescent-like" if cycle <= -0.5 else "indeterminate cycle state"
    stage_score = max(-3.0, min(3.0, plasticity * 0.65 + identity_z * 0.35))
    stage_confidence = _confidence(abs(stage_score), uncertainty)
    stage = "early / plastic" if stage_score < -0.25 else "transitional / committing" if stage_score < 0.75 else "late / fate-associated"
    available = sorted((item for item in scores.values() if item.z_score is not None), key=lambda item: float(item.z_score or 0.0), reverse=True)
    outputs = tuple(dict.fromkeys(item.output for item in available[:5]))
    if not outputs:
        outputs = ("global transcriptome activity (non-specific)", "directional RNA-velocity signal")
    compartments: dict[str, float] = {}
    for item in scores.values():
        compartments[item.compartment] = max(compartments.get(item.compartment, -math.inf), float(item.z_score or 0.0))
    evidence = (
        Evidence("identity", f"Highest available lineage program: {identity}.", _confidence(identity_z, uncertainty), tuple(item.name for item in scores.values() if item.key == identity_key), "RNA marker-program evidence; not a reference-atlas probability."),
        Evidence("stage", f"Plasticity and lineage scores place the cell in the {stage} region.", stage_confidence, (f"plasticity z={plasticity:.2f}", f"lineage z={identity_z:.2f}"), "A snapshot cannot prove developmental time or future fate."),
        Evidence("damage", f"The aggregate stress signature is {damage_state}.", damage_score, (f"stress z={stress:.2f}",), "Stress programs are not direct evidence of physical damage, apoptosis, or loss of viability."),
        Evidence("function", f"Dominant transcript-level outputs include: {', '.join(outputs[:3])}.", _confidence(float(available[0].z_score or 0.0) if available else 0.0, uncertainty), tuple(item.name for item in available[:3]) or ("global transcriptome amplitude", "RNA velocity"), "Specific organelle or functional assignment is not identifiable without canonical marker genes; global RNA activity is reported instead."),
    )
    caveats = (
        "This is an interpretable RNA-derived phenotype report, not a diagnosis of the cell.",
        "Organelle morphology, protein abundance, metabolite flux, subcellular localization, and ultrastructure are not directly observed.",
        "Identity and stage need a tissue-specific reference atlas, batch correction, and external validation for quantitative claims.",
        "The effective landscape coordinate is separate evidence and does not by itself establish biological commitment.",
    )
    return CellPhenotype(interpretation.cell_index, interpretation.cell_id, interpretation.label, identity, _confidence(identity_z, uncertainty), stage, stage_confidence, damage_state, damage_score, viability_state, viability_score, cycle_state, outputs, compartments, evidence, caveats)


def phenotype_population(expression: Sequence[Sequence[float]], gene_names: Sequence[str], cell_ids: Sequence[str], labels: Sequence[str], energies: Sequence[float] | None = None, uncertainties: Sequence[float] | None = None, lineage_outcomes: Sequence[str] | None = None) -> list[CellPhenotype]:
    interpretations = [interpret_cell(index, values, gene_names, cell_ids[index], labels[index], energies[index] if energies else None, uncertainties[index] if uncertainties else None, lineage_outcomes[index] if lineage_outcomes else None) for index, values in enumerate(expression)]
    return [interpret_phenotype(item, energies[item.cell_index] if energies else None, uncertainties[item.cell_index] if uncertainties else None) for item in interpretations]


def phenotype_summary(phenotypes: Sequence[CellPhenotype]) -> dict[str, Any]:
    def counts(values: Sequence[str]) -> dict[str, int]:
        result: dict[str, int] = {}
        for value in values:
            result[value] = result.get(value, 0) + 1
        return result
    return {"cells": len(phenotypes), "identities": counts([item.identity for item in phenotypes]), "stages": counts([item.developmental_stage for item in phenotypes]), "damage_states": counts([item.damage_state for item in phenotypes]), "viability_states": counts([item.viability_state for item in phenotypes]), "cycle_states": counts([item.cycle_state for item in phenotypes]), "caveats": ["All phenotype fields are expression-program inferences unless a modality is explicitly listed as measured."]}
