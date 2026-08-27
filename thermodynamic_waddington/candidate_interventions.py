from __future__ import annotations

"""Evidence-backed intervention hypotheses for hematopoietic fate studies.

These entries are literature-supported hypotheses, not measured effects in the
current Thermodynamic Waddington project. Numeric intervention vectors must be
learned from a matched perturbation dataset; this module never fabricates them.
"""

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence


@dataclass(frozen=True)
class InterventionHypothesis:
    name: str
    target_fates: tuple[str, ...]
    direction: str
    rationale: str
    evidence_level: str
    expected_readouts: tuple[str, ...]
    safety_flags: tuple[str, ...]
    sources: tuple[str, ...]
    measured_effect_status: str = "not_measured_in_this_project"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def hematopoietic_registry() -> list[InterventionHypothesis]:
    return [
        InterventionHypothesis(
            "SPI1_activation",
            ("Monocyte", "Neutrophil", "Myeloid"),
            "activation_hypothesis",
            "PU.1/SPI1 is a canonical myeloid regulator and has been experimentally linked to myeloid commitment; effect size and dose response must be measured in the selected system.",
            "literature_supported_hypothesis",
            ("target-fate fraction", "myeloid marker program", "viability", "clone-level fate entropy"),
            ("overactivation may alter proliferation or viability", "do not infer clinical suitability"),
            ("https://pubmed.ncbi.nlm.nih.gov/9694804", "https://pubmed.ncbi.nlm.nih.gov/10453070"),
        ),
        InterventionHypothesis(
            "IRF8_activation",
            ("Monocyte", "Dendritic"),
            "activation_hypothesis",
            "IRF8 has been reported to govern enhancer dynamics in mononuclear phagocyte progenitors and is required for monocyte/DC development.",
            "literature_supported_hypothesis",
            ("monocyte/DC fate", "enhancer-linked program", "viability", "off-target DC mass"),
            ("context and dosage dependence", "immune-state effects require orthogonal functional assays"),
            ("https://pubmed.ncbi.nlm.nih.gov/29514092"),
        ),
        InterventionHypothesis(
            "CEBPA_activation",
            ("Myeloid", "Neutrophil", "Monocyte"),
            "activation_hypothesis",
            "C/EBPα is a canonical hematopoietic lineage regulator and a mechanistic candidate for myeloid differentiation programs.",
            "literature_supported_hypothesis",
            ("myeloid fate", "granulocyte/monocyte balance", "cell-cycle state", "viability"),
            ("dose and temporal schedule can change outcome", "malignant-context effects are not transferable to normal progenitors"),
            ("https://pubmed.ncbi.nlm.nih.gov/21720200", "https://pubmed.ncbi.nlm.nih.gov/17890457"),
        ),
        InterventionHypothesis(
            "GATA1_activation",
            ("Erythroid", "Megakaryocyte"),
            "lineage_contrast_control",
            "GATA1 provides a biologically motivated contrast arm for the PU.1/GATA1 antagonistic hematopoietic decision axis.",
            "literature_supported_contrast",
            ("erythroid/megakaryocyte fate", "globin/platelet program", "off-target myeloid suppression", "viability"),
            ("not a monocyte-promoting intervention", "use as a contrast or positive-control arm only"),
            ("https://pubmed.ncbi.nlm.nih.gov/10364157", "https://pubmed.ncbi.nlm.nih.gov/24638828"),
        ),
        InterventionHypothesis(
            "non_targeting_control",
            ("All" ,),
            "negative_control",
            "Matched negative-control intervention for estimating handling, delivery, and perturbation-assignment effects.",
            "experimental_control_required",
            ("viability", "barcode recovery", "perturbation capture rate", "baseline fate distribution"),
            ("must be randomized and processed identically"),
            (),
        ),
    ]


def registry_payload() -> dict[str, Any]:
    return {
        "registry": "TW-HEM-FATE-CANDIDATES-v1",
        "status": "hypotheses_and_controls_only",
        "claim_boundary": "Literature support does not establish an effect in the proposed experiment; all effect sizes and causal conclusions require measured data.",
        "candidates": [item.to_dict() for item in hematopoietic_registry()],
    }


def write_registry(path: str | Path = "experiments/intervention_registry.json") -> dict[str, Any]:
    payload = registry_payload()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload
