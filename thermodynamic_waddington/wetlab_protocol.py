from __future__ import annotations

"""Experiment-ready, but not experiment-claiming, wet-lab validation plans."""

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Sequence


@dataclass(frozen=True)
class WetLabProtocol:
    protocol_id: str
    biological_system: str
    target_fate: str
    intervention_names: tuple[str, ...]
    biological_replicates: int
    donors_or_animals: int
    timepoints: tuple[str, ...]
    primary_endpoint: str
    secondary_endpoints: tuple[str, ...]
    negative_controls: tuple[str, ...]
    positive_controls: tuple[str, ...]
    randomization: str
    blinding: str
    acceptance_criteria: tuple[str, ...]
    exclusion_criteria: tuple[str, ...]
    required_assays: tuple[str, ...]
    preregistration_requirements: tuple[str, ...] = ()
    data_release_requirements: tuple[str, ...] = ()
    claim_status: str = "not_performed"
    safety_notes: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def build_hematopoietic_fate_protocol(target_fate: str, interventions: Sequence[str], biological_replicates: int = 3, donors_or_animals: int = 3) -> WetLabProtocol:
    if biological_replicates < 3:
        raise ValueError("use at least three biological replicates for the preregistered plan")
    if donors_or_animals < 2:
        raise ValueError("use at least two independent donors or animals")
    return WetLabProtocol(
        protocol_id="TW-HEM-FATE-001",
        biological_system="mouse hematopoietic progenitor differentiation with matched lineage/barcode readout",
        target_fate=target_fate,
        intervention_names=tuple(interventions),
        biological_replicates=biological_replicates,
        donors_or_animals=donors_or_animals,
        timepoints=("baseline", "early", "commitment", "endpoint"),
        primary_endpoint="pre-registered target-fate fraction among viable, barcode-resolved cells at endpoint",
        secondary_endpoints=("clone-level fate entropy", "target-fate probability calibration", "barrier-crossing timing", "cell-state viability", "off-target fate mass"),
        negative_controls=("vehicle or non-targeting guide", "sham handling", "label permutation analysis", "expression-only computational baseline"),
        positive_controls=("known fate-promoting condition", "known fate-inhibiting condition", "technical library control"),
        randomization="randomize wells, processing order, and sequencing lanes within donor and intervention blocks",
        blinding="pre-register cell and well IDs; blind outcome annotation and analysis to intervention identity until the locked analysis completes",
        acceptance_criteria=("all primary endpoints and exclusions are locked before unblinding", "target-fate gain exceeds the pre-registered minimum effect", "effect is directionally consistent across independent donors", "calibration interval contains the observed fate fraction", "negative controls do not reproduce the intervention effect"),
        exclusion_criteria=("failed library QC", "unresolved barcode-to-cell mapping", "pre-registered viability threshold failure", "cross-contamination or sample swap", "post hoc removal of unfavorable replicates"),
        required_assays=("single-cell RNA sequencing", "spliced/unspliced RNA when velocity is claimed", "lineage barcode capture", "perturbation identity capture", "viability and library QC", "orthogonal endpoint fate assay"),
        preregistration_requirements=("register protocol and analysis plan before data unblinding", "lock primary endpoint, exclusion rules, donor-level analysis, and stopping rules", "publish the preregistration hash and versioned code commit"),
        data_release_requirements=("release de-identified cell-level metadata and lineage/barcode mapping", "release raw or accession-linked counts, QC, perturbation identity, viability, and batch fields", "release analysis manifest, randomization map, and failed-sample accounting"),
        safety_notes=("This is a research design artifact, not a clinical protocol.", "Use institutionally approved biosafety, animal-care, and genetic-perturbation procedures.", "Do not infer therapeutic safety from computational rankings."),
    )


def write_protocol(protocol: WetLabProtocol, path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(protocol.to_dict(), indent=2, sort_keys=True), encoding="utf-8")

