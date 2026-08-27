from __future__ import annotations

"""Deterministic synthetic wet-lab study generator.

This creates a complete mock study with randomized donors, interventions, wells,
lineage barcodes, viability, orthogonal fate calls, and target probabilities. It
is useful for exercising the full pipeline, but it is never biological evidence.
"""

import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .causal_validation import StudyCell, StudyConfig, StudyReport, analyze_study


@dataclass(frozen=True)
class SimulationConfig:
    seed: int = 20260728
    donors: int = 8
    cells_per_group: int = 120
    target_fate: str = "Monocyte"
    control: str = "vehicle"
    interventions: tuple[str, ...] = ("SPI1_activation", "IRF8_activation", "CEBPA_activation", "GATA1_activation")
    target_rates: tuple[float, ...] = (0.48, 0.42, 0.45, 0.06)
    control_rate: float = 0.18
    viability: float = 0.94
    effect_noise: float = 0.018
    timepoints: tuple[str, ...] = ("baseline", "early", "commitment", "endpoint")
    baseline_fraction: float = 0.22
    modality_noise: float = 0.04
    barcode_capture: float = 0.88
    perturbation_capture: float = 0.93
    orthogonal_agreement: float = 0.96

    def validate(self) -> None:
        if self.donors < 3 or self.cells_per_group < 10:
            raise ValueError("simulation requires at least 3 donors and 10 cells per group")
        if len(self.interventions) != len(self.target_rates):
            raise ValueError("interventions and target_rates must have equal length")
        if not 0.0 <= self.control_rate <= 1.0:
            raise ValueError("control_rate must be in [0, 1]")
        if any(not 0.0 <= rate <= 1.0 for rate in self.target_rates):
            raise ValueError("target_rates must be in [0, 1]")
        if not self.timepoints or self.timepoints[-1] != "endpoint":
            raise ValueError("timepoints must end with endpoint")
        for name, value in (("viability", self.viability), ("baseline_fraction", self.baseline_fraction), ("barcode_capture", self.barcode_capture), ("perturbation_capture", self.perturbation_capture), ("orthogonal_agreement", self.orthogonal_agreement)):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be in [0, 1]")


def _fate(rng: random.Random, target: str, probability: float) -> str:
    if rng.random() < probability:
        return target
    other = ["Neutrophil", "Erythroid", "Megakaryocyte", "Lymphoid"]
    return other[rng.randrange(len(other))]


def generate_simulated_cells(config: SimulationConfig | None = None) -> list[StudyCell]:
    config = config or SimulationConfig()
    config.validate()
    rng = random.Random(config.seed)
    groups = [(config.control, config.control_rate)] + list(zip(config.interventions, config.target_rates))
    cells: list[StudyCell] = []
    for donor_index in range(config.donors):
        donor = f"sim_donor_{donor_index + 1:02d}"
        donor_shift = rng.gauss(0.0, config.effect_noise)
        for intervention, nominal_rate in groups:
            well = f"{donor}_{intervention}_well_{rng.randrange(1, 5):02d}"
            rate = max(0.0, min(1.0, nominal_rate + donor_shift))
            for cell_index in range(config.cells_per_group):
                cell_id = f"{donor}:{intervention}:{cell_index:04d}"
                lineage_id = f"{donor}:{intervention}:clone_{cell_index // 6:03d}" if rng.random() < config.barcode_capture else None
                viable = rng.random() < config.viability
                final_probability = max(0.0, min(1.0, rate + rng.gauss(0.0, config.effect_noise)))
                fate = _fate(rng, config.target_fate, final_probability) if viable else None
                perturbation_observed = intervention if rng.random() < config.perturbation_capture else "unassigned_perturbation"
                for timepoint_index, timepoint in enumerate(config.timepoints):
                    progress = timepoint_index / max(1, len(config.timepoints) - 1)
                    probability = config.baseline_fraction + progress * (final_probability - config.baseline_fraction)
                    probability = max(0.0, min(1.0, probability + rng.gauss(0.0, config.modality_noise)))
                    orthogonal_fate = _fate(rng, config.target_fate, probability) if viable and timepoint == "endpoint" and rng.random() < config.orthogonal_agreement else (fate if timepoint == "endpoint" else None)
                    cells.append(StudyCell(
                        cell_id=f"{cell_id}:{timepoint}",
                        donor=donor,
                        replicate=f"replicate_{donor_index + 1:02d}",
                        intervention=perturbation_observed,
                        fate=orthogonal_fate,
                        target_probability=probability,
                        viable=viable,
                        lineage_id=lineage_id,
                        timepoint=timepoint,
                        batch=f"batch_{(donor_index % 3) + 1}",
                        dose=0.0 if intervention == config.control else 1.0,
                        orthogonal_fate=orthogonal_fate if timepoint == "endpoint" else None,
                        barcode_captured=lineage_id is not None,
                        perturbation_captured=perturbation_observed == intervention,
                        expression_observed=True,
                        velocity_observed=True,
                    ))
    return cells


def cells_payload(cells: list[StudyCell], config: SimulationConfig) -> dict[str, Any]:
    return {
        "simulation": True,
        "simulation_status": "synthetic_not_biological_evidence",
        "generator": "thermodynamic_waddington.wetlab_simulation",
        "config": asdict(config),
        "cells": [cell.to_dict() for cell in cells],
        "claim_boundary": "Synthetic outcomes are generated from declared probabilities and cannot validate a biological intervention.",
    }


def write_simulated_study(path: str | Path = "data/simulated/wetlab/study.json", config: SimulationConfig | None = None) -> dict[str, Any]:
    config = config or SimulationConfig()
    cells = generate_simulated_cells(config)
    payload = cells_payload(cells, config)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def analyze_simulated_study(path: str | Path = "data/simulated/wetlab/study.json", output: str | Path = "experiments/simulated_wetlab_validation.json", config: SimulationConfig | None = None) -> StudyReport:
    config = config or SimulationConfig()
    payload = write_simulated_study(path, config)
    cells = [StudyCell.from_dict(row) for row in payload["cells"]]
    report = analyze_study(cells, StudyConfig(config.target_fate, config.control, minimum_donors=config.donors, minimum_cells_per_group=config.cells_per_group, bootstrap_rounds=1000, require_lineage=True, require_orthogonal_fate=True), "synthetic-wetlab-study")
    report.provenance.update({"simulation": True, "simulation_seed": config.seed, "source_path": str(path), "claim_boundary": "The study passes the software gate only because it is synthetic. It is not wet-lab validation."})
    report.save(output)
    return report
