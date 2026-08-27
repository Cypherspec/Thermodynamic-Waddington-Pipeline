from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any, Sequence

from .arrays import mean, variance


@dataclass(frozen=True)
class Program:
    key: str
    name: str
    compartment: str
    markers: tuple[str, ...]
    output: str
    description: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ProgramScore:
    key: str
    name: str
    compartment: str
    score: float | None
    z_score: float | None
    evidence: str
    matched_markers: tuple[str, ...]
    missing_markers: tuple[str, ...]
    output: str
    description: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "key": self.key,
            "name": self.name,
            "compartment": self.compartment,
            "score": self.score,
            "z_score": self.z_score,
            "evidence": self.evidence,
            "matched_markers": list(self.matched_markers),
            "missing_markers": list(self.missing_markers),
            "output": self.output,
            "description": self.description,
        }


@dataclass(frozen=True)
class CellInterpretation:
    cell_index: int
    cell_id: str
    label: str
    evidence_level: str
    measured_modalities: tuple[str, ...]
    inferred_modalities: tuple[str, ...]
    unavailable_modalities: tuple[str, ...]
    dominant_programs: tuple[str, ...]
    scores: tuple[ProgramScore, ...]
    key_insights: tuple[str, ...]
    cautions: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "cell_index": self.cell_index,
            "cell_id": self.cell_id,
            "label": self.label,
            "evidence_level": self.evidence_level,
            "measured_modalities": list(self.measured_modalities),
            "inferred_modalities": list(self.inferred_modalities),
            "unavailable_modalities": list(self.unavailable_modalities),
            "dominant_programs": list(self.dominant_programs),
            "scores": [score.to_dict() for score in self.scores],
            "key_insights": list(self.key_insights),
            "cautions": list(self.cautions),
        }


PROGRAMS: tuple[Program, ...] = (
    Program("mitochondrial_energy", "Mitochondrial energy", "mitochondrion", ("MT-CO1", "MT-CO2", "MT-CO3", "ATP5F1E", "NDUFA1", "COX4I1"), "oxidative phosphorylation and ATP production", "Transcript-level proxy for mitochondrial respiratory activity."),
    Program("mitochondrial_stress", "Mitochondrial stress", "mitochondrion", ("ATF4", "DDIT3", "HSPA5", "CLPP", "LONP1", "OMA1"), "mitochondrial proteotoxic stress", "Stress-response program associated with mitochondrial and ER quality control."),
    Program("nucleus_chromatin", "Nuclear chromatin", "nucleus", ("H2AFZ", "HIST1H4C", "HMGN1", "LMNB1", "SMARCA4", "HDAC1"), "chromatin organization and nuclear maintenance", "Expression proxy for nuclear structure and chromatin maintenance."),
    Program("transcriptional_plasticity", "Transcriptional plasticity", "nucleus", ("JUN", "FOS", "MYC", "EGR1", "KLF4", "NANOG"), "rapid-response transcription and plasticity", "Activation program often associated with state transitions and plasticity."),
    Program("ribosome_translation", "Ribosome / translation", "ribosome", ("RPLP0", "RPL3", "RPL10", "RPS3", "RPS6", "RPS18"), "ribosome biogenesis and translation", "Transcript-level proxy for protein synthesis capacity."),
    Program("endoplasmic_reticulum", "Endoplasmic reticulum", "ER", ("CANX", "CALR", "SEC61A1", "HSPA5", "PDIA3", "DERL1"), "protein folding and ER quality control", "Expression proxy for ER folding, translocation, and stress handling."),
    Program("golgi_trafficking", "Golgi / trafficking", "Golgi", ("GOLGA2", "GOLGB1", "COPB2", "VCP", "RAB5A", "RAB7A"), "vesicle trafficking and cargo processing", "Expression proxy for secretory-pathway trafficking."),
    Program("lysosome_autophagy", "Lysosome / autophagy", "lysosome", ("LAMP1", "LAMP2", "CTSD", "CTSB", "SQSTM1", "BECN1"), "degradation and autophagic recycling", "Expression proxy for lysosomal degradation and autophagy."),
    Program("peroxisome_lipid", "Peroxisome / lipid", "peroxisome", ("PEX14", "PEX19", "ABCD3", "ACOX1", "CAT", "HSD17B4"), "peroxisomal oxidation and lipid metabolism", "Expression proxy for peroxisome maintenance and lipid oxidation."),
    Program("cytoskeleton_motility", "Cytoskeleton / motility", "cytoskeleton", ("ACTB", "TUBB", "VIM", "PFN1", "MYH9", "RHOA"), "shape, adhesion, and migration", "Expression proxy for cytoskeletal architecture and motility."),
    Program("extracellular_matrix", "Extracellular matrix", "matrix", ("COL1A1", "COL1A2", "FN1", "SPARC", "VIM", "ITGB1"), "matrix deposition, adhesion, and remodeling", "Expression proxy for extracellular structural remodeling."),
    Program("secretory_output", "Secretory output", "secretory pathway", ("SEC11C", "SEC24C", "SRP72", "MALAT1", "RAB1A", "STX5"), "protein secretion and export", "Expression proxy for secretory throughput, not direct protein release."),
    Program("cell_cycle", "Cell cycle", "nucleus", ("MKI67", "TOP2A", "PCNA", "TUBA1B", "TYMS", "CDC20"), "proliferation and cell-cycle progression", "Expression proxy for proliferative state."),
    Program("apoptosis", "Apoptosis", "whole-cell", ("BAX", "BAK1", "CASP3", "CASP8", "BCL2", "MCL1"), "survival / programmed cell death balance", "Transcript-level proxy; cannot establish that apoptosis is occurring."),
    Program("hypoxia", "Hypoxia response", "whole-cell", ("HIF1A", "VEGFA", "LDHA", "SLC2A1", "EGLN3", "CA9"), "low-oxygen response and glycolytic adaptation", "Expression proxy for hypoxia-associated transcriptional response."),
    Program("dna_damage", "DNA damage response", "nucleus", ("GADD45A", "CDKN1A", "DDB2", "ATM", "ATR", "TP53"), "genome surveillance and repair", "Expression proxy for damage-response activation."),
    Program("antigen_presentation", "Antigen presentation", "secretory / membrane", ("H2-AA", "H2-AB1", "CD74", "CIITA", "B2M", "TAP1"), "MHC-II and antigen-presentation capacity", "Expression proxy for antigen processing and presentation."),
    Program("interferon_inflammation", "Interferon / inflammation", "whole-cell", ("IFIT1", "IFIT3", "ISG15", "STAT1", "IRF7", "TNF"), "innate immune activation", "Expression proxy for interferon and inflammatory signaling."),
    Program("hematopoietic_erythroid", "Erythroid commitment", "lineage", ("GATA1", "KLF1", "HBB-BT", "HBA-A1", "ALAS2", "EPOR"), "erythroid differentiation", "Lineage-associated program for hematopoietic erythroid fate."),
    Program("hematopoietic_myeloid", "Myeloid commitment", "lineage", ("SPI1", "CEBPA", "LYZ", "CSF1R", "FCER1G", "CTSS"), "myeloid differentiation", "Lineage-associated program for hematopoietic myeloid fate."),
    Program("hematopoietic_megakaryocyte", "Megakaryocyte commitment", "lineage", ("PPBP", "PF4", "NFE2", "ITGA2B", "GP9", "RGS18"), "megakaryocyte differentiation", "Lineage-associated program for megakaryocyte fate."),
)


def default_gene_names(n_features: int) -> list[str]:
    return [f"feature_{index}" for index in range(n_features)]


def _canonical_gene(name: str) -> str:
    value = str(name).strip().upper().replace("-", "").replace("_", "")
    if value.startswith("MT"):
        value = value[2:]
    return value


def _gene_index(gene_names: Sequence[str]) -> dict[str, int]:
    index: dict[str, int] = {}
    for position, raw_name in enumerate(gene_names):
        for alias in str(raw_name).replace(";", "|").split("|"):
            normalized = _canonical_gene(alias)
            if normalized:
                index.setdefault(normalized, position)
    return index


def _score(values: Sequence[float], indices: Sequence[int]) -> tuple[float, float]:
    selected = [float(values[index]) for index in indices]
    center = mean(selected)
    background = [float(value) for index, value in enumerate(values) if index not in set(indices)]
    background_mean = mean(background)
    background_sd = math.sqrt(max(1e-12, variance(background)))
    return center, (center - background_mean) / background_sd


def score_program(program: Program, values: Sequence[float], gene_names: Sequence[str]) -> ProgramScore:
    index = _gene_index(gene_names)
    matched = tuple(marker for marker in program.markers if _canonical_gene(marker) in index)
    missing = tuple(marker for marker in program.markers if marker.upper() not in index)
    if not matched:
        return ProgramScore(program.key, program.name, program.compartment, None, None, "unavailable", matched, missing, program.output, program.description)
    indices = [index[_canonical_gene(marker)] for marker in matched]
    score, z_score = _score(values, indices)
    coverage = len(matched) / max(1, len(program.markers))
    evidence = "measured-expression-program" if coverage >= 0.5 else "sparse-expression-program"
    return ProgramScore(program.key, program.name, program.compartment, score, z_score, evidence, matched, missing, program.output, program.description)


def _insights(scores: Sequence[ProgramScore], label: str, energy: float | None, uncertainty: float | None) -> tuple[str, ...]:
    available = [score for score in scores if score.z_score is not None]
    high = sorted(available, key=lambda score: float(score.z_score), reverse=True)
    insights = []
    for score in high[:3]:
        if score.z_score is not None and score.z_score >= 1.0:
            insights.append(f"Elevated {score.name.lower()} program; inferred output: {score.output}.")
    lineage = [score for score in available if score.compartment == "lineage"]
    if lineage:
        leader = max(lineage, key=lambda score: float(score.z_score or -math.inf))
        if leader.z_score is not None and leader.z_score >= 0.75:
            insights.append(f"The strongest lineage-associated program is {leader.name} ({leader.z_score:+.2f} z-score), consistent with {leader.output}.")
    if energy is not None:
        if energy <= 0:
            insights.append(f"The cell sits near the low end of the fitted effective-energy coordinate (F={energy:.3f}).")
        else:
            insights.append(f"The cell sits at F={energy:.3f} on the fitted effective-energy coordinate; compare with local barriers, not energy alone.")
    if uncertainty is not None and uncertainty > 1:
        insights.append(f"Local landscape uncertainty is high (σ={uncertainty:.3f}); interpret fate and organelle inferences conservatively.")
    return tuple(insights or [f"No program exceeded the interpretation threshold for {label or 'this cell'}."])


def interpret_cell(index: int, values: Sequence[float], gene_names: Sequence[str], cell_id: str, label: str, energy: float | None = None, uncertainty: float | None = None, lineage_outcome: str | None = None) -> CellInterpretation:
    scores = tuple(score_program(program, values, gene_names) for program in PROGRAMS)
    available = [score for score in scores if score.score is not None]
    dominant = tuple(score.name for score in sorted(available, key=lambda item: float(item.z_score or -math.inf), reverse=True)[:5])
    measured = ("RNA expression",)
    inferred = ("organelle programs", "functional outputs", "cell-state programs")
    unavailable = ("direct organelle morphology", "protein abundance", "metabolites", "subcellular localization")
    if lineage_outcome is not None:
        measured = measured + ("lineage outcome",)
    evidence_level = "program-level inference" if available else "unavailable"
    cautions = ["Organelle states are inferred from RNA marker programs, not directly observed.", "A program score is not proof that the corresponding organelle function is active.", "Correlated programs and ambient RNA can produce false positives."]
    if lineage_outcome is not None:
        cautions.append("The lineage outcome is an observed endpoint, not a causal explanation of the transcriptome.")
    return CellInterpretation(index, cell_id, label, evidence_level, measured, inferred, unavailable, dominant, scores, _insights(scores, label, energy, uncertainty), tuple(cautions))


def interpret_population(expression: Sequence[Sequence[float]], gene_names: Sequence[str] | None = None, cell_ids: Sequence[str] | None = None, labels: Sequence[str] | None = None, energies: Sequence[float] | None = None, uncertainties: Sequence[float] | None = None, lineage_outcomes: Sequence[str] | None = None) -> list[CellInterpretation]:
    if not expression:
        return []
    names = list(gene_names or default_gene_names(len(expression[0])))
    ids = list(cell_ids or [f"cell_{index:05d}" for index in range(len(expression))])
    cell_labels = list(labels or ["unlabeled"] * len(expression))
    return [interpret_cell(index, values, names, ids[index], cell_labels[index], energies[index] if energies and index < len(energies) else None, uncertainties[index] if uncertainties and index < len(uncertainties) else None, lineage_outcomes[index] if lineage_outcomes and index < len(lineage_outcomes) else None) for index, values in enumerate(expression)]


def atlas_summary(interpreted: Sequence[CellInterpretation]) -> dict[str, Any]:
    if not interpreted:
        return {"cells": 0, "programs": [], "evidence_levels": {}, "modalities": {}}
    program_counts: dict[str, int] = {}
    evidence_counts: dict[str, int] = {}
    for cell in interpreted:
        evidence_counts[cell.evidence_level] = evidence_counts.get(cell.evidence_level, 0) + 1
        for program in cell.dominant_programs:
            program_counts[program] = program_counts.get(program, 0) + 1
    return {"cells": len(interpreted), "programs": sorted(program_counts.items(), key=lambda item: item[1], reverse=True), "evidence_levels": evidence_counts, "modalities": {"measured": sorted({item for cell in interpreted for item in cell.measured_modalities}), "inferred": sorted({item for cell in interpreted for item in cell.inferred_modalities}), "unavailable": sorted({item for cell in interpreted for item in cell.unavailable_modalities})}, "scientific_boundary": "Cell and organelle interpretations are marker-program inferences from the supplied expression modality."}
