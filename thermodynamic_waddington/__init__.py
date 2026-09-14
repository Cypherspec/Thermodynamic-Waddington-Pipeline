"""Thermodynamic Waddington: effective free-energy inference for cell fate."""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

try:
    __version__ = _pkg_version("thermodynamic-waddington")
except PackageNotFoundError:
    __version__ = "0.3.0"

from .config import FitConfig
from .model import LandscapeFit, fit_landscape
from .synthetic import SyntheticDataset, make_synthetic_dataset
from .causal import InterventionEffect, intervention_scan, rank_targets
from .geometry import MetricSummary, summarize
from .streaming import OnlineLandscape

__all__ = [
    "FitConfig",
    "LandscapeFit",
    "SyntheticDataset",
    "fit_landscape",
    "make_synthetic_dataset",
    "InterventionEffect",
    "intervention_scan",
    "rank_targets",
    "MetricSummary",
    "summarize",
    "OnlineLandscape",
    "DiscoveryConfig",
    "DiscoveryLedger",
    "build_discovery_ledger",
    "write_discovery_ledger",
    "PredictionMetrics", "PairedComparison", "ChampionshipReport", "evaluate_methods", "paired_compare", "fingerprint_report",
    "NonequilibriumConfig", "NonequilibriumField", "estimate_field", "path_action", "compare_reversible_irreversible",
    "WetLabProtocol", "build_hematopoietic_fate_protocol", "write_protocol",
    "StudyCell", "StudyConfig", "StudyReport", "analyze_study", "create_preregistration",
]

from .discovery_engine import DiscoveryConfig, DiscoveryLedger, build_discovery_ledger, write_discovery_ledger

from .benchmark_championship import ChampionshipReport, PairedComparison, PredictionMetrics, evaluate_methods, fingerprint_report, paired_compare
from .latent_nonequilibrium import NonequilibriumConfig, NonequilibriumField, compare_reversible_irreversible, estimate_field, path_action

from .wetlab_protocol import WetLabProtocol, build_hematopoietic_fate_protocol, write_protocol
from .causal_validation import StudyCell, StudyConfig, StudyReport, analyze_study, create_preregistration, load_cells
from .closed_loop import ClosedLoopConfig, ClosedLoopEngine, ClosedLoopState, MeasuredOutcome, write_round_plan
__all__ += ["WetLabProtocol", "build_hematopoietic_fate_protocol", "write_protocol", "StudyCell", "StudyConfig", "StudyReport", "analyze_study", "create_preregistration", "load_cells", "ClosedLoopConfig", "ClosedLoopEngine", "ClosedLoopState", "MeasuredOutcome", "write_round_plan"]

from .candidate_interventions import InterventionHypothesis, hematopoietic_registry, registry_payload, write_registry
from .measured_study_manifest import empty_manifest, validate_manifest, write_template
from .larry_fate_model import FateModelConfig, run_real_larry_model
__all__ += ["FateModelConfig", "run_real_larry_model"]
from .flagship_benchmark import Observation, BenchmarkReport, evaluate as evaluate_flagship_benchmark, write_report as write_flagship_benchmark
__all__ += ["Observation", "BenchmarkReport", "evaluate_flagship_benchmark", "write_flagship_benchmark"]
from .real_intervention_bundle import IngestResult, ingest_measured_study
__all__ += ["IngestResult", "ingest_measured_study"]
from .claim_gate import ClaimDecision, EvidenceItem, decide as decide_claim
__all__ += ["ClaimDecision", "EvidenceItem", "decide_claim"]
from .award_capsule import CapsuleClaim, ResearchCapsule, build_capsule, write_capsule
__all__ += ["CapsuleClaim", "ResearchCapsule", "build_capsule", "write_capsule"]
from .review_packet import build_review_packet, write_review_packet
__all__ += ["build_review_packet", "write_review_packet"]

from .falsification import FalsificationCase, run_falsification_suite, write_falsification_report
__all__ += ["FalsificationCase", "run_falsification_suite", "write_falsification_report"]

from .evidence_audit import AuditCheck, build_evidence_audit, write_evidence_audit
__all__ += ["AuditCheck", "build_evidence_audit", "write_evidence_audit"]

from .replication import ReplicationStudy, ReplicationReport, evaluate_replication, load_studies, write_replication_report
__all__ += ["ReplicationStudy", "ReplicationReport", "evaluate_replication", "load_studies", "write_replication_report"]
from .intervention_trial import TrialObservation, TrialConfig, TrialReport, analyze_trial, load_observations, write_preregistration, write_trial_report
__all__ += ["TrialObservation", "TrialConfig", "TrialReport", "analyze_trial", "load_observations", "write_preregistration", "write_trial_report"]
from .transportability import StudyObservation, TransportabilityConfig, TransportabilityReport, evaluate_transportability, load_transportability_manifest, write_transportability_report
__all__ += ["StudyObservation", "TransportabilityConfig", "TransportabilityReport", "evaluate_transportability", "load_transportability_manifest", "write_transportability_report"]
from .robustness import evaluate_robustness, run_synthetic_robustness, write_robustness_report
__all__ += ["evaluate_robustness", "run_synthetic_robustness", "write_robustness_report"]
from .path_ensemble import PathEnsembleConfig, PathEnsembleReport, PathWorkSample, evaluate_path_ensemble, write_path_ensemble_report
__all__ += ["PathEnsembleConfig", "PathEnsembleReport", "PathWorkSample", "evaluate_path_ensemble", "write_path_ensemble_report"]
from .decision_engine import DecisionConfig, DecisionReport, EvidenceRecord, decide, decide_from_files, write_decision_report
__all__ += ["DecisionConfig", "DecisionReport", "EvidenceRecord", "decide", "decide_from_files", "write_decision_report"]

from .evidence_release import ReleaseCheck, build_release_report, write_release_report
__all__ += ["ReleaseCheck", "build_release_report", "write_release_report"]

from .evidence_graph import EvidenceGraph, build_evidence_graph, write_evidence_graph
__all__ += ["EvidenceGraph", "build_evidence_graph", "write_evidence_graph"]

from .translation_gate import EvidenceRecord, TranslationDecision, build_translation_report, evaluate_translation, write_translation_report
__all__ += ["EvidenceRecord", "TranslationDecision", "build_translation_report", "evaluate_translation", "write_translation_report"]
from .premium_platform import PlatformRun, run_premium_platform, write_premium_report
__all__ += ["PlatformRun", "run_premium_platform", "write_premium_report"]

from .optimal_transport import (
    SinkhornConfig, TransportPlan, BridgeChain,
    sinkhorn, squared_euclidean_cost, estimate_growth_weights,
    displacement_interpolate, fit_schrodinger_chain, chain_diagnostic_report,
)
__all__ += [
    "SinkhornConfig", "TransportPlan", "BridgeChain",
    "sinkhorn", "squared_euclidean_cost", "estimate_growth_weights",
    "displacement_interpolate", "fit_schrodinger_chain", "chain_diagnostic_report",
]

from .calibration import CalibrationReport, basin_barriers_kt, boltzmann_free_energy, calibrate, calibrate_fit
__all__ += ["CalibrationReport", "basin_barriers_kt", "boltzmann_free_energy", "calibrate", "calibrate_fit"]

from .developmental import CommitmentReport, FreeEnergyProfile, commitment_profile, committor_free_energy_profile, developmental_coordinate
__all__ += ["CommitmentReport", "FreeEnergyProfile", "commitment_profile", "committor_free_energy_profile", "developmental_coordinate"]

from .api import AnalysisReport, analyze
__all__ += ["AnalysisReport", "analyze"]

from .entropy_production import (
    EntropyProductionConfig, EntropyProductionReport, PairFlux,
    compute_pair_fluxes, estimate_entropy_production, stationarity_residual,
    cross_validate_with_jarzynski,
)
__all__ += [
    "EntropyProductionConfig", "EntropyProductionReport", "PairFlux",
    "compute_pair_fluxes", "estimate_entropy_production", "stationarity_residual",
    "cross_validate_with_jarzynski",
]
