from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class FitConfig:
    """Settings for a landscape fit run."""

    neighbors: int = 24
    dimensions: int = 6
    density_bandwidth: float = 0.65
    diffusion_floor: float = 1e-4
    temperature: float = 1.0
    velocity_scale: float = 1.0
    max_paths: int = 64
    bootstrap_replicates: int = 24
    bootstrap: int | None = None
    seed: int = 17
    reference_index: int = 0
    edge_alignment_threshold: float = 0.0
    barrier_quantile: float = 0.9
    regularization: float = 1e-6
    anisotropy_strength: float = 0.35
    current_strength: float = 0.5
    topology_persistence_threshold: float = 0.0
    trajectory_steps: int = 32
    trajectory_replicates: int = 64
    enable_counterfactuals: bool = True
    enable_continuation: bool = True
    enable_spectral_analysis: bool = True
    estimator: str = "effective_path_work"
    path_sampling: str = "edge_ensemble"
    diffusion_model: str = "local_scalar_residual"
    current_model: str = "directed_knn_proxy"
    bootstrap_edge_fraction: float = 0.8
    uncertainty_method: str = "edge_subsampling"
    audit_science: bool = True
    enable_report_payload: bool = True
    enable_entropy_production: bool = True
    enable_cycle_decomposition: bool = True
    entropy_production_bootstrap_replicates: int = 32
    entropy_production_permutation_replicates: int = 32

    def effective_bootstrap_replicates(self) -> int:
        return self.bootstrap if self.bootstrap is not None else self.bootstrap_replicates

    def validate(self, n_cells: int, n_features: int) -> None:
        if n_cells < 4:
            raise ValueError("need at least 4 cells")
        if n_features < 2:
            raise ValueError("need at least 2 features")
        if not 1 <= self.neighbors < n_cells:
            raise ValueError("neighbors must be between 1 and n_cells-1")
        if not 1 <= self.dimensions <= min(n_cells, n_features):
            raise ValueError("dimensions exceeds matrix rank")
        if self.density_bandwidth <= 0:
            raise ValueError("density_bandwidth must be positive")
        if self.temperature <= 0:
            raise ValueError("temperature must be positive")
        if self.velocity_scale <= 0:
            raise ValueError("velocity_scale must be positive")
        if self.max_paths < 1:
            raise ValueError("max_paths must be positive")
        if self.estimator not in {"effective_path_work", "density_only", "hybrid"}:
            raise ValueError("unknown estimator")
        if self.path_sampling not in {"edge_ensemble", "trajectory_ensemble", "protocol_ensemble"}:
            raise ValueError("unknown path_sampling")
        if self.diffusion_model not in {"local_scalar_residual", "local_diagonal_residual", "constant"}:
            raise ValueError("unknown diffusion_model")
        if not 0 < self.bootstrap_edge_fraction <= 1:
            raise ValueError("bootstrap_edge_fraction must be in (0, 1]")
        if self.uncertainty_method not in {"edge_subsampling", "none"}:
            raise ValueError("unknown uncertainty_method")
        if self.effective_bootstrap_replicates() < 1:
            raise ValueError("bootstrap replicates must be positive")
        if not 0 <= self.edge_alignment_threshold <= 1:
            raise ValueError("edge_alignment_threshold must be in [0, 1]")
        if self.trajectory_steps < 1 or self.trajectory_replicates < 1:
            raise ValueError("trajectory controls must be positive")
        if self.anisotropy_strength < 0 or self.current_strength < 0:
            raise ValueError("strength params must be non-negative")
