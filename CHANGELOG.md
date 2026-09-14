# Changelog

All notable changes to this project are documented here. Format loosely follows
Keep a Changelog. Versioning is semantic.

## [Unreleased]

### Added
- Physical free-energy calibration (`calibration.py`): density-anchored
  Boltzmann inversion expresses the landscape in kT, plus a fit of the path-work
  landscape to that reference (kT-per-work-unit factor and R^2). Exposed as
  `calibrate`, `calibrate_fit`, `boltzmann_free_energy`, `basin_barriers_kt`.
- `benchmarks/free_energy_calibration.py`: kT calibration on the pancreas branch
  (about 5 kT deep, path-vs-density R^2 ~ 0.14) with a figure.
- KD-tree `build_knn` (n log n, bit-identical to the old brute force); graph
  construction now scales to 50k cells in seconds.
- `enable_cycle_decomposition` flag to skip the combinatorial cycle work on
  large datasets while keeping the landscape and entropy production.
- `benchmarks/scaling.py`: graph and core-fit scaling, reported honestly (graph
  scales; the full fit is still super-linear past a few thousand cells).
- `benchmarks/predictive_ordering.py`: honest head-to-head on developmental
  ordering. Simple baselines (PC1) beat the thermodynamic signals; ordering is
  not what the pipeline is for, and the result is reported straight.
- `benchmarks/irreversibility_detection.py`: the benchmark the method is for.
  Entropy production separates directed differentiation from a velocity-shuffled
  control at AUROC 1.00 while the expression-based baselines are at chance.
- Developmental coordinate (`developmental.py`): the forward committor as a
  principled reaction coordinate, plus `commitment_profile` to order cell types
  and locate where the committor crosses 0.5. On the pancreas branch it recovers
  the lineage order and beats PC1 at ordering (Spearman 0.97 vs 0.89). Exposed as
  `developmental_coordinate` and `commitment_profile`.
- `benchmarks/statistical_validation.py`: multi-seed CIs, win-rate, and a paired
  Wilcoxon test. The committor beats the baselines 10/10 seeds (0.94 vs 0.89,
  p=0.002).
- `benchmarks/run_all.py` to reproduce every benchmark in one command, and
  `examples/tutorial.py` for an end-to-end run on synthetic data.
- Free-energy profile along the committor (`committor_free_energy_profile`): the
  potential of mean force, giving a commitment barrier of ~1.6 kT [1.4, 1.9].
- `benchmarks/commitment_barrier.py` and `benchmarks/second_dataset_validation.py`.

### Findings (reported straight)
- Entropy-production significance is robust to gene count (20-400) and to log
  normalization, and generalizes to a second independent dataset (gastrulation
  erythroid: directed p=0.01 vs shuffled 0.35).
- The committor's ordering advantage is dataset-dependent: it beats PC1 on
  pancreas but loses to PC1 on gastrulation erythroid (0.81 vs 0.94). Not a
  universal win.
- A variational (parametric) committor was prototyped and underperforms the
  exact iterative graph solve; not shipped.

### Fixed
- Iterative SCC decomposition (`topology/cycles.py`) so deep graphs no longer
  overflow the Python recursion limit.

## [0.2.0] - 2026-09-13

### Added
- Packaging for release: numpy runtime dependency, MIT license, classifiers,
  keywords, project URLs, and `real-data` / `viz` / `fast` / `dev` extras.
- `py.typed` marker and a single-sourced `__version__`.
- `benchmarks/runtime_scaling.py` and `benchmarks/plot_scaling.py`: reproducible
  runtime benchmark against the pure-Python baseline, with a figure.
- `benchmarks/pancreas_entropy_validation.py`: real-data entropy-production check
  on the endocrine branch with a shuffled-velocity negative control.
- `causal.counterfactual_ranking`: O(E) ranking used inside the fit.
- pytest config so `pytest` targets the package test suite directly.

### Changed
- Vectorized the numeric hotspots with numpy (spectral stationary distribution,
  counterfactual scan, provenance platform tag, gene-name canonicalization).
  Fit is 2.2x-6.2x faster over 120-400 cells.

### Guarantees
- The vectorization is output-preserving: energies max abs diff is 0.0 and
  attractors and rankings match the previous version exactly.

### Unchanged
- All scientific claims and their status, including the documented dentate gyrus
  negative result in `RESEARCH_NOTE_entropy_production_pancreas.md`.
- 180 unit tests pass.
