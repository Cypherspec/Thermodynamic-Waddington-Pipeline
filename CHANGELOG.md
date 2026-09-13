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
