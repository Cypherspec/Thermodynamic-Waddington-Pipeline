# Changelog

All notable changes to this project are documented here. Format loosely follows
Keep a Changelog. Versioning is semantic.

## [0.3.0] - 2026-09-14

### Added
- Calibrated irreversibility measure (`irreversibility.cyclic_irreversibility`): a
  discrete Hodge decomposition splits the velocity flow into a gradient
  (reversible) part and a cyclic (irreversible) part; the cyclic energy fraction
  is the irreversibility. Reported as `analyze().irreversibility_cyclic_fraction`.
  `benchmarks/synthetic_ground_truth.py` validates it on a linear diffusion with
  known entropy production: AUROC 1.0 separating equilibrium from non-equilibrium
  and monotonic recovery of the drive, versus AUROC 0.44 (below chance) and a 100%
  false-positive rate at equilibrium for the old density-based Seifert EP test.
- High-level `analyze()` API returning a compact, serializable `AnalysisReport`
  (irreversibility test, kT landscape depth, attractors, and the committor
  commitment coordinate and barrier when source/target labels are given).
- `analyze_adata()`: AnnData-native entry point that reads velocity/expression
  layers, runs the pipeline, and writes the committor and energy back into
  `adata.obs` (and the report into `adata.uns`).
- `plots` module: one-line `committor`, `landscape`, and `free_energy_profile`.
- A Jupyter tutorial notebook and a full mkdocs documentation site.
- `benchmarks/cellrank_comparison.py`: a fair head-to-head with CellRank on the
  same quantity (probability of reaching Beta before the Ductal progenitor: a
  two-boundary absorption for CellRank, the committor here). Runs CellRank under
  numpy 2 via a small pygpcca compatibility shim and single-process velocity
  graph. Result at 2,000 cells (3 seeds, the regime GPCCA is built for): CellRank
  0.99 > committor 0.95 > PC1 0.94. CellRank wins on this linear lineage, as it
  should; the committor beats PC1 and tracks CellRank closely. Reported straight:
  ordering is not this pipeline's claim, the entropy-production test and kT
  barrier are. (An earlier 150-cell single-terminal run was degenerate - one
  terminal on a linear lineage gives probability 1 everywhere - and is superseded
  by this two-boundary run.)
- Ruff lint/format configuration; version bumped to 0.3.0.
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

### Changed
- The CellFlux instrument design moved to its own repository
  (Cypherspec/cellflux); the README links to it. The software repo stays focused
  on the analysis pipeline.

### Findings (reported straight)
- Correction from the calibrated measure: the density-based Seifert EP permutation
  test over-detects. On synthetic ground truth it fires on conservative
  (equilibrium) fields (false-positive rate 1.0), because it scores against the
  density proxy rather than the stationary current, so its cross-tissue
  significance largely reflects velocity-position coupling, not broken detailed
  balance. Under the calibrated cyclic-fraction measure the three developmental
  lineages are close to gradient-like (reversible): pancreas carries a modest
  excess over the reversible floor, gastrulation and bone marrow essentially none.
  A linear lineage has an arrow of time but little circulation, so low cyclic
  entropy production is the expected, correct result. The committor still captures
  the directional progression; the strong "irreversible across three tissues"
  reading of the old EP test does not survive calibration.
- Entropy-production significance is robust to gene count (20-400) and to log
  normalization, and generalizes across three independent datasets from three
  tissues (pancreas, gastrulation erythroid, bone marrow). `generalization.py`
  now recomputes this from scratch at 1,500 cells, 100 genes, 1,000 permutations
  and 5 seeds per tissue (it previously redrew stored constants): directed
  entropy production beats every permutation in every seed (p=0.0010, the
  permutation floor) on all three, shuffled control never significant. Stronger
  than the earlier 100-permutation p=0.010 and fully reproduced, not stored.
- The committor's ordering advantage is dataset-dependent: it beats PC1 on
  pancreas but loses to PC1 on gastrulation erythroid (0.81 vs 0.94). Not a
  universal win.
- A variational (parametric) committor was prototyped and underperforms the
  exact iterative graph solve; not shipped.

### Fixed
- Iterative SCC decomposition (`topology/cycles.py`) so deep graphs no longer
  overflow the Python recursion limit.
- Diffusion-tensor estimation no longer materializes a dense genes-by-genes matrix
  per cell (it is scalar-diagonal and only its trace is used), so fit memory scales
  with cells instead of cells times genes^2. Outputs are identical; the full fit
  now runs at a few thousand cells without exhausting memory.

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
