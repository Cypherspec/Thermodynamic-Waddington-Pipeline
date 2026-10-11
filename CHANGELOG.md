# Changelog

All notable changes to this project are documented here. Format loosely follows
Keep a Changelog. Versioning is semantic.

## [1.3.0] - 2026-10-11

### Added
- The paper: `paper/manuscript.md` rewritten as a full article (six figures,
  Table 1, Methods with equations, data/code availability, numbered
  references) and typeset to `paper/manuscript.pdf` by `paper/build.py`
  (pandoc + Typst, `paper/template.typ`).
- `benchmarks/run_all.py` runs the cross-tissue benchmark with the published
  1,000 permutations (it used 400, so the published p-values did not reproduce
  through it).

### Changed
- The noise-matched readout now runs on all four lineages
  (`benchmarks/irreversibility_real.py`, new figure
  `figures/irreversibility_four_tissues.png`): pancreas -0.011, gastrulation
  -0.005, bone marrow -0.010, fly testis -0.007 (all |z| <= 1). None circulates
  beyond its own noise.

### Fixed
- The docs, README and manuscript described pancreas as carrying a "modest
  excess" (+0.05) over the reversible floor. Against the noise-matched null that
  excess is -0.011: it was noise. Corrected everywhere; the near-1D detection-power
  caveat is kept ("no circulation detected", not "none").

### Added
- `benchmarks/plot_overview.py` regenerates the pipeline overview diagram.

## [1.2.0] - 2026-10-11

### Added
- `irreversibility.noise_matched_irreversibility`: the cyclic fraction against a
  null built from the data's own noise (reversible velocity from the
  best-fitting potential plus neighbour-redrawn, sign-flipped residuals),
  bias-corrected by a parametric double bootstrap on reversible references;
  global and per-cell. `analyze()` reports
  `irreversibility_noise_matched_excess` / `_z`; `analyze_adata` writes
  `obs['tw_cyclic_excess']`.
- `benchmarks/noise_matched_ground_truth.py`: ground truth with measurement
  noise (none / uniform / uneven): false-positive rate at equilibrium 1.00 for
  the clean-floor excess vs 0.00 noise-matched; pooled AUROC 0.76 -> 0.96.
- Fly benchmark reports the noise-matched readout: germline excess -0.007
  (z = -0.6); the post-meiotic +0.12 over the clean floor is -0.010, and the
  per-stage excess no longer tracks depth.

## [1.1.0] - 2026-10-11

### Added
- Drosophila spermatogenesis benchmark (`benchmarks/fly_testis.py`) on the Fly
  Cell Atlas testis germline, 7,335 cells: velocity sanity check (unspliced
  fraction drops 5x after meiosis), calibrated irreversibility (+0.011 over the
  reversible floor), committor ordering (beats PC1 in 8/8 subsamples, 0.78 vs
  0.47), commitment barrier, and reproducible FlyBase-linked commitment-driver
  genes with a pre-registered regulator check (not enriched; reported).
- `benchmarks/fly_testis_prepare.py`: builds spliced/unspliced counts from the
  Fly Cell Atlas raw reads (unaligned 10x BAMs at ArrayExpress E-MTAB-10519)
  by streaming HTTP byte ranges and quantifying with kallisto | bustools (nac),
  with barcode matching to the atlas annotation. Parser tests in
  `tests/test_fly_prepare.py`.
- `irreversibility.cyclic_fraction_per_cell`: magnitude-invariant local cyclic
  fraction in [0, 1]. `analyze_adata` writes it to
  `adata.obs['tw_cyclic_fraction']`.

### Fixed
- The per-cell cyclic-energy map (`cyclic_flow_per_cell`) grows with velocity
  magnitude and so can track sequencing depth; on the fly testis it produced a
  spurious 7x late-spermatocyte peak (Spearman +0.84 with depth). Documented,
  and superseded for comparisons by `cyclic_fraction_per_cell`.

## [1.0.0] - 2026-10-11

First stable release.

### Fixed
- `analysis.report` called `basin_membership` and `stochastic_paths` with the wrong
  signatures and raised on any call; it now runs (regression test added).
- Removed shadowed duplicate definitions of `trajectory.stochastic_paths` and
  `uncertainty.bootstrap` (the effective, later definitions are kept).
- `thermodynamic_waddington.EvidenceRecord` was silently exported from two modules;
  it stays the translation-gate record, and the decision-engine record is now
  exported as `DecisionEvidenceRecord`.
- 29 modules had their docstring after `from __future__ import annotations`, so
  `__doc__` was `None`; docstrings now come first.
- The recorded note in `experiments/irreversibility_real.json` claimed injected
  rotation raises the cyclic fraction on the real pancreas geometry; its own
  numbers show it barely moves. The note now matches the benchmark script.

### Changed
- Documentation site and manuscript draft rewritten around the calibrated
  irreversibility result (Hodge cyclic fraction), matching the README; the legacy
  Seifert permutation test is documented as over-detecting and kept for comparison.
- Codebase is ruff-clean (unused imports and variables removed, imports sorted);
  CI now runs lint and a wheel/sdist build check alongside the tests.
- CITATION and REPRODUCIBILITY brought up to date (version, bone marrow dataset).

## [0.3.0] - 2026-09-14

### Added
- Per-cell irreversibility map (`irreversibility.cyclic_flow_per_cell`): the cyclic
  (non-gradient) flow energy incident on each cell, so you can see *where* along a
  lineage circulation concentrates. `analyze_adata` writes it to
  `adata.obs['tw_irreversibility']` next to the committor and energy.
- The cyclic-fraction irreversibility measure now scales: the gradient projection
  uses a sparse least-squares solve instead of a dense n-by-n pseudo-inverse, so it
  runs in O(nnz) rather than O(cells^3). Bit-identical to the old result (max diff
  1e-16), ~46x faster at 1,500 cells, and 8,000 cells in under a second where the
  dense version took minutes. Falls back to the dense path if SciPy is absent.
- Schnakenberg cycle-affinity entropy production (`irreversibility.cycle_affinities`):
  the thermodynamic affinity of every fundamental cycle of the kNN graph, in kT,
  computed from a spanning-tree potential in O(edges). Zero for a conservative
  field (Kolmogorov's criterion), rising with real circulation; validated on the
  synthetic ground truth (monotonic in the known drive). Reported as
  `analyze().irreversibility_cycle_affinity_kt`, a physical-units companion to the
  cyclic fraction.
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
