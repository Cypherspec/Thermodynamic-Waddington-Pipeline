# API reference

The public surface. Import everything from the top level:
`from thermodynamic_waddington import ...`.

## High level

### `analyze(expression, velocity, labels=None, source_labels=None, target_labels=None, config=None, significance=0.05) -> AnalysisReport`

Runs the pipeline end to end and returns the headline results. Pass `labels` plus
`source_labels` (progenitor) and `target_labels` (terminal) to also get the
committor commitment coordinate and barrier.

### `analyze_adata(adata, velocity_layer="velocity", expression_layer=None, label_key=None, source=None, target=None, config=None, n_top_genes=50, key_added="tw") -> AnalysisReport`

The same, directly on an AnnData. Reads the velocity layer (or derives an
unspliced-minus-spliced proxy), uses the top-`n_top_genes` by variance, and writes
the per-cell committor, energy, and irreversibility map to
`adata.obs['{key_added}_committor']`, `adata.obs['{key_added}_energy']`, and
`adata.obs['{key_added}_irreversibility']`, plus the magnitude-invariant
`adata.obs['{key_added}_cyclic_fraction']` and the noise-matched
`adata.obs['{key_added}_cyclic_excess']`, with the report in
`adata.uns['{key_added}']`.

### `AnalysisReport`

Dataclass with: `n_cells`, `n_edges`, `n_attractors`,
`irreversibility_cyclic_fraction`, `irreversibility_cycle_affinity_kt`,
`irreversibility_noise_matched_excess`, `irreversibility_noise_matched_z`,
`landscape_range_kt`, `path_vs_density_r2`, `committor_order`,
`commitment_label`, `commitment_barrier_kt`, `transition_state_q`, `warnings`,
and the legacy `entropy_production_rate`, `entropy_production_pvalue`,
`is_irreversible` (kept for comparison; that test over-detects). Methods:
`to_dict()`, `save(path)`.

## Plotting

`from thermodynamic_waddington import plots` (needs matplotlib):

- `plots.committor(values, coords, ax=None)` - scatter colored by the committor.
- `plots.landscape(energies, coords, ax=None)` - scatter colored by free energy.
- `plots.free_energy_profile(values, temperature=1.0, ax=None)` - the commitment
  barrier curve.

`coords` is any n-by-2 array: a UMAP, a PCA embedding, or `fit.embedding`.

## Fitting

### `fit_landscape(expression, velocity, config=None, labels=None, ...) -> LandscapeFit`

The full pipeline. Returns a `LandscapeFit` with `energies`, `attractors`,
`edges`, `barriers`, `embedding`, and a large `diagnostics` dict (entropy
production, Schnakenberg, spectral, and more). `LandscapeFit.save/load` persist
it as JSON.

### `FitConfig`

All numerical settings: `neighbors`, `dimensions`, `temperature`,
`velocity_scale`, `max_paths`, bootstrap and permutation replicate counts,
`enable_cycle_decomposition`, `enable_entropy_production`, `seed`, and more.
`FitConfig().validate(n_cells, n_features)` checks inputs.

## Developmental coordinate

### `developmental_coordinate(fit, source_labels, target_labels, iterations=400) -> list[float]`

Per-cell committor: 0 at progenitor states, 1 at terminal fates.

### `commitment_profile(fit, source_labels, target_labels) -> CommitmentReport`

Orders cell types by mean committor and locates commitment (`order`,
`mean_committor`, `commitment_label`, `peak_flux_label`).

### `committor_free_energy_profile(committor, temperature=1.0, grid=60) -> FreeEnergyProfile`

Potential of mean force along the committor. Returns `q_grid`, `free_energy_kt`,
`barrier_kt`, `barrier_q`.

## Thermodynamics

### `irreversibility.cyclic_irreversibility(pts, vels, graph, bootstrap=0, seed=0) -> IrreversibilityReport`

Calibrated irreversibility: the cyclic (non-gradient) energy fraction of the
velocity flow from a discrete Hodge decomposition (`cyclic_fraction`,
`bootstrap_ci`, `n_pairs`, `n_cells`). Sparse least-squares, O(edges).

### `irreversibility.cycle_affinities(pts, vels, graph, temperature=1.0, velocity_scale=1.0) -> CycleAffinityReport`

Schnakenberg cycle affinities in kT (`n_cycles`, `rms_affinity_kt`,
`max_affinity_kt`, `mean_abs_affinity_kt`).

### `irreversibility.noise_matched_irreversibility(pts, vels, graph, draws=10, references=5, seed=0, per_cell=False)`

Cyclic fraction against a null carrying the data's own noise (reversible
velocity from the best-fitting potential, plus residuals re-drawn among
neighbours with random signs), bias-corrected on reversible-by-construction
references. Returns `NoiseMatchedReport` (`cyclic_fraction`, `null_mean`,
`bias`, `excess`, `z`, and `per_cell_excess` when `per_cell=True`). Use it
whenever velocity noise or sequencing depth varies across the data.

### `irreversibility.cyclic_fraction_per_cell(pts, vels, graph)`

Per-cell cyclic fraction in [0, 1]: on each cell's edges, cyclic energy over
cyclic plus gradient energy. Invariant to velocity magnitude, so regions with
different sequencing depth compare on one scale. Read it against the same
function applied to a gradient field on the same cells (the local floor).

### `irreversibility.cyclic_flow_per_cell(pts, vels, graph)`

Per-cell share of the cyclic flow energy. It also grows with velocity
magnitude, so it can track sequencing depth; prefer `cyclic_fraction_per_cell`
to compare regions.

### `estimate_entropy_production(points, velocity, densities, diffusions, graph, config) -> EntropyProductionReport`

Legacy Seifert entropy-production rate with bootstrap CIs and a permutation null
(`entropy_production_rate`, `permutation_p_value`, `bootstrap_ci`, `null_mean`).
Kept for comparison: it is significant on reversible ground truth, so do not use
it alone for irreversibility claims.

### `calibrate_fit(fit) -> CalibrationReport`  ·  `calibrate(path_energies, coords)`

Density-anchored kT calibration (`boltzmann_energy_range_kt`, `r_squared`,
`kt_per_work_unit`). `boltzmann_free_energy(coords)` returns the kT landscape;
`basin_barriers_kt(energies, labels)` gives per-type depths.

### `estimate_mfpt(edges, energies, attractors, temperature=1.0) -> MFPTResult`

Mean first-passage times and per-cell commitment times.

## Optimal transport

`sinkhorn`, `fit_schrodinger_chain`, `displacement_interpolate`,
`estimate_growth_weights` - work-cost optimal transport chained across stages.

## Data

`make_synthetic_dataset(cells, genes, seed) -> SyntheticDataset` for a labeled
synthetic dataset with expression, velocity, and lineage labels.

## Command line

Console scripts installed with the package include `twaddington` (fit / report /
audit), `tw-real` (real-data ingest and validation), and `tw-benchmark`. Run any
with `--help`.
