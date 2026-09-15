# API reference

The public surface. Import everything from the top level:
`from thermodynamic_waddington import ...`.

## High level

### `analyze(expression, velocity, labels=None, source_labels=None, target_labels=None, config=None, significance=0.05) -> AnalysisReport`

Runs the pipeline end to end and returns the headline results. Pass `labels` plus
`source_labels` (progenitor) and `target_labels` (terminal) to also get the
committor commitment coordinate and barrier.

### `AnalysisReport`

Dataclass with: `n_cells`, `n_edges`, `n_attractors`, `entropy_production_rate`,
`entropy_production_pvalue`, `is_irreversible`, `landscape_range_kt`,
`path_vs_density_r2`, `committor_order`, `commitment_label`,
`commitment_barrier_kt`, `transition_state_q`, `warnings`. Methods: `to_dict()`,
`save(path)`.

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

### `estimate_entropy_production(points, velocity, densities, diffusions, graph, config) -> EntropyProductionReport`

Seifert entropy-production rate with bootstrap CIs and a permutation null
(`entropy_production_rate`, `permutation_p_value`, `bootstrap_ci`, `null_mean`).

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
