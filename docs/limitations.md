# Limitations and scope

A method is trustworthy in proportion to how clearly it states where it works and
where it does not. This page is that statement. None of it is hidden elsewhere.

## What is well supported

- **Irreversibility detection.** Entropy production distinguishes directed
  differentiation from an equilibrium (velocity-shuffled) control on three
  independent datasets from three tissues (pancreas, gastrulation erythroid, bone
  marrow), at chance for direction-blind baselines. Robust to gene count and
  normalization. This is the core, defensible claim.
- **The software.** Deterministic given a seed, 203 unit tests, reproducible in
  one command, numeric core pinned to bit-identical outputs.

## What is dataset-dependent

- **Committor ordering.** The committor recovers developmental order better than
  PCA pseudotime on pancreas but not on the gastrulation erythroid lineage, where
  plain PC1 wins. Do not present the committor as a universally better pseudotime;
  present it as a commitment coordinate tied to the thermodynamics.

## What did not replicate

- **Dentate gyrus neurogenesis.** The entropy-production signal did not replicate
  across seven independent tests, including a preregistered confirmatory test.
  This negative result is retained in `RESEARCH_NOTE_entropy_production_pancreas.md`.
  It defines a boundary of the method, not a bug.

## What is not claimed

- **The free energy is not a molecular free energy.** It is a density-based
  pseudopotential expressed in kT. Rare or under-sampled states read as high
  energy. Read barriers as a pseudopotential and use the permutation p-value for
  non-equilibrium claims.
- **RNA velocity is an input, not a truth.** Results depend on the velocity
  model. Most benchmarks use an unspliced-minus-spliced proxy; one dynamical-
  velocity check is included. Broader velocity models remain to be tested.
- **Validation is not exhaustive.** Two positive datasets and one documented
  negative. More datasets, and ideally wet-lab confirmation, are needed before
  strong biological conclusions.

## How to use it responsibly

- Report the entropy-production p-value, not the raw rate.
- Treat kT barriers as pseudopotential barriers.
- State the velocity model and preprocessing.
- If applying to a new system, run the shuffled-velocity control first; if the
  observed value is not separated from the null, the system is not measurably
  irreversible under this method, and that is a valid result to report.

## A preregistered test

A falsifiable prediction and its decision rule are written in advance in
`PREREGISTRATION_commitment.md`, so the commitment claim can be confirmed or
refuted by experiment rather than argued.
