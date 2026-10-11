# Limitations and scope

A method is trustworthy in proportion to how clearly it states where it works and
where it does not. This page is that statement. None of it is hidden elsewhere.

## What is well supported

- **A calibrated irreversibility measure.** The Hodge cyclic fraction separates
  equilibrium from non-equilibrium ground truth (AUROC 1.00) and tracks a known
  drive monotonically, with a correct false-positive rate. This is the core,
  defensible methodological claim.
- **Real lineages are close to reversible.** On pancreas, gastrulation
  erythroid, and bone marrow erythroid, the cyclic fraction sits at or modestly
  above a matched gradient floor (pancreas +0.05, the others ~0).
- **The commitment barrier.** 1.6 kT [1.4, 1.9] on pancreas across 8 subsamples.
- **Committor ordering on a long staged lineage.** On Drosophila spermatogenesis
  the committor beats PC1 in 8/8 subsamples (0.78 vs 0.47).
- **The software.** Deterministic given a seed, 225 unit tests, lint-clean CI,
  reproducible in one command, numeric core pinned to bit-identical outputs.

## What is dataset-dependent

- **Committor ordering.** The committor recovers developmental order better than
  PCA pseudotime on pancreas but not on the gastrulation erythroid lineage, where
  plain PC1 wins (on Drosophila spermatogenesis it wins clearly). Do not present the committor as a universally better pseudotime;
  present it as a commitment coordinate tied to the thermodynamics.

## What did not survive calibration

- **"Differentiation is strongly irreversible."** An earlier version reported
  significant entropy production across three tissues from the Seifert
  permutation test. That test false-positives on equilibrium ground truth (FPR
  1.0), so the reading is withdrawn; the calibrated measure shows near-reversible
  lineages instead. The legacy test remains in the package for comparison.
- **A depth-confounded per-cell map.** The per-cell cyclic-energy map
  (`cyclic_flow_per_cell`) grows with velocity magnitude; on the fly testis it
  produced a 7x "peak" that tracked sequencing depth. Use the magnitude-invariant
  `cyclic_fraction_per_cell` against a local floor instead.
- **Detection power on real geometry.** Injecting a rotation into the real
  near-1D pancreas geometry barely moves the cyclic fraction, because a 1D
  structure has few cycles to carry circulation. Power is established on the 2D
  synthetic ground truth, not on these lineages.

## What did not replicate

- **Dentate gyrus neurogenesis.** The legacy entropy-production signal did not replicate
  across seven independent tests, including a preregistered confirmatory test.
  This negative result is retained in `RESEARCH_NOTE_entropy_production_pancreas.md`.
  It defines a boundary of the method, not a bug.

## What is not claimed

- **The free energy is not a molecular free energy.** It is a density-based
  pseudopotential expressed in kT. Rare or under-sampled states read as high
  energy. Read barriers as a pseudopotential, and use the calibrated cyclic
  fraction (not the legacy p-value) for non-equilibrium claims.
- **RNA velocity is an input, not a truth.** Results depend on the velocity
  model. Most benchmarks use an unspliced-minus-spliced proxy; one dynamical-
  velocity check is included. Broader velocity models remain to be tested.
- **Validation is not exhaustive.** Three vertebrate tissues, Drosophila
  spermatogenesis, synthetic ground truth, and one documented negative. The fly
  analysis uses one library at ~10% read depth, and its commitment-driver genes
  are a ranking for follow-up, not causal evidence. More datasets, and ideally wet-lab confirmation, are needed before
  strong biological conclusions.

## How to use it responsibly

- Report the cyclic fraction next to its matched gradient floor, not the legacy
  entropy-production p-value.
- Treat kT barriers as pseudopotential barriers.
- State the velocity model and preprocessing.
- If applying to a new system, compute the gradient floor at the same sample
  size first; if the observed cyclic fraction is not above it, the system is not
  measurably irreversible under this method, and that is a valid result to report.

## A preregistered test

A falsifiable prediction and its decision rule are written in advance in
`PREREGISTRATION_commitment.md`, so the commitment claim can be confirmed or
refuted by experiment rather than argued.
