# Benchmarks and validation

Every number here reproduces with `python benchmarks/run_all.py`. Results are
reported straight, negatives included.

## Irreversibility detection (the robust result)

Can the method separate a directed differentiation from a velocity-shuffled
equilibrium control? Entropy production reaches **AUROC 1.00** while expression-
based baselines (first principal component, distance from progenitor) sit at
**0.50 (chance)**, because they never use velocity direction.

This **generalizes** across datasets: significant on pancreatic endocrinogenesis
(directed p ~ 0.005-0.01) and on an independent mouse gastrulation erythroid
lineage (directed p ~ 0.01 vs shuffled ~ 0.35). It is robust to gene count
(20-400) and to log normalization.

`benchmarks/irreversibility_detection.py`, `benchmarks/pancreas_entropy_validation.py`,
`benchmarks/second_dataset_validation.py`.

## Commitment barrier

The potential of mean force along the committor gives a commitment barrier of
**1.6 kT [1.4, 1.9]** on the pancreas branch, with the transition state mid-
trajectory. `benchmarks/commitment_barrier.py`.

## Physical calibration

The landscape is **~5 kT** deep on a kT scale; the low path-work-vs-density R^2
(~0.14) quantifies its non-equilibrium departure. `benchmarks/free_energy_calibration.py`.

## Developmental ordering (dataset-dependent)

On pancreas the committor beats PC1 and an endpoint-informed axis at recovering
developmental stage: **0.94 +/- 0.02 vs 0.89**, winning **10/10** subsamples,
paired Wilcoxon **p = 0.002**. On the cleaner gastrulation erythroid lineage,
**plain PC1 wins** (0.94 vs 0.81). The committor's ordering advantage is
therefore **not universal**; its value is as a commitment coordinate tied to the
thermodynamics. `benchmarks/predictive_ordering.py`,
`benchmarks/statistical_validation.py`.

## CellRank comparison

`benchmarks/cellrank_comparison.py` compares the committor with CellRank's fate
probability to the Beta terminal on identical proxy velocity, scored on ordering.
On the tested machine CellRank's GPCCA step could not run (CellRank 2.0 via
pygpcca is incompatible with numpy 2.x), so its number is absent and recorded as
such; the committor scored 0.957 and PC1 0.937 on that run. Run the script in a
numpy<2 environment for CellRank's value. The point stands regardless: CellRank
does fate mapping, and this pipeline adds the entropy-production test and kT
barrier it does not compute.

## Performance and scale

Vectorizing the core made fits **2.2x-6.2x faster** with **bit-identical**
outputs, and a KD-tree graph handles **50,000 cells in seconds**. The full
diagnostic fit remains super-linear past a few thousand cells; the graph core is
what scales. `benchmarks/runtime_scaling.py`, `benchmarks/scaling.py`.

## Dynamical velocity

The headline results were re-run on real scVelo dynamical velocity: the committor
ordering held (0.98 on pancreas), and one entropy-production run came back non-
significant under scVelo's Ms-smoothed representation specifically, which the
gene-count and normalization sweeps then bounded. `benchmarks/dynamical_velocity_validation.py`.
