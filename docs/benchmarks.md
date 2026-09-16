# Benchmarks and validation

Every number here reproduces with `python benchmarks/run_all.py`. Results are
reported straight, negatives included.

## Irreversibility detection (the robust result)

Can the method separate a directed differentiation from a velocity-shuffled
equilibrium control? Entropy production reaches **AUROC 1.00** while expression-
based baselines (first principal component, distance from progenitor) sit at
**0.50 (chance)**, because they never use velocity direction.

This **generalizes across three tissues**: pancreatic endocrinogenesis,
mouse gastrulation erythroid, and human bone marrow. Recomputed from scratch at
**1,500 cells, 100 genes, 400 permutations, 5 seeds each**, directed entropy
production beats **every one of the 400 permutations in every seed** (p = 0.0025,
the permutation floor) on all three tissues, while the shuffled-velocity control
is never significant (best-case shuffled p = 0.14, 0.32, 0.17). It is also robust
to gene count (20-400) and to log normalization.

`benchmarks/generalization.py`, `benchmarks/irreversibility_detection.py`,
`benchmarks/pancreas_entropy_validation.py`, `benchmarks/second_dataset_validation.py`.

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

`benchmarks/cellrank_comparison.py` puts the committor and CellRank on the same
footing: both compute the probability of reaching the Beta terminal before falling
back to the Ductal progenitor (a two-boundary absorption for CellRank's
VelocityKernel + GPCCA, the transition-path committor here) on identical proxy
velocity, scored on recovering developmental stage. Run at **2,000 cells** (3
seeds), the regime GPCCA is built for:

| method | abs Spearman with stage |
|---|---|
| CellRank two-boundary absorption | **0.99 +/- 0.00** |
| TW committor | 0.95 +/- 0.00 |
| PC1 | 0.94 +/- 0.00 |

**CellRank wins**, and it should: on a linear lineage with enough cells its
absorption probability is exactly the tool for the job. The committor beats PC1
and tracks CellRank closely, but ordering was never this pipeline's claim. A
single terminal is degenerate on a linear lineage (CellRank returns probability 1
everywhere), so the comparison uses two boundaries; an earlier 150-cell,
single-terminal run gave a noisy, non-meaningful CellRank number and is
superseded. The script includes a numpy-2 compatibility shim (for pygpcca) and
runs CellRank single-process to avoid a Windows multiprocessing hang. The durable
point: CellRank does fate mapping well, and this pipeline adds the
entropy-production test and kT barrier it does not compute.

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
