# Benchmarks and validation

Every number here reproduces with `python benchmarks/run_all.py`. Results are
reported straight, negatives included.

## Irreversibility, measured correctly

Two ways to ask "is the differentiation irreversible", with different answers.

**The naive test over-detects.** The density-based Seifert entropy-production
permutation test (shuffle velocity across cells) reaches AUROC 1.00 against that
control, but so does a one-line velocity-coherence score, because the shuffle
destroys all velocity structure, not just irreversibility. On a synthetic
equilibrium field whose true entropy production is zero it fires every time
(false-positive rate 1.0; AUROC 0.44 for equilibrium vs non-equilibrium). It is
detecting velocity-position coupling, not broken detailed balance.

**The calibrated measure.** `irreversibility.cyclic_irreversibility` splits the
velocity flow (discrete Hodge decomposition) into a gradient (reversible) part and
a cyclic (irreversible) part and reports the cyclic energy fraction. On the same
ground truth it is calibrated: an equilibrium floor of 0.12 rising monotonically
to 0.85 under a known rotational drive, **AUROC 1.00** separating equilibrium from
non-equilibrium where coherence sits at 0.53. `benchmarks/synthetic_ground_truth.py`.

**What the real tissues show (reported straight).** On real proxy velocity the
three developmental lineages are close to gradient-like:

| tissue | cyclic fraction | reversible floor | excess |
|---|---|---|---|
| pancreas endocrine | 0.151 | 0.099 | **+0.051** |
| gastrulation erythroid | 0.080 | 0.085 | -0.005 |
| bone marrow erythroid | 0.108 | 0.097 | +0.012 |

Only pancreas carries a modest excess over the reversible floor; gastrulation and
bone marrow essentially none. A linear lineage has an arrow of time but little
circulation, so low cyclic entropy production is the expected, correct result. The
measure's power is established on the synthetic geometry above (AUROC 1.00); on
these near-1D real lineages, even injecting a rotation barely moves the fraction,
because a 1D structure has no cycles to carry circulation -- which is itself why
the real result sits at the floor. `benchmarks/irreversibility_real.py`.

The directional progression is captured by the committor. The strong "irreversible
across three tissues" reading of the old permutation test does not survive
calibration, and is corrected here. `benchmarks/generalization.py`.

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
and tracks CellRank closely, but ordering was never this pipeline's claim. This
holds in the low-data regime too: at **150 cells** (5 seeds) the two-boundary
CellRank is still stable and ahead (0.99 +/- 0.00 vs committor 0.94 +/- 0.03 vs
PC1 0.92) -- so there is no small-n advantage for the committor either. A single
terminal is degenerate on a linear lineage (CellRank returns probability 1
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
