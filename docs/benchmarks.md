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

| tissue | cyclic fraction | clean floor | excess over clean floor | noise-matched excess |
|---|---|---|---|---|
| pancreas endocrine | 0.151 | 0.099 | +0.051 | **-0.011** (z = -1.0) |
| gastrulation erythroid | 0.080 | 0.085 | -0.005 | **-0.005** (z = -0.5) |
| bone marrow erythroid | 0.108 | 0.097 | +0.012 | **-0.010** (z = -0.8) |
| fly testis germline | 0.097 | 0.086 | +0.011 | **-0.007** (z = -0.6) |

Against the clean floor pancreas looked like it carried a modest excess; against
the noise-matched null (below) it does not, and no lineage does: all four are
reversible within their own noise. A linear lineage has an arrow of time but little
circulation, so low cyclic entropy production is the expected, correct result. The
measure's power is established on the synthetic geometry above (AUROC 1.00); on
these near-1D real lineages, even injecting a rotation barely moves the fraction,
because a 1D structure has no cycles to carry circulation -- which is itself why
the real result sits at the floor. The same holds for the noise-matched readout:
injecting rotation into the pancreas geometry moves it only at the strongest
drive (+0.010 at w = 2). `benchmarks/irreversibility_real.py`.

![four tissues](https://raw.githubusercontent.com/Cypherspec/Thermodynamic-Waddington-Pipeline/master/figures/irreversibility_four_tissues.png)

The directional progression is captured by the committor. The strong "irreversible
across three tissues" reading of the old permutation test does not survive
calibration, and is corrected here. `benchmarks/generalization.py`.

## Noise-matched irreversibility (robust to sequencing depth)

The matched gradient floor is a clean radial field: it has no noise. Real
velocity does, and noise looks cyclic, so shallow-sequenced or noisy cells read
as circulating. `noise_matched_irreversibility` builds the null from the data
itself: a reversible velocity (the local gradient of the potential that best
explains the work along each edge) plus the observed residual noise (each cell's
velocity minus its neighbourhood mean), re-drawn from neighbouring cells with
random signs. Running the identical procedure on reference datasets that are
reversible by construction estimates the method's own bias, which is subtracted
(a parametric double bootstrap).

`benchmarks/noise_matched_ground_truth.py` adds measurement noise to the linear
diffusion ground truth -- none, uniform, or uneven (half the cells 10x noisier,
like shallow cells) -- over 5 drives x 4 noise settings x 10 seeds:

| | clean-floor excess | noise-matched excess |
|---|---|---|
| false-positive rate at equilibrium, any noise | **1.00** | **0.00** |
| AUROC equilibrium vs driven, pooled over noise | 0.76 | **0.96** |
| noisy equilibrium ranked above a clean weak drive | 1.00 | **0.00** |
| per-cell excess in the noisy half, at equilibrium | +0.46 | **-0.001** |

Power falls with noise, as it should: under heavy uneven noise a weak drive is
no longer detectable, while w >= 1 still is. `analyze()` reports
`irreversibility_noise_matched_excess` and `..._z`; `analyze_adata` writes the
per-cell version to `obs['tw_cyclic_excess']`.

On the fly germline it settles the spermatid question: the post-meiotic excess
of +0.12 over the clean floor becomes **-0.010** against the noise-matched null,
the per-stage excess no longer tracks depth (Spearman +0.06, was -0.77), and the
whole germline is **-0.007** (z = -0.6): reversible within its own noise.

![noise-matched ground truth](https://raw.githubusercontent.com/Cypherspec/Thermodynamic-Waddington-Pipeline/master/figures/noise_matched_ground_truth.png)

## Drosophila spermatogenesis (a second phylum)

`benchmarks/fly_testis.py` runs the same calibrated readouts on the male germline
of the Fly Cell Atlas testis (Li et al., *Science* 2022): spermatogonia, the
spermatocyte program (stages 0-7a), and elongating spermatids, 7,335 cells with
the atlas's own annotation as stage labels. The atlas ships no spliced/unspliced
layers, so `benchmarks/fly_testis_prepare.py` builds them from the raw reads: it
streams an evenly spaced sample (61.8M reads, about 10% of library FCA59) out of
the unaligned 10x BAMs at ArrayExpress E-MTAB-10519 over HTTP range requests,
and quantifies nascent and mature RNA with kallisto | bustools against Ensembl
BDGP6.46 (88.6% pseudoaligned). All 15,292 atlas cells in the library match by
barcode. Genes carry FlyBase IDs.

| readout | result |
|---|---|
| velocity sanity check | unspliced fraction **0.19 in spermatocytes -> 0.04 in spermatids** (5x drop): the known post-meiotic transcription shutdown, recovered from reads quantified here |
| calibrated irreversibility | cyclic fraction **0.097** vs reversible floor 0.086 (excess **+0.011**), shuffled ceiling 0.52; noise-matched excess **-0.007** (z = -0.6) -- reversible within its noise, like the vertebrate tissues |
| committor ordering vs stage | Spearman **0.78 +/- 0.05 vs PC1 0.47**, committor wins **8/8** subsamples |
| where commitment happens | inside the spermatocyte program in 7/8 subsamples (spermatocyte 4 to 7a) |
| commitment barrier | 2.1 kT [1.6, 2.6], 8 subsamples |
| commitment-driver genes | top-50 list overlap between independent fits **0.36 vs 0.011 by chance**; a regulator list fixed in advance (tMAC, tTAF, bam/bgcn, spermatid genes) is **not** enriched (0/14 in the top 5%) |

The committor result is the clearest gain: on this long, staged lineage it beats
PC1 decisively, unlike gastrulation where PC1 wins. The driver genes (FlyBase
links in `experiments/fly_testis.json`) are a reproducible ranking for follow-up,
not causal evidence; the canonical regulators switch on when the spermatocyte
program begins, not at its end where commitment sits.

**A second artifact, caught.** The per-cell cyclic-energy map first showed a
7x peak of circulation in late spermatocytes. That map is an energy, so it grows
with velocity magnitude, and spermatocytes are the deepest-sequenced cells: across
stages it tracks sequencing depth (Spearman +0.84). The new magnitude-invariant
`cyclic_fraction_per_cell`, read against a per-stage reversible floor on the same
cells, removes it: the pre-meiotic germline sits within +0.014 of reversible. The
residual post-meiotic excess (+0.12) runs opposite to depth (-0.77) and to the
unspliced fraction (-0.75): it appears exactly where the velocity signal vanishes.
The noise-matched null (above) confirms it is noise: against it the spermatids
read -0.010.

![fly testis](https://raw.githubusercontent.com/Cypherspec/Thermodynamic-Waddington-Pipeline/master/figures/fly_testis.png)

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
point: CellRank does fate mapping well, and this pipeline adds the calibrated
irreversibility measure and kT barrier it does not compute.

## Performance and scale

Vectorizing the core made fits **2.2x-6.2x faster** with **bit-identical**
outputs, and a KD-tree graph handles **50,000 cells in seconds**. The full
diagnostic fit remains super-linear past a few thousand cells; the graph core is
what scales. `benchmarks/runtime_scaling.py`, `benchmarks/scaling.py`.

The irreversibility measures scale too. The Hodge cyclic-fraction projection moved
from a dense pseudo-inverse (cubic in cells) to a sparse least-squares solve
(bit-identical output): **88x faster at 2,000 cells**, and the cyclic fraction and
cycle affinity both run on **50,000 cells in about six seconds**, where the dense
projection would take tens of hours. `benchmarks/irreversibility_scaling.py`.

## Dynamical velocity

The headline results were re-run on real scVelo dynamical velocity: the committor
ordering held (0.98 on pancreas), and one entropy-production run came back non-
significant under scVelo's Ms-smoothed representation specifically, which the
gene-count and normalization sweeps then bounded. `benchmarks/dynamical_velocity_validation.py`.

The calibrated measure gives the same answer under both velocity models: on 800
pancreas cells the cyclic-fraction excess over the gradient floor is -0.003 for
the proxy and -0.001 for scVelo dynamical velocity, so the near-reversible
reading is not an artifact of the proxy. `benchmarks/dynamical_cyclic.py`.
