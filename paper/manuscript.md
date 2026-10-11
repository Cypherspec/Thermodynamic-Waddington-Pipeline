# Non-equilibrium thermodynamics of single-cell fate: calibrated irreversibility and commitment from RNA velocity

**Nischay Kommisetty**

*Affiliation to be confirmed by the author before submission.*

Correspondence: nischay.kommisetty@gmail.com
Code: https://github.com/Cypherspec/Thermodynamic-Waddington-Pipeline

---

## Abstract

Trajectory-inference tools reconstruct the order and branching of single-cell
differentiation, but they do not quantify how far the process is from
thermodynamic equilibrium. We present a pipeline that treats differentiation as a
stochastic process on a kNN graph annotated with a Jarzynski path-work functional
derived from RNA velocity, and that computes stochastic-thermodynamic observables
directly: a calibrated irreversibility measure, Schnakenberg cycle affinities in
kT, an effective free-energy landscape, a transition-path-theory committor as a
commitment coordinate, and the potential of mean force along it. We first show
that a natural irreversibility test - a Seifert entropy-production rate with a
velocity-shuffle permutation null - is miscalibrated: on synthetic equilibrium
ground truth it is significant in every run (false-positive rate 1.0; AUROC 0.44
for equilibrium versus driven), because the shuffle destroys velocity-position
coupling rather than detailed balance. We replace it with a discrete Hodge
decomposition of the velocity flow into gradient (reversible) and cyclic
(irreversible) parts. The cyclic fraction is calibrated on the same ground truth
(AUROC 1.00, monotonic in a known drive, equilibrium floor 0.12). Applied to three
developmental lineages (pancreatic endocrinogenesis, gastrulation erythropoiesis,
bone marrow erythropoiesis), it shows they are close to reversible: pancreas
carries a modest excess of 0.05 over a matched gradient floor, the other two
essentially none, under both a velocity proxy and scVelo dynamical velocity. The
same holds in a second phylum: on Drosophila spermatogenesis (Fly Cell Atlas),
with velocity built from the raw reads, the germline is near-reversible, and the
committor orders it far better than principal-component pseudotime. On
pancreas the committor locates commitment at the Pre-endocrine stage with a
commitment barrier of 1.6 kT [1.4, 1.9]; we report that its ordering advantage
over principal-component pseudotime does not generalize to gastrulation, and that
CellRank orders the pancreas lineage better. We release a tested,
one-command-reproducible package and a preregistered falsifiable test. The
contribution is a correct thermodynamic readout of differentiation, including the
correction of an artifact that a naive analysis would report as a discovery.

## 1. Introduction

C. H. Waddington's epigenetic landscape is a metaphor of cells rolling downhill
into fates. Modern single-cell methods have made the landscape quantitative
through pseudotime and RNA velocity: scVelo (Bergen et al., 2020), CellRank
(Lange et al., 2022), Palantir (Setty et al., 2019), and Waddington-OT
(Schiebinger et al., 2019) estimate ordering, fate probabilities, and transport
between stages. These tools answer *where* a cell is going. They do not answer a
distinct, physically meaningful question: *how far from equilibrium is the
journey?* A process can be ordered yet close to equilibrium, or strongly driven
and dissipative; the difference is the entropy production, a central quantity of
stochastic thermodynamics (Seifert, 2012).

We bring the non-equilibrium toolkit to single-cell data, and find that the
measurement itself is the hard part. Because RNA velocity supplies a per-cell
drift, it is tempting to estimate entropy production from it and test against a
null that shuffles velocity across cells. We show that this test fires on
reversible flow, so a "significant" result is not evidence of irreversibility.
Our central claim is methodological and testable: the cyclic (non-gradient)
fraction of the velocity flow from a discrete Hodge decomposition is a
calibrated irreversibility measure, and under it real developmental lineages are
close to gradient-like. We add a transition-path-theory committor (E and
Vanden-Eijnden, 2006) as a principled commitment coordinate, and a potential of
mean force that turns the landscape into a commitment barrier in units of kT.

## 2. Methods

**Graph and path work.** Cells are embedded by PCA; a k-nearest-neighbor graph
is built with an exact KD-tree. Each velocity-aligned directed edge is annotated
with a path-work value combining a local density ratio, the velocity drift along
the edge, and a diffusion noise term, scaled by an effective temperature.

**Free-energy landscape.** Effective free energies are propagated by multi-source
Jarzynski exponential averaging, assigning every cell a value from the
minimum-work path in its weakly connected component. Attractor basins and barrier
heights follow.

**Calibrated irreversibility.** The velocity field induces an edge flow
f_ij = (v_i + v_j)/2 . (x_j - x_i) on the undirected kNN graph. Any edge flow
splits uniquely into a gradient component B psi, the flow of a node potential,
which is reversible, and an orthogonal cyclic remainder, which circulates and
breaks detailed balance (discrete Hodge decomposition). The potential is the
least-squares solution of B psi = f, solved sparsely in O(edges). The
irreversibility is the cyclic energy fraction ||f - B psi||^2 / ||f||^2, with an
edge-bootstrap confidence interval. Because a finite graph gives a small non-zero
fraction even for a pure gradient field, real data are compared with a matched
gradient floor on the same cells. A per-cell map assigns each cell the cyclic
flow energy on its incident edges.

**Cycle affinities.** Following Schnakenberg (1976), the thermodynamic affinity of
each fundamental cycle is computed in kT from a spanning-tree potential; it is
zero for a conservative field (Kolmogorov's criterion) and grows with
circulation.

**Legacy entropy production.** For comparison we retain Seifert's discrete-state
pair-flux entropy-production rate with a label-permutation null that shuffles
velocity across cells.

**Ground-truth validation.** A linear diffusion dx = -Bx dt + sqrt(2D) dW with
B = aI + wJ has zero entropy production at w = 0 (equilibrium) and entropy
production growing as w^2 otherwise. We sample 400 cells per field over 15 seeds
and drives w in {0, 0.25, 0.5, 1, 2, 4}, and score each method on separating
equilibrium from driven fields.

**Commitment coordinate.** Given progenitor (source) and terminal (target) cell
sets, the forward committor q(x) - the probability of reaching the terminal fate
before returning to the progenitor - is obtained by fixed-point iteration of the
harmonic equation on the directed graph, weighted by the landscape. The
potential of mean force F(q) = -kT ln P(q), estimated by kernel density along q,
gives the commitment barrier; its maximum is the transition state. Mean
first-passage times (Kemeny and Snell, 1960) give commitment timescales.

**Physical calibration.** The landscape is placed on a kT scale by a
density-anchored Boltzmann inversion, F = -kT ln P over a kernel density of the
embedding, the standard Waddington pseudopotential. We report the agreement
(R^2) between the path-work landscape and this equilibrium reference as a
diagnostic of non-equilibrium departure.

**Implementation.** The numeric core is vectorized (numpy) with an exact KD-tree
graph; results are bit-identical to a pure-Python reference. Graph construction
and both irreversibility measures run on 50,000 cells in seconds. The package is
typed, lint-clean, tested (230 unit tests), and reproducible in one command.

## 3. Results

**The naive entropy-production test is miscalibrated.** On the synthetic ground
truth, the Seifert permutation test is significant at alpha = 0.05 in every
equilibrium run (false-positive rate 1.0) and separates equilibrium from driven
fields no better than chance (AUROC 0.44). A direction-only velocity-coherence
score also fails (AUROC 0.53), since both kinds of field are coherent. On real
pancreas data the same test separates true velocity from shuffled velocity with
AUROC 1.00 - but so would any statistic sensitive to velocity-position coupling,
which the shuffle destroys whether or not detailed balance is broken.

**The cyclic fraction is calibrated.** On the same ground truth the cyclic
fraction is low and flat at equilibrium (0.12) and rises monotonically with the
drive (0.15, 0.25, 0.47, 0.72, 0.85 for w = 0.25 to 4), separating equilibrium
from driven fields with AUROC 1.00. Cycle affinities in kT are likewise monotonic
in the drive.

**Real lineages are close to reversible.** At 1,500 cells over five subsamples,
the cyclic fraction of real proxy velocity against a matched gradient floor is
0.151 vs 0.099 on pancreatic endocrinogenesis (Bastidas-Ponce et al., 2019;
excess +0.051), 0.080 vs 0.085 on gastrulation erythropoiesis (Pijuan-Sala et
al., 2019; -0.005), and 0.108 vs 0.097 on bone marrow erythropoiesis (Setty et al., 2019; +0.012),
all far below a shuffled ceiling of 0.29-0.35. On 800 pancreas cells, scVelo
dynamical velocity gives the same reading as the proxy (excess -0.001 vs -0.003).
A linear lineage has an arrow of time, carried by the committor, but little
circulation. We note a limit on this null: injecting a rotation into the real
near-1D pancreas geometry barely moves the fraction, because a 1D structure has
few cycles to carry circulation, so detection power rests on the synthetic
benchmark rather than on these geometries.

**A commitment coordinate and barrier.** On pancreas, the committor recovers the
lineage order unsupervised (Ductal to Ngn3 to Pre-endocrine to Beta), places the
commitment transition (q crosses 0.5) at the Pre-endocrine stage, and yields a
potential-of-mean-force barrier of 1.6 kT [1.4, 1.9] across eight subsamples, with
the transition state at q ~ 0.43. The density-anchored landscape is about 5 kT
deep, and its low agreement with the path-work landscape (R^2 ~ 0.14) quantifies
the flow structure a density landscape cannot see.

**An honest limit on ordering.** On pancreas the committor orders cells better
than principal-component pseudotime (Spearman 0.94 +/- 0.02 vs 0.89, winning
10/10 subsamples, paired p = 0.002). This advantage does not generalize: on the
gastrulation erythroid lineage, principal-component pseudotime wins (0.94 vs
0.81). In a like-for-like comparison at 2,000 pancreas cells, CellRank's
two-boundary absorption probability orders best (0.99), ahead of the committor
(0.95) and PC1 (0.94). We therefore do not claim the committor as a better
ordering method; its value is as a commitment coordinate tied to the
thermodynamics.

**A noise-matched null.** The matched gradient floor is noise-free, but
measured velocity is not, and noise looks cyclic. Adding measurement noise to
the ground truth (uniform, or uneven with half the cells 10x noisier), the
clean-floor excess false-positives in every equilibrium run. We therefore build
the null from the data: a reversible velocity, the local gradient of the
potential that best fits the work along each edge, plus the observed residuals
(velocity minus its neighbourhood mean) re-drawn among neighbours with random
signs; the method's own bias is estimated on reference datasets that are
reversible by construction and subtracted. Across 200 ground-truth runs the
false-positive rate at equilibrium falls from 1.00 to 0.00 at every noise level,
and pooled equilibrium-vs-driven AUROC rises from 0.76 to 0.96, with power
decreasing with noise as expected.

**A second phylum: Drosophila spermatogenesis.** The Fly Cell Atlas testis (Li
et al., 2022) provides no spliced/unspliced layers, so we built them from the
raw reads: 61.8M reads (about 10% of library FCA59) streamed from the unaligned
10x BAMs and quantified with kallisto | bustools (88.6% pseudoaligned), matched to
the atlas annotation by barcode (7,335 germline cells). The velocity passes a
biological check: the unspliced fraction falls from 0.19 in spermatocytes to 0.04
in spermatids, the known post-meiotic transcription shutdown. The cyclic fraction
is 0.097 against a 0.086 floor (excess +0.011). The committor orders the germline
stages with Spearman 0.78 +/- 0.05 against 0.47 for PC1, winning 8/8 subsamples,
and places commitment within the spermatocyte program in 7/8; the barrier is 2.1
kT [1.6, 2.6]. Genes that switch across q = 0.5 form a reproducible ranking (top-50
overlap 0.36 between independent fits against 0.011 by chance), but a regulator
list fixed in advance is not enriched in it. The analysis also exposed an
artifact in a per-cell map: cyclic flow *energy* tracks sequencing depth (Spearman
+0.84 across stages) and produced a spurious late-spermatocyte peak. A
magnitude-invariant per-cell cyclic fraction removes it. Against the clean
floor, post-meiotic spermatids still showed an excess (+0.12) that ran opposite
to depth; against the noise-matched null it is -0.010, and the germline as a
whole is -0.007 (z = -0.6): reversible within its own noise.

**Prior negative result.** The legacy entropy-production signal did not replicate
on dentate gyrus neurogenesis across seven independent tests, including a
preregistered confirmatory test; this is retained in the research note.

## 4. Discussion

The pipeline occupies a niche the trajectory tools leave open: it measures the
non-equilibrium character of differentiation rather than its geometry. Its
strongest result is a calibration result. A permutation test that shuffles
velocity looks like a test of irreversibility and reports striking significance
on real tissues, but it is equally significant on reversible ground truth. The
Hodge cyclic fraction is calibrated, and it reads real developmental lineages as
close to gradient-like. That reading is consistent with the picture of
differentiation as descent on a landscape, and it is also the reading a
miscalibrated test would have hidden.

Limitations are explicit. The free-energy magnitude is a density-based
pseudopotential, not a molecular free energy. The cyclic fraction is a property
of the velocity field on a finite graph, so it must be read against a matched
floor, and its power on near-1D real geometries is low. Most results use an
unspliced-minus-spliced velocity proxy, with one dynamical-velocity check;
broader velocity models remain to be tested. Validation spans three vertebrate
lineages, Drosophila spermatogenesis (one library at reduced depth), synthetic
ground truth, and one documented negative. The committor
ordering advantage is dataset-dependent. These are stated because a method that
reports where it fails is more useful than one that does not.

## 5. Reproducibility

All results reproduce from the released code with `python benchmarks/run_all.py`;
individual benchmarks and figures are listed in the README. The package installs
with `pip install -e .`, passes 230 unit tests, and pins its numerical core to
bit-identical outputs across versions. Public datasets are obtained through
scVelo. An end-to-end tutorial runs on synthetic data with no download.

## 6. Preregistration

We preregister a falsifiable test of the commitment prediction
(`PREREGISTRATION_commitment.md`): a perturbation applied at the transition
state (Pre-endocrine, q ~ 0.5) should change terminal-fate outcome more than the
same perturbation at the progenitor state, with a pre-specified effect size,
paired test, and stated falsification criterion.

## References

- Waddington, C. H. (1957). *The Strategy of the Genes.*
- Jarzynski, C. (1997). Nonequilibrium equality for free energy differences. *Phys. Rev. Lett.*
- Schnakenberg, J. (1976). Network theory of microscopic and macroscopic behavior of master equation systems. *Rev. Mod. Phys.*
- Kemeny, J. G., Snell, J. L. (1960). *Finite Markov Chains.*
- Seifert, U. (2012). Stochastic thermodynamics, fluctuation theorems and molecular machines. *Rep. Prog. Phys.*
- E, W., Vanden-Eijnden, E. (2006). Towards a theory of transition paths. *J. Stat. Phys.*
- La Manno, G. et al. (2018). RNA velocity of single cells. *Nature.*
- Bergen, V. et al. (2020). Generalizing RNA velocity to transient cell states (scVelo). *Nat. Biotechnol.*
- Lange, M. et al. (2022). CellRank for directed single-cell fate mapping. *Nat. Methods.*
- Setty, M. et al. (2019). Characterization of cell fate probabilities (Palantir). *Nat. Biotechnol.*
- Schiebinger, G. et al. (2019). Optimal-transport analysis of developmental trajectories (Waddington-OT). *Cell.*
- Bastidas-Ponce, A. et al. (2019). Comprehensive single-cell mRNA profiling of pancreatic endocrinogenesis. *Development.*
- Pijuan-Sala, B. et al. (2019). A single-cell molecular map of mouse gastrulation. *Nature.*
- Li, H. et al. (2022). Fly Cell Atlas: a single-nucleus transcriptomic atlas of the adult fruit fly. *Science.*
- Bray, N. L. et al. (2016). Near-optimal probabilistic RNA-seq quantification (kallisto). *Nat. Biotechnol.*
