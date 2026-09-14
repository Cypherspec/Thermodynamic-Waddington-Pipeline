# Non-equilibrium thermodynamics of single-cell fate: detecting irreversibility and locating commitment from RNA velocity

**Nischay Kommisetty**

*Affiliation to be confirmed by the author before submission.*

Correspondence: nischay.kommisetty@gmail.com
Code: https://github.com/Cypherspec/Thermodynamic-Waddington-Pipeline

---

## Abstract

Trajectory-inference tools reconstruct the order and branching of single-cell
differentiation, but they do not quantify whether the process is
thermodynamically irreversible. We present a pipeline that treats differentiation
as a non-equilibrium stochastic process on a kNN graph annotated with a
Jarzynski path-work functional derived from RNA velocity, and that computes
stochastic-thermodynamic observables directly: a Seifert entropy-production rate
with a label-permutation null, a Schnakenberg cycle decomposition, an effective
free-energy landscape, a transition-path-theory committor as a commitment
coordinate, and the potential of mean force along it. Across two independent
datasets (pancreatic endocrinogenesis and mouse gastrulation erythropoiesis),
entropy production reliably distinguishes directed differentiation from a
velocity-shuffled equilibrium control (directed p ~ 0.01; shuffled not
significant), and does so where expression-based baselines are at chance
(AUROC 1.00 vs 0.50). This irreversibility signal is robust to gene count and
normalization. On pancreas, the committor recovers the lineage order, locates
commitment at the Pre-endocrine stage, and yields a commitment free-energy
barrier of 1.6 kT [1.4, 1.9]; we report honestly that the committor's ordering
advantage over principal-component pseudotime does not generalize to the
gastrulation lineage. We document a prior negative result on dentate gyrus
neurogenesis, provide a falsifiable preregistration of the commitment
prediction, and release a tested, one-command-reproducible package. The
contribution is a validated way to ask a question existing tools do not: is a
differentiation process out of equilibrium, and by how much?

## 1. Introduction

C. H. Waddington's epigenetic landscape is a metaphor of cells rolling downhill
into fates. Modern single-cell methods have made the landscape quantitative
through pseudotime and RNA velocity: scVelo (Bergen et al., 2020), CellRank
(Lange et al., 2022), Palantir (Setty et al., 2019), and Waddington-OT
(Schiebinger et al., 2019) estimate ordering, fate probabilities, and transport
between stages. These tools answer *where* a cell is going. They do not answer a
distinct, physically meaningful question: *is the journey irreversible?* A
process can be ordered yet close to equilibrium, or strongly driven and
dissipative; the difference is the entropy production, a central quantity of
stochastic thermodynamics (Seifert, 2012).

We bring the non-equilibrium toolkit to single-cell data. Building on the
observation that RNA velocity supplies a per-cell drift, we annotate a cell-state
graph with a path-work functional in the spirit of the Jarzynski relation
(Jarzynski, 1997), and compute stochastic-thermodynamic observables on it. Our
central, validated claim is narrow and testable: the entropy production of a
differentiating population is measurable and statistically distinguishable from
an equilibrium null. We add a transition-path-theory committor (E and
Vanden-Eijnden, 2006) as a principled commitment coordinate, connecting to a
line of work that has recently reached top venues for molecular systems (Nat.
Comput. Sci., 2025), and a potential of mean force that turns the landscape into
a commitment barrier in units of kT.

## 2. Methods

**Graph and path work.** Cells are embedded by PCA; a k-nearest-neighbor graph
is built with an exact KD-tree. Each velocity-aligned directed edge is annotated
with a path-work value combining a local density ratio, the velocity drift along
the edge, and a diffusion noise term, scaled by an effective temperature.

**Free-energy landscape.** Effective free energies are propagated by multi-source
Jarzynski exponential averaging, assigning every cell a value from the
minimum-work path in its weakly connected component (removing a single-reference
coverage artifact). Attractor basins and barrier heights follow.

**Entropy production.** We estimate the entropy-production rate with Seifert's
discrete-state pair-flux formula, with bootstrap confidence intervals and a
label-permutation null that shuffles the velocity assignment across cells. The
resulting permutation p-value is the test for departure from equilibrium. A
Schnakenberg cycle decomposition attributes the global rate to fundamental
graph cycles; velocity-divergence and gradient-alignment estimators provide
independent proxies.

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
scales to 50,000 cells in seconds. The package is typed, tested (193 unit
tests), and reproducible in one command.

## 3. Results

**Entropy production detects irreversibility, and generalizes.** On the
pancreatic endocrine branch (Bastidas-Ponce et al., 2019), observed entropy
production is significant against the shuffled-velocity null (p ~ 0.005-0.01)
while the control is not (p ~ 0.3-0.75). The same holds on an independent mouse
gastrulation erythroid lineage (Pijuan-Sala et al., 2019): directed p ~ 0.01,
shuffled ~ 0.35, across five subsamples. In a paired directed-vs-control
discrimination task, entropy production reaches AUROC 1.00 while
expression-based baselines (first principal component, distance from progenitor)
sit at chance (0.50), because they do not use velocity direction. The signal is
robust to gene count (20-400) and to log normalization; a single non-significant
run traced to scVelo's Ms-smoothed representation specifically.

**A commitment coordinate and barrier.** On pancreas, the committor recovers the
lineage order unsupervised (Ductal to Ngn3 to Pre-endocrine to Beta), places the
commitment transition (q crosses 0.5) at the Pre-endocrine stage, and yields a
potential-of-mean-force barrier of 1.6 kT [1.4, 1.9] across eight subsamples.

**An honest limit on ordering.** On pancreas the committor orders cells better
than principal-component pseudotime (Spearman 0.94 +/- 0.02 vs 0.89, winning
10/10 subsamples, paired p = 0.002). This advantage does not generalize: on the
gastrulation erythroid lineage, principal-component pseudotime wins (0.94 vs
0.81). We therefore do not claim the committor as a universally better ordering
method; its value is as a commitment coordinate tied to the thermodynamics.

**Prior negative result.** The entropy-production signal did not replicate on
dentate gyrus neurogenesis across seven independent tests, including a
preregistered confirmatory test; this is retained in the research note.

## 4. Discussion

The pipeline occupies a niche the trajectory tools leave open: it measures the
non-equilibrium character of differentiation rather than its geometry. Its
strongest, cross-dataset result is that differentiation is measurably
irreversible and that this is not reproduced by direction-blind baselines. The
committor connects our single-cell setting to transition-path theory, where
learning committors is an active area for molecular systems (Nat. Comput. Sci.,
2025); we use the exact discrete solve, which for a fixed cell graph outperforms
a parametric approximation we tested.

Limitations are explicit. The free-energy magnitude is a density-based
pseudopotential, not a molecular free energy. Results here use an
unspliced-minus-spliced velocity proxy and one dynamical-velocity check; broader
velocity models remain to be tested. Validation spans two positive datasets and
one documented negative; more are needed. The committor ordering advantage is
dataset-dependent. These are stated because a method that reports where it fails
is more useful than one that does not.

## 5. Reproducibility

All results reproduce from the released code with `python benchmarks/run_all.py`;
individual benchmarks and figures are listed in the README. The package installs
with `pip install -e .`, passes 193 unit tests, and pins its numerical core to
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
- (2025). Iterative variational learning of committor-consistent transition pathways. *Nat. Comput. Sci.*
