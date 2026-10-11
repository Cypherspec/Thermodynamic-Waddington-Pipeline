# Concepts

The science behind each number, in plain terms.

## The Waddington landscape

C. H. Waddington pictured cells as balls rolling downhill into fates. The modern
version is quantitative: from single-cell data we estimate an effective
free-energy landscape whose valleys are attractor states and whose ridges are
barriers between them. This pipeline builds that landscape and, crucially, asks a
question the picture leaves out: is the ball actually being *pushed*, or could it
roll back up just as easily?

## RNA velocity

RNA velocity estimates, for each cell, the direction it is moving in gene-
expression space, from the ratio of unspliced to spliced transcripts. It turns a
static snapshot into a vector field of motion. The pipeline accepts velocity from
scVelo's dynamical model, or a transparent unspliced-minus-spliced proxy.

## Path work and the free-energy landscape

Cells are placed on a k-nearest-neighbor graph in PCA space. Each velocity-
aligned edge is annotated with a *path-work* value combining three physical
terms: a density ratio, the velocity drift along the edge, and a diffusion noise
term. Free energies are then propagated by multi-source Jarzynski exponential
averaging, so every cell receives a value from the minimum-work path in its
component. This is an *effective* landscape, not a molecular free energy.

## Irreversibility - the central quantity

For a system in equilibrium, forward and reverse transitions balance (detailed
balance) and no entropy is produced. A driven system breaks that balance.

**Why the naive test fails.** A natural estimator is Seifert's pair-flux entropy
production with a null built by shuffling velocity across cells. But a shuffle
destroys *all* velocity structure, not only broken detailed balance, so the test
also fires on a perfectly reversible gradient flow: on synthetic equilibrium
ground truth it reports significance every time. It detects velocity-position
coupling, not irreversibility. It is kept in the package for comparison only.

**The calibrated measure.** Any edge flow on a graph splits uniquely (discrete
Hodge decomposition) into a gradient part, which is the flow of a potential and is
reversible, and a cyclic part, which circulates and is what breaks detailed
balance. The **cyclic fraction** - the share of flow energy in the cyclic part -
is ~0 for a conservative field and rises with real circulation. It is validated
on a linear diffusion with a known rotational drive (AUROC 1.00, monotonic in the
drive). Compare it with a matched gradient floor at the same sample size, since a
finite kNN graph gives a small non-zero floor.

## Schnakenberg cycles

Schnakenberg's network theory attributes entropy production to the fundamental
cycles of the graph. The pipeline computes the thermodynamic affinity of every
cycle, in kT, from a spanning-tree potential: zero for a conservative field
(Kolmogorov's criterion), non-zero with circulation. It is a physical-units
companion to the cyclic fraction.

## The committor and commitment

The forward committor q(x) is the probability that a cell at state x reaches the
terminal fate before returning to the progenitor. It is the principled reaction
coordinate of transition-path theory: 0 at the progenitor, 1 at the terminal
fate. Where it crosses 0.5 is the point of commitment. The potential of mean
force along it, F(q) = -kT ln P(q), gives the **commitment barrier** in kT.

## Physical calibration

The landscape is placed on a kT scale by a density-anchored Boltzmann inversion,
F = -kT ln P over a kernel density of the embedding - the standard Waddington
pseudopotential. The agreement (R^2) between the path-work landscape and this
equilibrium reference is itself informative: a low value means the flow-based
landscape carries non-equilibrium structure the density landscape cannot see.

## Mean first-passage time

MFPT gives the expected number of graph steps for a cell to first reach an
attractor (Kemeny-Snell theory), a kinetic commitment timescale that complements
the thermodynamic irreversibility measures.

## References

Waddington 1957; Jarzynski 1997; Schnakenberg 1976; Kemeny & Snell 1960;
Seifert 2012; E & Vanden-Eijnden 2006; La Manno et al. 2018; Bergen et al. 2020.
