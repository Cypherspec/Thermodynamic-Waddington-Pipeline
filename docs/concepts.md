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

## Entropy production - the central quantity

For a system in equilibrium, forward and reverse transitions balance and no
entropy is produced. A driven system breaks that balance. Following Seifert's
stochastic thermodynamics, the pipeline estimates the entropy-production rate
from the net probability fluxes between cell states.

The key is the **null test**: the velocity assignment is shuffled across cells
many times to build a distribution of entropy production under "no real
direction". The observed value's permutation p-value is the test for departure
from equilibrium. A small p-value means the process is measurably irreversible;
a large one means it is not distinguishable from equilibrium. Always read the
p-value, not the raw rate.

## Schnakenberg cycles

The global entropy production can be decomposed exactly into contributions from
the fundamental cycles of the graph (Schnakenberg's network theory). This
attributes irreversibility to specific loops in state space.

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
the thermodynamic entropy-production rate.

## References

Waddington 1957; Jarzynski 1997; Schnakenberg 1976; Kemeny & Snell 1960;
Seifert 2012; E & Vanden-Eijnden 2006; La Manno et al. 2018; Bergen et al. 2020.
