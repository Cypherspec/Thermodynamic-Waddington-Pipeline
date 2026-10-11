# Thermodynamic Waddington Pipeline

Non-equilibrium thermodynamics of single-cell fate.

Trajectory tools tell you *where* a cell is going. This pipeline asks how far
from equilibrium the journey is, with an irreversibility measure that is
**calibrated on ground truth** (it does not false-positive on reversible flow), a
committor commitment coordinate, and a free-energy landscape in kT.

## Why it exists

Differentiation is a driven, non-equilibrium process. The standard single-cell
toolkit (scVelo, CellRank, Palantir, Waddington-OT) measures its *geometry* -
order, fate probabilities, transport. None of them measures its
*irreversibility*. A process can be well-ordered yet close to equilibrium, or
strongly driven and dissipative; the difference is the entropy production, a
central quantity of stochastic thermodynamics. Measuring it correctly from a
velocity field is harder than it looks: the naive estimator over-detects, and this
package replaces it with a calibrated one.

## What you get

- **A calibrated irreversibility measure.** A discrete Hodge decomposition splits
  the velocity flow into a gradient (reversible) and a cyclic (irreversible) part;
  the cyclic fraction is ~0 for reversible flow and rises with real circulation.
  AUROC 1.00 on equilibrium-vs-driven ground truth, where the naive
  entropy-production permutation test false-positives every time.
- **Robust to noise and sequencing depth.** A noise-matched null built from the
  data's own noise: false positives 100% -> 0% on noisy ground truth.
- **Validated across two phyla.** Pancreas, gastrulation and bone marrow
  (vertebrate) and Drosophila spermatogenesis (Fly Cell Atlas, velocity built
  from the raw reads, FlyBase-linked genes).
- **A free-energy landscape** on a physical kT scale, from a Jarzynski path-work
  functional and a density-anchored Boltzmann inversion.
- **A commitment coordinate.** The transition-path-theory committor, which
  locates where fate commits and gives a commitment free-energy barrier in kT.
- **Cycle-resolved irreversibility** (Schnakenberg cycle affinities in kT), MFPT
  commitment timescales, and optimal-transport chaining across stages.

## Thirty-second example

```python
from thermodynamic_waddington import analyze

report = analyze(expression, velocity, labels=labels,
                 source_labels=["Ductal"], target_labels=["Beta"])

report.irreversibility_cyclic_fraction   # ~0 reversible, higher = circulation
report.irreversibility_cycle_affinity_kt # the same in kT (Schnakenberg)
report.commitment_label         # where fate commits
report.commitment_barrier_kt    # how hard, in kT
```

## Where to go next

- [Installation](installation.md)
- [Getting started](getting-started.md) - a full walkthrough
- [Concepts](concepts.md) - the science, explained
- [API reference](api.md)
- [Benchmarks and validation](benchmarks.md) - what holds and how well
- [Limitations and scope](limitations.md) - where it works and where it does not

## Honesty note

This package is deliberate about its domain. Its methodological result is a
*correct* irreversibility readout: the calibrated measure shows real
developmental lineages are close to reversible (gradient-like), which corrects an
earlier reading of the naive entropy-production test. Its committor ordering
advantage is real on some datasets and not others. All of that is in
[Limitations and scope](limitations.md), because a method you can trust is one
that tells you where it fails.
