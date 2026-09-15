# Thermodynamic Waddington Pipeline

Non-equilibrium thermodynamics of single-cell fate.

Trajectory tools tell you *where* a cell is going. This pipeline tells you
whether the journey can run backwards, by measuring the **entropy production** of
differentiation from RNA velocity, with a permutation test rather than a headline
number.

## Why it exists

Differentiation is a driven, non-equilibrium process. The standard single-cell
toolkit (scVelo, CellRank, Palantir, Waddington-OT) measures its *geometry* -
order, fate probabilities, transport. None of them measures its
*irreversibility*. A process can be well-ordered yet close to equilibrium, or
strongly driven and dissipative; the difference is the entropy production, a
central quantity of stochastic thermodynamics that this package estimates
directly.

## What you get

- **An irreversibility test.** A Seifert entropy-production rate with a
  label-permutation null: a p-value for "is this process out of equilibrium".
  Validated on two independent datasets.
- **A free-energy landscape** on a physical kT scale, from a Jarzynski path-work
  functional and a density-anchored Boltzmann inversion.
- **A commitment coordinate.** The transition-path-theory committor, which
  locates where fate commits and gives a commitment free-energy barrier in kT.
- **Cycle-resolved irreversibility** (Schnakenberg), MFPT commitment timescales,
  and optimal-transport chaining across stages.

## Thirty-second example

```python
from thermodynamic_waddington import analyze

report = analyze(expression, velocity, labels=labels,
                 source_labels=["Ductal"], target_labels=["Beta"])

report.is_irreversible          # True / False
report.entropy_production_pvalue
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

This package is deliberate about its domain. Its robust, cross-dataset result is
irreversibility detection. Its committor ordering advantage is real on some
datasets and not others. It has a documented negative result. All of that is in
[Limitations and scope](limitations.md), because a method you can trust is one
that tells you where it fails.
