# FAQ

**Is this a replacement for scVelo or CellRank?**
No. Those estimate pseudotime, fate probabilities, and transport, and do it well.
This adds a complementary layer they do not have: it measures whether the process
is thermodynamically irreversible. Use them together.

**Do I need RNA velocity?**
Yes, in some form. Velocity supplies the direction that makes the entropy-
production and path-work computations meaningful. You can supply scVelo dynamical
velocity or a transparent unspliced-minus-spliced proxy. Without measured
velocity the velocity-dependent outputs are disclosed as unavailable.

**What does a significant entropy-production p-value mean?**
That the observed irreversibility is not reproducible by shuffling the velocity
assignment, i.e. the process is measurably out of equilibrium. It does not by
itself measure a physical dissipation rate in joules.

**Why is my landscape's raw free energy not a good pseudotime?**
Because it is not meant to be. The raw energy is the landscape; the commitment
coordinate is the committor. Use `developmental_coordinate`, not `fit.energies`,
for ordering.

**Why does the committor beat PCA pseudotime on one dataset but not another?**
Ordering advantage depends on the geometry of the lineage. On a clean linear
lineage, PCA already orders nearly perfectly and is hard to beat. See
[Limitations](limitations.md).

**How big a dataset can it handle?**
Graph construction scales to tens of thousands of cells in seconds. The full
diagnostic fit is practical to a few thousand cells; disable the cycle
decomposition (`enable_cycle_decomposition=False`) for larger runs.

**Is it deterministic?**
Yes, given a seed. Benchmarks report means and confidence intervals over seeded
random subsamples, so they reproduce exactly.

**How do I cite it?**
See [`CITATION.cff`](https://github.com/Cypherspec/Thermodynamic-Waddington-Pipeline/blob/main/CITATION.cff).

**Something failed on my data. What now?**
Check `report.warnings` from `analyze`. Common causes: too few cells (need at
least four), fewer than two features, or no cells matching your `source_labels` /
`target_labels`. Open an issue with the shapes and the config.
