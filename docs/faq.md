# FAQ

**Is this a replacement for scVelo or CellRank?**
No. Those estimate pseudotime, fate probabilities, and transport, and do it well.
This adds a complementary layer they do not have: a calibrated measure of how
far from equilibrium the process is, and a commitment coordinate with a kT
barrier. Use them together.

**Do I need RNA velocity?**
Yes, in some form. Velocity supplies the direction that makes the
irreversibility and path-work computations meaningful. You can supply scVelo dynamical
velocity or a transparent unspliced-minus-spliced proxy. Without measured
velocity the velocity-dependent outputs are disclosed as unavailable.

**Which irreversibility number should I report?**
The calibrated cyclic fraction (`report.irreversibility_cyclic_fraction`),
compared with a matched gradient floor, and optionally the cycle affinity in kT.
Not the legacy entropy-production p-value: that test is significant even on
reversible ground truth, so a small p-value does not mean the process is out of
equilibrium. Neither number is a physical dissipation rate in joules.

**My real lineage has a cyclic fraction near the floor. Is the method broken?**
No. A linear lineage has an arrow of time (captured by the committor) but little
circulation, so a near-reversible reading is the expected result. The measure's
detection power is established on synthetic ground truth with a known drive.

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
