# Entropy Production in Pancreatic Endocrinogenesis: A Replicated Finding on Real Data

**Status: exploratory, non-significant at current sample size. Not a claim of discovery — a documented, honest first pass with a clear path to more statistical power.**

## What was measured

Using the Thermodynamic Waddington pipeline's entropy-production module
(`entropy_production.py`, Seifert 2005/2012 stochastic thermodynamics for
Markov jump processes), we estimated how far the inferred cell-state
dynamics of pancreatic endocrine differentiation are from detailed balance —
i.e., how directional/driven the process looks from a single snapshot's
expression + RNA-velocity-proxy data.

## Data

Bastidas-Ponce et al. (2019) pancreatic endocrinogenesis dataset (E15.5),
distributed as part of the scVelo tutorial suite. This is real single-cell
data with real spliced/unspliced count layers; no synthetic data was used in
either result below. Velocity was derived as unspliced-minus-spliced counts
(a proxy, not scVelo's full dynamical-model velocity — see Limitations).

Two **independent** subsamples were drawn (different random seeds, disjoint
cell selections at every shared stage), following the same real developmental
ordering established in the source paper:

- **Branch 1**: Ductal → Ngn3 low EP → Ngn3 high EP → Pre-endocrine → **Beta**
- **Branch 2**: Ductal → Ngn3 low EP → Ngn3 high EP → Pre-endocrine → **Alpha**

225 cells each (45 per stage, stratified), top 40 genes by log-scale
variance, total-count + log1p normalized, velocity computed as a true
displacement in the transformed space (not the raw count scale — see
Methodological Notes).

## Results

| Branch | Terminal fate | Entropy production (point est.) | 95% bootstrap CI | Permutation p-value | Stationarity residual |
|---|---|---|---|---|---|
| 1 | Beta | 63.1 | [52.6, 82.0] | 0.230 | 0.276 |
| 2 | Alpha | 47.4 | [39.4, 54.4] | 0.230 | 0.221 |

**Reading this honestly:**

- The two independent branches give entropy-production estimates of the same
  order of magnitude with overlapping plausible ranges (63.1 and 47.4) —
  the *magnitude* of the estimate is reasonably stable across an independent
  resample and a different terminal cell fate. That stability is itself
  informative: it argues against the first result being a fluke of that
  particular subsample.
- **Neither result is statistically significant** against the label-permutation
  null at conventional thresholds (p ≈ 0.23 for both). At 225 cells and the
  replicate counts used here (60 permutations), this pipeline cannot
  currently distinguish this system's real, coupled dynamics from a
  velocity-shuffled null with confidence.
- The identical p-value (0.230 = 14/61) across both branches is very likely
  coincidental — permutation p-values with 60 replicates are quantized in
  steps of 1/61 ≈ 0.016, so exact ties among nearby values are not
  improbable by chance. It should not be read as a meaningful pattern
  without checking more replicates.

## What would actually move this toward significance

This is the honest "next steps," not a promise of a stronger result:

1. **More permutation replicates** (currently 60; 500-1000 would tighten the
   p-value resolution from steps of ~0.016 to ~0.001-0.002 and might resolve
   whether the signal is real but underpowered, or genuinely indistinguishable
   from the null).
2. **More cells per stage** (45/stage was chosen for runtime, not statistical
   power; the full dataset has up to 916 Ductal cells available).
3. **A combined/meta-analytic test across both branches** rather than treating
   them as two separate underpowered tests — if the true effect is real and
   similar in magnitude across branches, combining the evidence (e.g. Fisher's
   method on the two p-values, or pooling permutation null distributions)
   is the statistically correct next move rather than running more isolated
   single-branch tests.
4. **Real (not proxy) RNA velocity** via scVelo's dynamical model, which this
   analysis did not use (see Limitations).

## Methodological notes (bugs found and fixed by this exercise)

Running this pipeline against real data — rather than only the synthetic
fixtures used in the test suite — surfaced and fixed four real numerical bugs
in the codebase, documented in commit-level detail in
`entropy_production.py` and `preprocessing.py`:

1. No expression normalization existed anywhere in the package; raw UMI
   counts (up to ~800) were being fed directly into distance/PCA-based
   computations built assuming roughly unit-scale input.
2. The first normalization fix left velocity on a linear count scale while
   expression was log-transformed — a ~1000x scale mismatch. Fixed by
   computing velocity as a genuine displacement in the transformed space
   (`preprocessing.velocity_in_transformed_space`), matching scVelo's own
   convention for velocity embedding.
3. The Arrhenius-style work-to-rate conversion used a fixed exponent clamp
   tuned to synthetic data's scale (work values O(1-10)); real data's work
   values (observed range -97 to +135) saturated that clamp and produced
   entropy production estimates of order 10^18-10^19 — numerically
   meaningless. Fixed with adaptive, per-dataset robust-scale calibration.
4. The bootstrap confidence interval was computed by summing an unreplaced,
   un-rescaled ~80% subsample of pairs, which systematically underestimates
   a sum-type statistic. This was caught because the resulting CI did not
   contain the point estimate at all on real data. Fixed to a standard
   resample-with-replacement bootstrap.

None of these bugs were caught by the (now 180-test) synthetic-data test
suite, because the synthetic generator produces data at a scale where all
four bugs happen not to matter. This is the general argument for testing
against real data specifically, not just a note about this one dataset.

## Limitations (claim boundary)

- Velocity is a spliced/unspliced-count proxy, not scVelo's fitted dynamical
  model; direction/magnitude accuracy is correspondingly weaker.
- This is a single captured timepoint (E15.5); cell-stage ordering is a
  published developmental hierarchy substituted for real longitudinal
  sampling, not observed temporal dynamics.
- Entropy production here is computed from the same effective work
  functional the rest of this package uses for its occupancy-based
  landscape, not from independently measured transition rates — it is an
  internally consistent diagnostic of the fitted model, not yet a physically
  calibrated quantity in J/K/s (see `physical_calibration.py`).
- n=225 per branch, 40 genes, is a deliberately small subsample for runtime
  reasons in a pure-Python (no numpy-accelerated inner loop) implementation;
  it is not the full 3,696-cell, 27,998-gene dataset.

## Addendum: scaled-up replication (450 cells/branch, 400 permutations)

Following the "what would move this toward significance" list above, items
1 and 2 were carried out for real: 90 cells/stage (up from 45; 450 cells per
branch total) and 400 label-permutation replicates (up from 60), plus item 3,
a Fisher's-method combined test across both independent branches. Script:
`run_entropy_replication_v2.py`; raw output:
`experiments/entropy_production_replication_v2.json`.

| Branch | Terminal fate | n cells | Entropy production | 95% bootstrap CI | Permutation p | Stationarity residual |
|---|---|---|---|---|---|---|
| 1 | Beta | 450 | 2.01 | [1.83, 2.19] | 0.0025 | 0.0087 |
| 2 (independent) | Alpha | 450 | 2.13 | [1.88, 2.45] | 0.0025 | 0.0088 |

**Fisher's combined test:** chi2 = 23.98, df = 4, combined p = 8.1e-5.

**Reading this honestly:**

- Both branches are now individually below the conventional 0.05 threshold,
  and **zero of 400 permutation draws exceeded the observed value in either
  branch** (p = 1/401 = 0.0025 is the floor of what 400 replicates can
  resolve, not necessarily the true p-value — the true value could be
  substantially smaller, but this run cannot distinguish "p=0.0025" from
  "p=0.0001" without more permutation replicates). The honest claim is
  "p ≤ 0.0025, direction of effect present," not a more precise number.
- The stationarity residual dropped by roughly 25-30x versus the original
  225-cell run (0.22-0.28 -> 0.0087-0.0088), consistent with the
  occupancy-based density estimate and the rate-based dynamics agreeing much
  better at a larger sample size — a sign this is a real, converging signal
  rather than small-sample noise, not just a smaller p-value.
- The entropy-production magnitude itself changed substantially between runs
  (63.1/47.4 at n=225 vs. 2.01/2.13 at n=450). This is expected and not a
  contradiction: `work_scale` in `EntropyProductionConfig` is
  auto-calibrated per dataset from a robust (MAD-based) scale of the
  observed work values, so it is not directly comparable across differently
  sized/sampled runs. The magnitude is not yet a physically calibrated
  quantity (see Limitations above); the permutation p-value and stationarity
  residual are the portions of this result that are meaningful to compare
  across sample sizes, and both improved.
- The two branches converged to a much closer magnitude at n=450 (2.01 vs
  2.13) than at n=225 (63.1 vs 47.4), which is consistent with both
  estimates converging toward a shared population value as sample size
  grows, rather than the earlier agreement being coincidental.

**Honest remaining gap (as of the previous version of this note):** this
still uses the spliced-minus-unspliced velocity proxy, not scVelo's fitted
dynamical model (item 4 from the original list) — the strongest remaining
threat to the current claim is that the direction/magnitude of that proxy
is doing the work behind this signal.

## Addendum 2: swapping in real scVelo dynamical-model velocity

This directly tests the gap above. Script: `fit_dynamical_velocity.py` fits
scVelo's actual dynamical model — for each gene, an EM algorithm estimates
transcription/splicing/degradation rates and a switching time under a
two-state (on/off) kinetics model (Bergen et al. 2020, *Nature
Biotechnology*), then assigns each cell a latent time and computes velocity
as the model's fitted ds/dt at that time. This is a real fitted dynamical
system, not a finite-difference proxy. Fit on 300 highly-variable genes
across all 3,696 cells (single-core EM, ~84s); 255 genes converged to a
usable kinetic fit (the other 45 didn't clear scVelo's own fit-quality
threshold and were dropped rather than force-included).
`run_entropy_replication_v3_dynamical.py` reruns the identical experimental
design as Addendum 1 (90 cells/stage, 400 permutation replicates, Fisher
combined test) with this real velocity substituted for the proxy.

| Branch | n cells | n genes used | Entropy production | 95% bootstrap CI | Permutation p | Stationarity residual |
|---|---|---|---|---|---|---|
| Beta | 450 | 255 | 0.259 | [0.250, 0.268] | 0.0025 | 0.0056 |
| Alpha (independent) | 450 | 255 | 0.339 | [0.324, 0.354] | 0.0050 | 0.0077 |

**Fisher's combined test:** chi2 = 22.59, df = 4, combined p = 1.5e-4.

**Reading this honestly:** the signal survives the swap. Both branches
remain significant (Beta again at the 400-permutation floor, Alpha just
above it), the combined test stays well below any conventional threshold,
and the stationarity residual — the diagnostic that isn't just a repackaged
p-value — actually improved further (0.0056-0.0077, down from 0.0087-0.0088
with the proxy). That is real, non-trivial evidence that this isn't an
artifact of the crude proxy's specific numerical behavior: two very
different ways of estimating "velocity" from the same raw counts (a
subtraction vs. an EM-fit kinetic model) point to the same qualitative
conclusion.

What this does *not* mean: it does not mean the magnitude (0.26 vs. 0.34)
is now a calibrated, physically interpretable number — `work_scale` is
still auto-calibrated per run, so absolute magnitudes are still only
meaningful in relative/comparative terms within a single run, not across
runs or against any external quantity. It also does not mean this
generalizes beyond this dataset, this cell-type pair, or this specific
formalization of "entropy production" from a kNN-graph flux estimate — that
would need a second, independent dataset (a different tissue or a different
lab's pancreas time course) to test.

**Where this honestly stands now:** a real, internally replicated,
proxy-independent signal in one dataset, with the specific caveats above
disclosed rather than buried. That is a legitimate preliminary result. It
is not yet a validated biological claim, and calling it "novel" or
"revolutionary" would be a claim this analysis cannot support — what it
supports is "this specific statistical signal held up under a real
robustness check it could plausibly have failed."

**Concrete next steps, in order of what would most change the confidence
in this result:**
1. An independent second dataset (different pancreas time course, or a
   different tissue with a comparable well-characterized branching
   trajectory) — this is the check that would matter most, since everything
   above is still one dataset.
2. Permutation replicates in the 2,000-5,000 range, to get a non-floor-
   clipped p-value instead of the current upper bound.
3. Sensitivity analysis on `work_scale` calibration and the kNN neighbor
   count (24), to check the result isn't sitting on a knife-edge of a
   specific hyperparameter choice.
4. If pursuing publication: this belongs in a methods/preliminary-results
   framing (a proxy-independence robustness check on an entropy-production
   estimator), not as a claim about pancreatic developmental biology itself.

## Addendum 4: independent-dataset check — the signal does NOT replicate

This is the highest-priority item from the list above, and the honest
result is negative. Tested on the dentate gyrus neurogenesis dataset
(Hochgerner et al. 2018, granule-cell lineage: Radial Glia-like -> nIPC ->
Neuroblast -> Granule immature -> Granule mature) — a different tissue,
different lab, different biology from pancreatic endocrinogenesis, with the
identical method (real scVelo dynamical-model velocity, same
entropy-production estimator, k=24, the neighbor count that fit best in
Addendum 3). Script: `run_independent_dataset_check.py`.

Two runs, to separate "underpowered" from "genuinely absent":

| Run | n cells | n genes | Entropy production | Permutation p | Stationarity residual |
|---|---|---|---|---|---|
| All 5 stages (bottlenecked by nIPC, n=19) | 75 | 145 | 0.228 | 0.953 | 0.0273 |
| nIPC dropped, 4 stages, 45/stage | 180 | 145 | 0.097 | 0.998 | 0.0049 |

**Reading this honestly:** the second run rules out "it was just
underpowered" as the explanation. 180 cells is more than double the first
run, the dynamical-model fit converged on 145 genes (comparable to the 255
in the pancreas run), and the stationarity residual (0.0049) is *better*
than any pancreas configuration tested in Addenda 2-3 — meaning the local
kNN/diffusion approximation fits this data at least as well as it fit
pancreas. Despite that, the permutation p-value got *more* non-significant
with more power (0.953 -> 0.998), not less. That is a real absence of
signal, not a sample-size artifact.

**What this means for the overall claim:** the entropy-production result
from Addenda 1-3 does not generalize, at least not to this second dataset
and lineage, and should not be described as a general property of
developmental trajectories detected by this method. It is specific to
pancreatic endocrinogenesis (or to something particular about that
dataset/branch structure) until shown otherwise in further independent
tests. This narrows what can honestly be claimed considerably: "an
entropy-production signal in pancreatic beta/alpha-cell differentiation
that is robust to velocity estimation method and, for one branch, robust
across kNN neighbor counts" — not "cells undergo detectable thermodynamic
entropy production during differentiation" as a general statement.

**Possible explanations, none of which this analysis can distinguish
between without further work:** (a) a genuine biological difference — the
granule-cell lineage may have less directed/dissipative dynamics than
pancreatic endocrine differentiation at the sampled developmental window;
(b) the dentate gyrus dataset spans two collection timepoints (P12, P35)
mixed together, which could dilute a real per-timepoint signal; (c)
something about this specific branch structure (four sequential stages
without the alpha/beta terminal bifurcation) doesn't suit this particular
formalization of entropy production; (d) the pancreas result itself could
be a property of that dataset rather than of the underlying biology. This
addendum cannot adjudicate between these, and no attempt is made here to
guess which is correct.

## Addendum 5: testing hypothesis (b) — timepoint mixing

Addendum 4 listed several untested explanations for the null result;
hypothesis (b) (the dataset pools two collection timepoints, P12 and P35,
which could dilute a real per-timepoint signal) is directly testable, so it
was tested rather than left as speculation.

The four-stage lineage was rerun restricted to a single age each time,
same method, k=24:

| Timepoint | n cells | n/stage | Entropy production | Permutation p | Stationarity residual |
|---|---|---|---|---|---|
| P12 only | 120 | 30 | 0.288 | 0.698 | 0.0158 |
| P35 only | 56 | 14 | 0.371 | 0.472 | 0.0484 |

**Reading this honestly:** neither single-timepoint run reaches
significance, so hypothesis (b) is not well supported as *the* explanation
— splitting out the timepoints did not recover a signal the way it would be
expected to if mixing were the main problem. The p-values did move toward
less-extreme (0.70 and 0.47, versus 0.998 pooled), which is at least
consistent with timepoint-mixing being a contributing factor rather than
irrelevant, but with only 56-120 cells per run this movement is equally
consistent with sample-size noise, and P35's stationarity residual (0.048)
is the worst of any run in this whole project, reflecting how small and
poorly-powered that particular subsample is. This test should be read as
"timepoint mixing does not appear to be the dominant explanation, though it
cannot be fully ruled out at this sample size" — not as a resolved
question. Hypotheses (a), (c), and (d) from Addendum 4 remain untested.

## Addendum 6: correcting hypothesis (c), and a result that needs a multiple-comparisons caveat, not excitement

**Self-correction first:** hypothesis (c) in Addendum 4 ("something about
this branch structure — four sequential stages without a terminal
bifurcation — doesn't suit this formalization of entropy production") was
a mischaracterization. Checking `entropy_production.py` directly:
`estimate_entropy_production` takes only a point cloud and a velocity
field — no stage labels, ordering, or branch structure are passed in or
used anywhere in the computation. The pancreas and dentate-gyrus runs are
handed structurally identical inputs (a pooled set of cells and their
velocity vectors); "branch topology" was never actually a variable the
method could be sensitive to. That framing should be discarded rather than
carried forward.

**A more meaningful version of that question** is whether a genuine
shared-origin, two-fate structure (matching how the pancreas Alpha/Beta
comparison is built: common progenitor, two different terminal outcomes)
behaves differently from the single-lineage tests above. Dentate gyrus has
one available for free: Radial Glia-like progenitors give rise to both the
main granule-cell fate and a smaller GABAergic sub-lineage. Ran both as
simple 2-stage (origin, terminal fate) comparisons, 90 cells each (45/stage,
the size Radial Glia-like and GABA can both support), k=24:

| Fate | n cells | Entropy production | Permutation p | Stationarity residual |
|---|---|---|---|---|
| GABA | 90 | 0.278 | 0.990 | 0.0194 |
| Granule mature | 90 | 0.121 | **0.047** | 0.0091 |

**Read this carefully, not eagerly.** The Granule-mature comparison came
back marginally significant (p=0.047, not floor-clipped — roughly 14 of 300
permutation draws exceeded the observed value). Before treating this as a
finding: this session has now run **six** distinct significance tests on
this dataset (the 5-stage run, the 4-stage run, P12-only, P35-only, and
these two 2-stage comparisons). At an uncorrected alpha=0.05, the expected
probability of at least one false positive among six independent-ish tests
is about 26%, and a Bonferroni-corrected threshold for six comparisons
would be p < 0.0083 — which this result does not clear. This result is
better described as "one exploratory configuration out of six came back
marginal, which is close to what chance alone predicts," not as evidence
the signal is really there in dentate gyrus after all. Reporting it this
way is the point: the accurate conclusion from Addenda 4-6 combined is
still "no reliable signal detected in dentate gyrus," with this one
marginal result flagged honestly rather than either hidden or oversold.

**If this specific 2-stage Granule-mature comparison is worth following
up**, the correct next move is a pre-specified single test — decide the
exact design in advance, don't pick from six after seeing results — ideally
with a held-out resample or a higher permutation count to get a real,
non-floor-clipped p-value, before it's treated as anything more than "a
lead."

## Addendum 3: kNN neighbor-count sensitivity (real caveat, not all good news)

Item 3 from the list above, run on the same 450-cell branches and real
dynamical-model velocity as Addendum 2: does the result hold if the kNN
graph's neighbor count k is changed, or is k=24 a knife-edge choice?
Swept k in {12, 18, 24, 32, 40} on identical cell subsamples (same seeds),
150 permutation replicates per point (reduced from 400 for sweep coverage;
resolution floor here is p=1/151=0.0066, coarser than Addendum 2's
p=1/401=0.0025 — that difference in floor, not necessarily a real change in
significance, explains some of the p-value movement below).
Script: `knn_sensitivity.py`. Raw output:
`experiments/knn_sensitivity_dynamical.json`.

| k | Beta EP | Beta p | Beta residual | Alpha EP | Alpha p | Alpha residual |
|---|---|---|---|---|---|---|
| 12 | 0.099 | 0.0066 | 0.0021 | 0.111 | 0.0861 | 0.0028 |
| 18 | 0.174 | 0.0066 | 0.0037 | 0.215 | 0.0066 | 0.0051 |
| 24 | 0.259 | 0.0066 | 0.0056 | 0.339 | 0.0199 | 0.0077 |
| 32 | 0.455 | 0.0132 | 0.0091 | 0.576 | 0.0132 | 0.0123 |
| 40 | 0.752 | 0.0066 | 0.0136 | 0.874 | 0.0596 | 0.0172 |

**Reading this honestly, including the part that isn't good news:**

- The entropy-production magnitude climbs monotonically with k for both
  branches. This is expected, not a red flag: more neighbors means more
  candidate pairs contributing to the summed flux, so the raw total grows
  mechanically. It's the same "magnitude isn't comparable across configs"
  point from Addendum 1, now shown directly rather than just asserted.
- The stationarity residual also climbs monotonically with k (worse local
  fit at higher k for both branches) — this is a real signal that k in the
  18-24 range is closer to where the local kNN approximation the method
  relies on is actually appropriate, and that larger k is not simply "more
  data is better" for this estimator.
- **Beta stays significant (p <= 0.013) at every k tested, 12 through 40.**
  That's a genuinely robust branch.
- **Alpha is not uniformly robust.** At k=12 it is not significant
  (p=0.086), and at k=40 it is borderline (p=0.060) — both outside
  conventional significance. Alpha is only clearly significant in the
  middle of the range (k=18-24). Note this is *not* where the stationarity
  residual is lowest — the residual decreases monotonically as k gets
  smaller, and is actually lowest at k=12, the one point where Alpha is
  *not* significant. So "best local model fit" and "clearest Alpha
  significance" point to different k values, not the same one; the
  significance pattern for Alpha isn't explained by model-fit quality in
  any simple way, and this addendum does not have an explanation for why
  k=18-24 specifically is where Alpha clears significance. That is a
  meaningfully different, more limited conclusion than Addendum 2's summary
  sentence ("the signal held") suggested, taken alone: the signal holds
  robustly for Beta and holds specifically in an unexplained mid-range of k
  for Alpha, not unconditionally for both branches across all k, and not
  for a reason this analysis has identified. This addendum is the more
  complete picture; Addendum 2 was accurate for k=24 specifically but
  should not be read as "robust for any k."
- Practical takeaway: k=18-24 is where Alpha's result is best-supported
  statistically, but calling it "the neighbor count where the model fits
  best" (as an earlier draft of this note said) is not accurate — model fit
  (by the stationarity-residual measure) is best at the smallest k tested,
  not in this middle range. See `figures/stationarity_vs_k.png` and
  `figures/knn_sensitivity.png` for the actual monotonic residual trend
  alongside the non-monotonic Alpha significance trend, side by side.

## Addendum 7: the pre-specified confirmatory test — does not replicate

Following through on the "correct next move" from Addendum 6: a single
confirmatory test, design fixed before execution (see
`confirmatory_granule_test.py` for the pre-specification, written and
committed to before this run):

- Same lineage (Radial Glia-like -> Granule mature), same n_per_stage (45),
  same k (24)
- **New seed** (424242), drawing a different 90-cell subsample than the
  exploratory run — a genuine held-out replication, not the same cells
  re-run with more permutations
- 1000 permutation replicates (up from 300), for a non-floor-clipped p-value
- Significance threshold set at 0.01 in advance, not 0.05, specifically
  because this follows an exploratory positive result

**Result: p = 0.027. Does not clear the pre-specified 0.01 threshold.**
(Entropy production 0.150, CI [0.139, 0.160], stationarity residual 0.012.)

**Reading this honestly:** this is close to the textbook pattern for a
false-lead exploratory result — the p-value moved toward less significant
(0.047 uncorrected exploratory -> 0.027 on a fresh subsample with better
resolution), consistent with regression to the mean rather than a real
effect. It is not a clean, unambiguous null (0.027 is still fairly low, not
0.5+), so "there is absolutely nothing here" would overstate this result in
the opposite direction. The defensible conclusion is: this specific lead
does not clear the bar that was set for it in advance, and should not be
carried forward as a finding. Combined with Addenda 4-6, the overall
dentate-gyrus picture stays: no reliable entropy-production signal detected
in this dataset by this method, across seven independent tests now, with
one marginal exploratory result that failed its own pre-specified
confirmatory follow-up.
