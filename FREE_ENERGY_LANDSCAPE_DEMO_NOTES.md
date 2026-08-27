# Free-energy landscape demo: what it is, and what testing it actually found

## Short answer to "can this pipeline make a free-energy landscape of a dataset of cells"

Yes — `fit_landscape()` in `thermodynamic_waddington/model.py` is real, working code,
separate from the `entropy_production.py` pipeline used throughout the rest of this
project. Given expression + RNA velocity for a set of cells, it:

1. PCA-reduces expression to a low-dimensional space
2. Builds a kNN graph and estimates local density/diffusion
3. Annotates directed graph edges with "work" (Jarzynski-style path integration,
   using velocity alignment along each edge)
4. Propagates an effective free energy outward from one reference cell, along
   directed paths, with bootstrap uncertainty
5. Detects attractor basins and barrier heights between them
6. Also runs persistent homology, spectral analysis, and several other
   diagnostics on the same fitted structure

It ran successfully on 200 real pancreatic cells (Ductal through Beta) in ~10-12
seconds and produced a full `LandscapeFit` object with real per-cell energies,
a 2D embedding, and ~90 diagnostic fields.

## What actually testing it (not just running it once) found

This feature has **not** been through the kind of repeated stress-testing the
entropy-production pipeline got earlier in this project (independent dataset
checks, sensitivity sweeps, pre-specified confirmatory tests). Running it here
and actually inspecting the output surfaced two real problems within about
fifteen minutes of checking:

### 1. Nearly half the cells get a fake placeholder value, not a real free energy

`propagate_free_energy` (in `jarzynski.py`) computes energy by propagating
outward from **one single reference cell** along **directed** edges. Any cell
with no directed path back to the reference — common in a small, sparse graph
— never gets a real value. The code silently fills it in with
`fallback = max(finite)`, i.e. the largest real value found elsewhere, dressed
up to look like a normal (if extreme) energy.

In the 200-cell demo: **97 of 200 cells (48.5%) were unreachable** and got this
fallback value. The `diagnostics['coverage']` field, which you'd expect to
catch this, reported `1.0` — it does not measure reachability and would not
have caught this on its own.

### 2. Attractor detection is highly sensitive to path-count settings, and doesn't check whether a cell's energy is real before flagging it

Two configs, otherwise identical (same 200 cells, same seed):

| Config | Attractors detected |
|---|---|
| `max_paths=32, bootstrap_replicates=8` (reduced, for demo speed) | 17 |
| `max_paths=64, bootstrap_replicates=24` (package defaults) | 94 |

That's confirmed reproducible (reran the default config three times, got 94
each time — this is parameter sensitivity, not run-to-run randomness). Of the
94 attractors detected under default settings, **77 were cells with no real
computed energy** (the fallback-filled cells from problem #1) — the attractor
detector doesn't check `path_count > 0` before flagging something as an
attractor.

## What this means

- The underlying math (PCA, kNN graph, path-integrated work, propagation) is
  real and runs correctly.
- The specific outputs — which cells are "attractors," what a given cell's
  free energy is — are **not currently trustworthy** in this small/sparse
  regime without a fix: either (a) using more reference cells / an
  undirected or bidirectional propagation, (b) explicitly marking
  unreachable cells as missing rather than fallback-filling them, and
  (c) having attractor detection check `path_count > 0` first.
- This is a legitimate, real bug to fix in the package, not a reason to
  distrust the entropy-production findings elsewhere in this project — that
  pipeline (`entropy_production.py`) is structurally different (undirected
  kNN flux estimation, no single-reference propagation) and was tested far
  more thoroughly.

## Files

- `figures/free_energy_landscape_demo.png` — the demo visualization, with
  fallback (non-real) cells shown separately in gray from cells with real
  computed energy, and both distinguished from detected attractors.
- Regenerate with: the scripts pieced together in this session (see chat
  history) using `fit_landscape` from `model.py` on
  `data/real/endocrinogenesis_day15_dynamical.h5ad`.
