# Benchmarks

Reproducible performance benchmarks for the pipeline.

## Runtime scaling

`runtime_scaling.py` times the full `fit_landscape` call across a range of cell
counts and compares against a recorded pure-Python baseline (commit `40d1fcd`,
before the numpy vectorization of the spectral, causal, and provenance paths).

```bash
python benchmarks/runtime_scaling.py --sizes 120 200 300 400
python benchmarks/plot_scaling.py
```

Outputs:

- `experiments/runtime_scaling.json` - timing table and speedups
- `figures/runtime_speedup.png` - runtime and speedup plots

The speedup grows with cell count (2.2x at 120 cells to 6.2x at 400) because the
vectorized paths removed the dominant O(n^2) pure-Python cost. The fit outputs
are bit-identical between the two versions: energies max abs diff is 0.0, and
attractors and the counterfactual ranking match exactly, so this is pure speed,
not a quality tradeoff.

## Drosophila spermatogenesis

`fly_testis_prepare.py` builds spliced/unspliced counts for the Fly Cell Atlas
testis from its raw reads (ArrayExpress E-MTAB-10519) and attaches the atlas
annotation; `fly_testis.py` runs the calibrated readouts on the germline.

```bash
python benchmarks/fly_testis_prepare.py --mb-per-lane 1200   # needs kb-python
python benchmarks/fly_testis.py
```

Outputs: `experiments/fly_testis.json`, `figures/fly_testis.png`. Results and
caveats are in `docs/benchmarks.md`.

## Resolution-matched commitment barriers

`barrier_resolution.py` asks whether a commitment barrier is comparable across
datasets of different sizes. It compares Silverman's bandwidth with a fixed
resolution (h = 0.08, reflected at q = 0 and 1, bootstrap bias-corrected) on a
double well with an exact committor, on independent refits of pancreas and fly
spermatogenesis at three cell numbers, and across bandwidths.

```bash
python benchmarks/barrier_resolution.py --seeds 5
```

Outputs: `experiments/barrier_resolution.json`, `figures/barrier_resolution.png`.

## Commitment ensemble

`commitment_ensemble.py` refits the landscape on ten random 80% subsets of 1,500
cells per lineage and pools the committors: where the committor is uncertain,
whether every fit names the same commitment stage, and how the spread of the
barrier between fits compares with one fit's bootstrap error.

```bash
python benchmarks/commitment_ensemble.py --fits 10 --cells 1500
```

Outputs: `experiments/commitment_ensemble.json`, `figures/commitment_ensemble.png`.
