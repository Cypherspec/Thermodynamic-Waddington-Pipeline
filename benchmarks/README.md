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
