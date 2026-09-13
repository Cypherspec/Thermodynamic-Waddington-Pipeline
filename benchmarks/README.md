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
