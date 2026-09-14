# Reproducibility

Everything in the README, the manuscript draft, and the benchmarks reproduces
from this repository.

## Environment

Tested on Python 3.11.5. Exact package versions are pinned in
`requirements-tested.txt`. The numeric core needs only numpy; scVelo, anndata,
scanpy, h5py are used for real-data loading, and matplotlib/pandas for figures.

```bash
pip install -e ".[all]"
```

## Tests

```bash
pytest        # 193 unit tests
```

The vectorized numeric core is pinned to bit-identical outputs against a
pure-Python reference (energies max abs diff 0.0); see the runtime benchmark.

## Reproduce every result

```bash
python benchmarks/run_all.py
```

This regenerates the JSON in `experiments/` and the figures in `figures/`.
Individual benchmarks are listed in the README. Synthetic benchmarks need no
data; the real-data benchmarks need the datasets below.

## Data availability

All datasets are public and obtained through scVelo; none are redistributed here.

- Pancreatic endocrinogenesis (E15.5): Bastidas-Ponce et al., *Development* 2019,
  via `scvelo.datasets.pancreas()`. Place at
  `data/real/endocrinogenesis_day15.h5ad`.
- Mouse gastrulation erythroid: Pijuan-Sala et al., *Nature* 2019, via
  `scvelo.datasets.gastrulation_erythroid()`. Place at
  `data/real/gastrulation_erythroid.h5ad`.
- Dentate gyrus (documented negative control): Hochgerner et al., 2018, via
  `scvelo.datasets.dentategyrus()`.

## Determinism

All benchmarks take a `--seed`. Reported means and confidence intervals are over
independent random subsamples with fixed seeds, so they reproduce exactly.
