# Thermodynamic Waddington Pipeline

Estimates an effective free-energy landscape from single-cell RNA expression and velocity data, using a path-work Jarzynski-style approach on a kNN graph.

This is a research codebase, not a published package. The core claim it supports: cell differentiation trajectories show detectable entropy production that is statistically distinguishable from a null of randomly shuffled velocity vectors. That claim has been validated on pancreatic endocrinogenesis data and tested (unsuccessfully) on dentate gyrus neurogenesis -- the negative result is documented in the research note, not hidden.

## What it does

- builds a kNN graph over cells in PCA space
- annotates each edge with a path-work value combining density, velocity drift, and a diffusion noise term
- propagates effective free energies outward from a reference cell using Jarzynski exponential averaging
- estimates entropy production using Seifert's discrete-state formula, with bootstrap CIs and a label-permutation null
- fits attractor basins and computes barrier heights between them

No GPU needed. Runs on real scRNA-seq h5ad files (scVelo dynamical model output or a spliced/unspliced proxy). No external dependencies beyond numpy/scipy for optional preprocessing.

## Quick start

```python
from thermodynamic_waddington.model import fit_landscape
from thermodynamic_waddington.config import FitConfig

cfg = FitConfig(neighbors=24, dimensions=6, seed=42)
fit = fit_landscape(expression, velocity, config=cfg, labels=labels)
print(fit.energies[:5])  # per-cell effective free energy
print(fit.attractors)    # indices of detected attractor cells
```

For real scRNA-seq data, normalize first:

```python
from thermodynamic_waddington.preprocessing import velocity_in_transformed_space

pts, vel, report = velocity_in_transformed_space(expression, velocity)
fit = fit_landscape(pts, vel, config=cfg)
```

## Key modules

- `graph.py` -- kNN graph, edge work computation, local density/diffusion
- `entropy_production.py` -- Seifert EP estimator with permutation test
- `jarzynski.py` -- free-energy propagation and path ensemble
- `preprocessing.py` -- size normalization and log1p transform for real counts
- `model.py` -- top-level `fit_landscape` call, LandscapeFit output struct
- `config.py` -- all numerical settings in one place

## Running experiments from the paper

```bash
# fit scVelo dynamical model (requires scvelo, anndata)
python fit_dynamical_velocity.py

# entropy production replication (scaled up, dynamical velocity)
python run_entropy_replication_v3_dynamical.py

# kNN sensitivity sweep
python knn_sensitivity.py

# independent dataset check (dentate gyrus)
python run_independent_dataset_check.py

# pre-specified confirmatory test
python confirmatory_granule_test.py
```

## Tests

```bash
python -m unittest discover -s tests
```

180 unit tests. Should finish in under a minute.

## Data

Raw h5ad files are not included (50-200MB each). The pipeline downloads them automatically the first time if you have internet access, or you can point it at local files directly.

Pancreas data: Bastidas-Ponce et al. 2019, via the scVelo tutorial mirror.
Dentate gyrus data: Hochgerner et al. 2018, same mirror.

## Scientific status

The entropy production signal in pancreas (Beta and Alpha branches, scVelo dynamical velocity) replicated across proxy and dynamical velocity methods and held up in a kNN sensitivity sweep. It did not replicate in dentate gyrus across 7 independent tests, including one pre-registered confirmatory test. Full results are in `RESEARCH_NOTE_entropy_production_pancreas.md`.

The free-energy magnitude is not physically calibrated (no unit conversion, no temperature calibration). Use the permutation p-value for "is this nonequilibrium" claims, not the raw number.
