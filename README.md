# Thermodynamic Waddington Pipeline

Estimates an effective free-energy landscape from single-cell RNA expression and
velocity data using a path-work Jarzynski-style approach on a kNN graph, and
measures whether a differentiation process is thermodynamically irreversible.

Installable Python package with a pure-Python core (numpy is the only required
dependency). The core scientific claim it supports: cell differentiation
trajectories show detectable entropy production that is statistically
distinguishable from a null of randomly shuffled velocity vectors. That claim
has been validated on pancreatic endocrinogenesis data and tested
(unsuccessfully) on dentate gyrus neurogenesis -- the negative result is
documented in the research note, not hidden.

## What it does

- builds a kNN graph over cells in PCA space
- annotates each edge with a path-work value combining density, velocity drift, and a diffusion noise term
- propagates effective free energies with multi-source Jarzynski averaging so every cell gets a real value
- estimates entropy production using Seifert's discrete-state formula, with bootstrap CIs and a label-permutation null
- decomposes the entropy production into Schnakenberg cycle contributions
- computes mean first-passage times to attractor basins and velocity-divergence / gradient-alignment EP proxies
- fits attractor basins and computes barrier heights between them

No GPU needed. Runs on real scRNA-seq h5ad files (scVelo dynamical model output
or a transparent spliced/unspliced proxy).

## What makes it different

Trajectory tools such as scVelo, CellRank, Palantir, and Waddington-OT estimate
pseudotime, fate probabilities, or optimal-transport maps between stages. This
pipeline answers a different, complementary question: is the process measurably
out of equilibrium, and by how much? It quantifies entropy production,
cycle-resolved irreversibility (Schnakenberg), and gives a permutation test for
the "is this nonequilibrium" claim. Those diagnostics are not produced by the
trajectory tools above.

It is deliberate about limits. The free-energy magnitude is not physically
calibrated (no unit conversion, no temperature calibration), so the permutation
p-value, not the raw number, is what supports nonequilibrium claims.

## Install

```bash
pip install -e .                 # core, requires numpy
pip install -e ".[real-data]"    # anndata, h5py, scvelo, scanpy
pip install -e ".[viz]"          # matplotlib, plotly, pandas
pip install -e ".[dev]"          # pytest, build, twine
```

Requires Python 3.10+.

## Quick start

```python
from thermodynamic_waddington.model import fit_landscape
from thermodynamic_waddington.config import FitConfig

cfg = FitConfig(neighbors=24, dimensions=6, seed=42)
fit = fit_landscape(expression, velocity, config=cfg, labels=labels)
print(fit.energies[:5])  # per-cell effective free energy
print(fit.attractors)    # indices of detected attractor cells
print(fit.diagnostics["entropy_production"]["permutation_p_value"])
```

For real scRNA-seq data, normalize first:

```python
from thermodynamic_waddington.preprocessing import velocity_in_transformed_space

pts, vel, report = velocity_in_transformed_space(expression, velocity)
fit = fit_landscape(pts, vel, config=cfg)
```

## Performance

The numeric core was vectorized with numpy (spectral power iteration,
counterfactual scan, cached provenance and gene lookups). Fit outputs are
bit-identical to the previous pure-Python version -- energies max abs diff is
0.0, attractors and rankings match exactly -- so only the runtime changed.

| cells | before | after | speedup |
|------:|-------:|------:|--------:|
| 120   | 3.6s   | 1.6s  | 2.2x    |
| 200   | 8.3s   | 2.8s  | 2.9x    |
| 300   | 15.2s  | 3.9s  | 4.0x    |
| 400   | 28.1s  | 4.6s  | 6.2x    |

The speedup grows with cell count because the removed costs were O(n^2).

```bash
python benchmarks/runtime_scaling.py
python benchmarks/plot_scaling.py
```

![runtime speedup](figures/runtime_speedup.png)

## Real-data validation

```bash
python benchmarks/pancreas_entropy_validation.py
```

On the endocrine differentiation branch (Ductal -> Ngn3 -> Pre-endocrine ->
Beta) of the public scVelo endocrinogenesis dataset, built with a transparent
unspliced-minus-spliced velocity proxy, the observed entropy production is
significant against the label-permutation null (p ~ 0.01), while a
shuffled-velocity negative control is not (p ~ 0.1). The dynamical-velocity
result on the Beta and Alpha branches is in the research note.

## Key modules

- `graph.py` -- kNN graph, edge work computation, local density/diffusion
- `jarzynski.py`, `multi_source_propagation.py` -- free-energy propagation
- `entropy_production.py` -- Seifert EP estimator with permutation test
- `cycle_decomposition.py` -- Schnakenberg cycle decomposition of the EP
- `mfpt.py`, `divergence.py`, `gradient_alignment.py` -- attractor timescales and EP proxies
- `preprocessing.py` -- size normalization and log1p transform for real counts
- `model.py` -- top-level `fit_landscape` call, LandscapeFit output struct
- `config.py` -- all numerical settings in one place

## Running experiments from the paper

```bash
python fit_dynamical_velocity.py                  # fit scVelo dynamical model
python run_entropy_replication_v3_dynamical.py    # entropy production replication
python knn_sensitivity.py                         # kNN sensitivity sweep
python run_independent_dataset_check.py           # independent dataset check (dentate gyrus)
python confirmatory_granule_test.py               # pre-specified confirmatory test
```

## Tests

```bash
pytest
```

180 unit tests. The suite runs many full fits, so it takes a few minutes.

## Data

Raw h5ad files are not included (50-200MB each). Point the real-data commands at
local files, or download the public datasets first.

Pancreas data: Bastidas-Ponce et al. 2019, via the scVelo tutorial mirror.
Dentate gyrus data: Hochgerner et al. 2018, same mirror.

## Scientific status

The entropy production signal in pancreas (Beta and Alpha branches, scVelo
dynamical velocity) replicated across proxy and dynamical velocity methods and
held up in a kNN sensitivity sweep. It did not replicate in dentate gyrus across
7 independent tests, including one pre-registered confirmatory test. Full results
are in `RESEARCH_NOTE_entropy_production_pancreas.md`.

The free-energy magnitude is not physically calibrated (no unit conversion, no
temperature calibration). Use the permutation p-value for "is this
nonequilibrium" claims, not the raw number.
