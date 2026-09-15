# Getting started

This walkthrough runs the whole pipeline and explains each output. It uses
synthetic data so it needs no download; the same calls work on real data.

## 1. Make or load data

The pipeline takes an expression matrix and a velocity matrix of the same shape
(cells x features), and optional per-cell labels.

```python
from thermodynamic_waddington import make_synthetic_dataset

ds = make_synthetic_dataset(cells=200, genes=16, seed=7)
ds.expression, ds.velocity, ds.labels
```

For real data, load an h5ad and build a spliced/unspliced velocity proxy, or use
scVelo dynamical velocity. See [Concepts](concepts.md#rna-velocity).

## 2. One call

`analyze` runs everything and returns a compact report.

```python
from thermodynamic_waddington import analyze

labels = ds.labels
report = analyze(ds.expression, ds.velocity, labels=labels,
                 source_labels=[sorted(set(labels))[0]],
                 target_labels=[sorted(set(labels))[-1]])

report.is_irreversible            # from the permutation test
report.entropy_production_pvalue
report.landscape_range_kt         # landscape depth in kT
report.committor_order            # cell types ordered by commitment
report.commitment_label           # where the committor crosses 0.5
report.commitment_barrier_kt      # commitment barrier in kT
report.warnings                   # anything the analysis had to skip
```

`report.to_dict()` and `report.save("out.json")` serialize it.

## 3. Or drive the pipeline directly

For the full landscape object and all diagnostics:

```python
from thermodynamic_waddington import fit_landscape, FitConfig

cfg = FitConfig(neighbors=20, dimensions=6, seed=7)
fit = fit_landscape(ds.expression, ds.velocity, config=cfg, labels=ds.labels)

fit.energies         # per-cell effective free energy
fit.attractors       # detected attractor cells
fit.diagnostics["entropy_production"]["permutation_p_value"]
fit.diagnostics["schnakenberg"]   # cycle decomposition
```

## 4. The committor and the commitment barrier

```python
from thermodynamic_waddington import (
    developmental_coordinate, commitment_profile, committor_free_energy_profile,
)

q = developmental_coordinate(fit, source_labels=["HSC"], target_labels=["erythroid"])
profile = commitment_profile(fit, ["HSC"], ["erythroid"])
pmf = committor_free_energy_profile(q)

profile.order            # cell types by mean committor
profile.commitment_label # where commitment happens
pmf.barrier_kt           # commitment barrier
pmf.barrier_q            # committor value at the transition state
```

## 5. Calibrate to kT

```python
from thermodynamic_waddington import calibrate_fit

cal = calibrate_fit(fit)
cal.boltzmann_energy_range_kt   # landscape depth in kT
cal.r_squared                   # path-work vs density agreement
```

## Next

- [Concepts](concepts.md) explains what each number means.
- [API reference](api.md) lists every public function.
- [Benchmarks](benchmarks.md) shows how well it works, and where it does not.
