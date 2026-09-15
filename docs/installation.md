# Installation

Requires Python 3.10 or newer. The numeric core needs only numpy; everything else
is optional.

## From source

```bash
git clone https://github.com/Cypherspec/Thermodynamic-Waddington-Pipeline
cd Thermodynamic-Waddington-Pipeline
pip install -e .
```

## Optional extras

```bash
pip install -e ".[real-data]"   # anndata, h5py, scvelo, scanpy - load real h5ad
pip install -e ".[viz]"         # matplotlib, plotly, pandas - figures
pip install -e ".[fast]"        # scipy, scikit-learn - KD-tree graph acceleration
pip install -e ".[dev]"         # pytest, build, twine, ruff
pip install -e ".[all]"         # everything, the tested environment
```

## Verify

```bash
python -c "import thermodynamic_waddington as t; print(t.__version__)"
pytest        # 197 unit tests
python examples/tutorial.py   # end-to-end on synthetic data, no download
```

## Reproducible environment

Exact tested versions are pinned in
[`requirements-tested.txt`](https://github.com/Cypherspec/Thermodynamic-Waddington-Pipeline/blob/main/requirements-tested.txt)
(Python 3.11.5, numpy 2.4.6, scipy 1.17.1, scvelo 0.3.4, and so on). See
[REPRODUCIBILITY.md](https://github.com/Cypherspec/Thermodynamic-Waddington-Pipeline/blob/main/REPRODUCIBILITY.md).
