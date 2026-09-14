# Contributing

Thanks for your interest. This is a research package; contributions that improve
correctness, tests, documentation, or reproducibility are very welcome.

## Setup

```bash
pip install -e ".[dev]"
```

## Before opening a pull request

```bash
pytest            # all tests must pass (193+)
ruff check .      # lint
ruff format .     # format
```

## Guidelines

- Keep the numeric core output-preserving: changes that alter fitted energies,
  attractors, or entropy-production values must be justified and, where possible,
  verified against the previous version (the vectorization commits show the
  pattern).
- New scientific claims need a benchmark under `benchmarks/` and a figure, with
  the result reported straight, negatives included.
- Add or update unit tests for any new public function.
- Public API changes should be reflected in `README.md` and `CHANGELOG.md`.

## Reporting results

Benchmarks write JSON to `experiments/` and figures to `figures/`. Reproduce
everything with `python benchmarks/run_all.py`.
