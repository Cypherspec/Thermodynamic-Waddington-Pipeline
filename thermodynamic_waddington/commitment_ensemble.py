"""How reproducible is commitment? The committor over an ensemble of independent fits.

One fit of the landscape gives one committor. To see how much of it is data and how much is the
particular sample, refit on random subsets of the cells (80% by default) and compare. Every cell
falls in several fits, so it gets a committor mean and spread; every fit gives its own commitment
stage and its own resolution-matched barrier, so the barrier gets an uncertainty across fits rather
than only across resamples of one fit's cells (which understates it).

    from thermodynamic_waddington import commitment_ensemble
    ens = commitment_ensemble(expression, velocity, labels, ["Ductal"], ["Beta"], fits=8)
    ens.q_sd                      # per-cell spread of the committor across fits
    ens.barrier_mean_kt, ens.barrier_sd_kt
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Sequence

import numpy as np

from .config import FitConfig
from .developmental import developmental_coordinate, resolved_barrier
from .model import fit_landscape


@dataclass
class CommitmentEnsemble:
    fits: int
    fraction: float
    q_mean: list[float] = field(repr=False)     # per cell, over the fits that included it
    q_sd: list[float] = field(repr=False)       # per cell; NaN for cells in fewer than two fits
    n_fits: list[int] = field(repr=False)       # per cell, how many fits included it
    commitment_stages: list[str | None]         # per fit: first stage with mean committor > 0.5
    commitment_agreement: float                 # share of fits naming the most common stage
    barriers_kt: list[float]                    # per fit: resolution-matched, bias-corrected barrier
    barrier_se_kt: list[float]                  # per fit: bootstrap standard error over that fit's cells
    barrier_mean_kt: float
    barrier_sd_kt: float
    sd_by_phase: dict[str, float]               # median committor spread in each commitment phase
    pairwise_spearman: list[float]              # per pair of fits, on the shared cells outside source and target
    rank_agreement: float                       # mean of pairwise_spearman: do the fits order the cells alike?

    def to_dict(self) -> dict:
        out = asdict(self)
        for k in ("q_mean", "q_sd", "n_fits"):
            out.pop(k)
        return out


def _average_ranks(a: np.ndarray) -> np.ndarray:
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a))
    ranks[order] = np.arange(len(a), dtype=float)
    _, inv = np.unique(a, return_inverse=True)
    return (np.bincount(inv, weights=ranks) / np.bincount(inv))[inv]


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra, rb = _average_ranks(a), _average_ranks(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


PHASES = ((0.0, 0.2, "progenitor"), (0.2, 0.4, "early"), (0.4, 0.6, "transition state"),
          (0.6, 0.8, "late"), (0.8, 1.0001, "committed"))


def commitment_ensemble(expression: Sequence[Sequence[float]], velocity: Sequence[Sequence[float]],
                        labels: Sequence[str], source_labels: Sequence[str], target_labels: Sequence[str],
                        fits: int = 8, fraction: float = 0.8, config: FitConfig | None = None,
                        seed: int = 0) -> CommitmentEnsemble:
    """Refit the landscape on `fits` random subsets of the cells and pool the committors."""
    x = np.asarray(expression, dtype=float)
    v = np.asarray(velocity, dtype=float)
    lab = np.asarray([str(s) for s in labels])
    n = len(lab)
    if fits < 2:
        raise ValueError("an ensemble needs at least two fits")
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    m = max(3, int(round(fraction * n)))
    rng = np.random.default_rng(seed)
    total = np.zeros(n)
    total_sq = np.zeros(n)
    count = np.zeros(n, dtype=int)
    stages, barriers, ses, per_fit = [], [], [], []
    interior = ~np.isin(lab, list(source_labels) + list(target_labels))
    for k in range(fits):
        idx = np.sort(rng.choice(n, size=m, replace=False))
        cfg = config or FitConfig(neighbors=20, dimensions=6, bootstrap_replicates=2, enable_cycle_decomposition=False,
                                  entropy_production_bootstrap_replicates=2,
                                  entropy_production_permutation_replicates=2, seed=seed + k)
        fit = fit_landscape(x[idx].tolist(), v[idx].tolist(), config=cfg, labels=lab[idx].tolist(),
                            metadata={"velocity_status": "provided"})
        q = np.asarray(developmental_coordinate(fit, source_labels, target_labels))
        total[idx] += q
        full = np.full(n, np.nan)
        full[idx] = q
        per_fit.append(full)
        total_sq[idx] += q * q
        count[idx] += 1
        means = {s: float(q[lab[idx] == s].mean()) for s in set(lab[idx].tolist())}
        stages.append(next((s for s in sorted(means, key=means.get) if means[s] > 0.5), None))
        rb = resolved_barrier(q, seed=seed + k)
        barriers.append(rb.barrier_kt)
        ses.append(rb.se_kt)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(count > 0, total / np.maximum(count, 1), np.nan)
        var = np.where(count > 1, (total_sq - count * mean * mean) / np.maximum(count - 1, 1), np.nan)
    sd = np.sqrt(np.clip(var, 0.0, None))
    named = [s for s in stages if s is not None]
    top = Counter(named).most_common(1)
    by_phase = {}
    for lo, hi, name in PHASES:
        sel = (mean >= lo) & (mean < hi) & np.isfinite(sd)
        by_phase[name] = float(np.median(sd[sel])) if sel.any() else float("nan")
    b = np.asarray(barriers)
    pairs = []
    for i in range(fits):
        for j in range(i + 1, fits):
            shared = interior & np.isfinite(per_fit[i]) & np.isfinite(per_fit[j])
            if shared.sum() >= 3:
                pairs.append(_spearman(per_fit[i][shared], per_fit[j][shared]))
    return CommitmentEnsemble(
        fits=fits, fraction=fraction, q_mean=mean.tolist(), q_sd=sd.tolist(), n_fits=count.tolist(),
        commitment_stages=stages, commitment_agreement=(top[0][1] / fits) if top else 0.0,
        barriers_kt=b.tolist(), barrier_se_kt=ses, barrier_mean_kt=float(b.mean()), barrier_sd_kt=float(b.std(ddof=1)),
        sd_by_phase=by_phase, pairwise_spearman=pairs,
        rank_agreement=float(np.nanmean(pairs)) if pairs else float("nan"),
    )
