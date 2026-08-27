"""Entropy production estimator for single-cell RNA velocity data.

Uses Seifert's discrete-state formula:
    sigma = 0.5 * sum_pairs (p_i*k_ij - p_j*k_ji) * log(p_i*k_ij / p_j*k_ji)
where k_ij come from the same work functional used for the landscape fit.

Not a physically calibrated rate (no unit conversion). Use the permutation
p-value for the "is this nonequilibrium" claim, not the raw magnitude.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Sequence

from .arrays import mean, quantile, safe_log
from .graph import Edge, NeighborGraph, edge_work
from .observability import fingerprint


@dataclass(frozen=True)
class EntropyProductionConfig:
    temperature: float = 1.0
    velocity_scale: float = 1.0
    rate_floor: float = 1e-9
    bootstrap_replicates: int = 64
    permutation_replicates: int = 64
    seed: int = 17
    work_scale: float | None = None  # None = auto from data
    max_scaled_work: float = 30.0

    def normalized(self) -> "EntropyProductionConfig":
        return EntropyProductionConfig(
            temperature=max(1e-9, float(self.temperature)),
            velocity_scale=float(self.velocity_scale),
            rate_floor=max(1e-15, float(self.rate_floor)),
            bootstrap_replicates=max(0, int(self.bootstrap_replicates)),
            permutation_replicates=max(0, int(self.permutation_replicates)),
            seed=int(self.seed),
            work_scale=None if self.work_scale is None else max(1e-9, float(self.work_scale)),
            max_scaled_work=max(1.0, float(self.max_scaled_work)),
        )


@dataclass(frozen=True)
class PairFlux:
    i: int
    j: int
    rate_ij: float
    rate_ji: float
    work_ij: float
    work_ji: float
    contribution: float  # (p_i*k_ij - p_j*k_ji) * log(...)
    net_flux: float      # signed: positive = net i->j current


def _candidate_pairs(graph: NeighborGraph) -> list[tuple[int, int]]:
    # all unordered pairs that appear in either direction of the kNN graph
    seen: set[tuple[int, int]] = set()
    pairs: list[tuple[int, int]] = []
    for i, nbrs in enumerate(graph.neighbors):
        for j in nbrs:
            key = (i, j) if i < j else (j, i)
            if key not in seen:
                seen.add(key)
                pairs.append(key)
    return pairs


def _robust_scale(vals: Sequence[float], floor: float = 1.0) -> float:
    # MAD * 1.4826, consistent estimator of std under normality (Huber 1981)
    # much less sensitive to outlier work values than plain std
    if not vals:
        return floor
    s = sorted(vals)
    n = len(s)
    med = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    devs = sorted(abs(v - med) for v in vals)
    mad = devs[n // 2] if n % 2 else 0.5 * (devs[n // 2 - 1] + devs[n // 2])
    return max(floor, 1.4826 * mad)


def _rate(w: float, floor: float, scale: float, cap: float) -> float:
    # Arrhenius form, work normalized by robust scale to avoid exp overflow
    scaled = w / max(scale, 1e-9)
    return max(floor, math.exp(-max(-cap, min(cap, scaled))))


def compute_pair_fluxes(
    pts: Sequence[Sequence[float]],
    vels: Sequence[Sequence[float]],
    dens: Sequence[float],
    diff: Sequence[float],
    graph: NeighborGraph,
    config: EntropyProductionConfig | None = None,
) -> tuple[list[PairFlux], float]:
    """Compute rates and EP contribution for every candidate pair.

    Returns (fluxes, work_scale) -- scale is exposed so callers can see
    what normalization was applied.
    """
    cfg = (config or EntropyProductionConfig()).normalized()
    total_d = sum(dens) or 1.0
    p = [max(1e-15, d / total_d) for d in dens]

    pairs = _candidate_pairs(graph)
    raw: list[tuple[int, int, object, object]] = []
    for i, j in pairs:
        fwd = edge_work(Edge(i, j, 0.0, 0.0), pts, vels, dens, diff, cfg.temperature, cfg.velocity_scale)
        bwd = edge_work(Edge(j, i, 0.0, 0.0), pts, vels, dens, diff, cfg.temperature, cfg.velocity_scale)
        raw.append((i, j, fwd, bwd))

    scale = cfg.work_scale if cfg.work_scale is not None else \
        _robust_scale([r[2].work for r in raw] + [r[3].work for r in raw])

    fluxes: list[PairFlux] = []
    for i, j, fwd, bwd in raw:
        k_ij = _rate(fwd.work, cfg.rate_floor, scale, cfg.max_scaled_work)
        k_ji = _rate(bwd.work, cfg.rate_floor, scale, cfg.max_scaled_work)
        a = p[i] * k_ij
        b = p[j] * k_ji
        contrib = (a - b) * safe_log(a / max(b, 1e-300)) if b > 0 else 0.0
        fluxes.append(PairFlux(i, j, k_ij, k_ji, fwd.work, bwd.work, contrib, a - b))
    return fluxes, scale


@dataclass
class EntropyProductionReport:
    entropy_production_rate: float
    n_pairs: int
    n_cells: int
    fluxes: list[PairFlux]
    stationarity_residual: float
    bootstrap_ci: tuple[float, float] | None
    permutation_p_value: float | None
    null_mean: float | None
    config: EntropyProductionConfig
    audit: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "entropy_production_rate": self.entropy_production_rate,
            "n_pairs": self.n_pairs,
            "n_cells": self.n_cells,
            "stationarity_residual": self.stationarity_residual,
            "bootstrap_ci": list(self.bootstrap_ci) if self.bootstrap_ci else None,
            "permutation_p_value": self.permutation_p_value,
            "null_mean": self.null_mean,
            "audit": self.audit,
            "claim_boundary": (
                "Entropy production rate is computed from the fitted work functional, "
                "not independently measured transition rates. Finite-sample EP estimators "
                "are upward-biased -- use the permutation p-value for directional claims, "
                "not the raw magnitude."
            ),
        }


def stationarity_residual(fluxes: Sequence[PairFlux], n: int) -> float:
    """RMS per-cell net probability flux. Zero when density proxy is truly
    stationary under the rate matrix. Large values mean density and dynamics
    are inconsistent with a steady state."""
    if n == 0:
        return 0.0
    net = [0.0] * n
    for f in fluxes:
        net[f.i] += f.net_flux
        net[f.j] -= f.net_flux
    return math.sqrt(mean(v * v for v in net))


def _total_from_fluxes(fluxes: Sequence[PairFlux]) -> float:
    # 0.5 matches Seifert (2012) eq.27 convention; each pair appears once here
    return 0.5 * sum(f.contribution for f in fluxes) if fluxes else 0.0


def estimate_entropy_production(
    pts: Sequence[Sequence[float]],
    vels: Sequence[Sequence[float]],
    dens: Sequence[float],
    diff: Sequence[float],
    graph: NeighborGraph,
    config: EntropyProductionConfig | None = None,
) -> EntropyProductionReport:
    """Run full pipeline: pair fluxes -> EP rate -> bootstrap CI -> permutation test.

    Bootstrap resamples WITH replacement (same size as original) -- resampling
    without replacement biases a sum-type statistic downward.
    Permutation null shuffles velocity vectors across cells, destroying any
    real position-velocity coupling while keeping marginal distributions intact.
    """
    cfg = (config or EntropyProductionConfig()).normalized()
    fluxes, ws = compute_pair_fluxes(pts, vels, dens, diff, graph, cfg)
    total = _total_from_fluxes(fluxes)
    resid = stationarity_residual(fluxes, len(pts))

    rng = random.Random(cfg.seed)

    ci = None
    if cfg.bootstrap_replicates > 1 and fluxes:
        n = len(fluxes)
        samples = [_total_from_fluxes([fluxes[rng.randrange(n)] for _ in range(n)])
                   for _ in range(cfg.bootstrap_replicates)]
        ci = (quantile(samples, 0.025), quantile(samples, 0.975))

    p_val = None
    null_mu = None
    if cfg.permutation_replicates > 0 and len(pts) > 1:
        idx = list(range(len(pts)))
        nulls = []
        for _ in range(cfg.permutation_replicates):
            sh = idx[:]
            rng.shuffle(sh)
            pv = [vels[sh[i]] for i in range(len(pts))]
            nf, _ = compute_pair_fluxes(pts, pv, dens, diff, graph, cfg)
            nulls.append(_total_from_fluxes(nf))
        null_mu = mean(nulls)
        exceed = sum(1 for v in nulls if v >= total)
        p_val = (exceed + 1) / (len(nulls) + 1)  # add-one so p never equals 0

    audit = {
        "method": "seifert_markov_jump_entropy_production",
        "n_pairs": len(fluxes),
        "n_cells": len(pts),
        "work_scale": ws,
        "work_scale_source": "user_specified" if cfg.work_scale is not None else "auto_mad",
        "fingerprint": fingerprint({"n": len(pts), "pairs": len(fluxes), "seed": cfg.seed}),
    }

    return EntropyProductionReport(
        entropy_production_rate=total,
        n_pairs=len(fluxes),
        n_cells=len(pts),
        fluxes=fluxes,
        stationarity_residual=resid,
        bootstrap_ci=ci,
        permutation_p_value=p_val,
        null_mean=null_mu,
        config=cfg,
        audit=audit,
    )


def cross_validate_with_jarzynski(report: EntropyProductionReport, jarzynski_records: Sequence[object]) -> dict[str, object]:
    """Compare graph-based EP estimate against per-node Jarzynski dissipation.
    The two use different representations of the same dynamics; large disagreement
    is itself diagnostic information, not necessarily an error."""
    dissipations = []
    for rec in (jarzynski_records or []):
        v = rec.get("dissipation") if isinstance(rec, dict) else getattr(rec, "dissipation", None)
        if v is not None:
            dissipations.append(float(v))
    jmu = mean(dissipations) if dissipations else None
    agree = None
    if jmu is not None and report.entropy_production_rate >= 0:
        denom = max(abs(jmu), report.entropy_production_rate, 1e-9)
        agree = 1.0 - abs(jmu - report.entropy_production_rate) / denom
    return {
        "graph_based_entropy_production": report.entropy_production_rate,
        "jarzynski_mean_dissipation": jmu,
        "relative_agreement": agree,
    }
