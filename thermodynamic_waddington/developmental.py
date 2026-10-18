"""Developmental coordinate and commitment profile from the committor.

The forward committor q(x) is the probability that a cell at state x reaches the
terminal fate before returning to the progenitor state. It is the principled
reaction coordinate for a non-equilibrium process: 0 at the progenitor, 1 at the
terminal fate, monotone along commitment. Because it uses the velocity-derived
directed flow and the free-energy landscape, not just expression geometry, it
orders cells along differentiation better than a linear pseudotime axis, and the
place where it crosses 0.5 is where fate commitment happens.

    from thermodynamic_waddington import developmental_coordinate, commitment_profile
    q = developmental_coordinate(fit, source_labels=["Ductal"], target_labels=["Beta"])
    report = commitment_profile(fit, ["Ductal"], ["Beta"])
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .expansion import committor_iteration, reactive_flux
from .graph import Edge


@dataclass(frozen=True)
class FreeEnergyProfile:
    q_grid: list[float]          # committor values along the reaction coordinate
    free_energy_kt: list[float]  # potential of mean force at each, in kT
    barrier_kt: float            # height of the commitment barrier
    barrier_q: float             # committor value at the barrier top (transition state)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class CommitmentReport:
    order: list[str]                 # cell-type labels sorted by mean committor
    mean_committor: dict[str, float] # per label
    commitment_label: str | None     # first label whose mean committor exceeds 0.5
    n_reactive_edges: int
    peak_flux_label: str | None      # label at the reactive-flux bottleneck (transition state)

    def to_dict(self) -> dict:
        return asdict(self)


def _edges_from_fit(fit) -> list[Edge]:
    out = []
    for e in fit.edges:
        if isinstance(e, dict):
            out.append(Edge(e["source"], e["target"], e["distance"], e["alignment"], e["work"], e.get("action", "forward")))
        else:
            out.append(e)
    return out


def _temperature(fit) -> float:
    cfg = getattr(fit, "config", {})
    if isinstance(cfg, dict):
        return float(cfg.get("temperature", 1.0))
    return float(getattr(cfg, "temperature", 1.0))


def developmental_coordinate(fit, source_labels, target_labels, iterations: int = 400) -> list[float]:
    """Per-cell committor: 0 at the progenitor states, 1 at the terminal fates."""
    labels = list(fit.labels)
    src = {i for i, l in enumerate(labels) if l in set(source_labels)}
    tgt = {i for i, l in enumerate(labels) if l in set(target_labels)}
    if not src:
        raise ValueError("no cells matched source_labels")
    if not tgt:
        raise ValueError("no cells matched target_labels")
    edges = _edges_from_fit(fit)
    return committor_iteration(edges, fit.energies, src, tgt, temperature=_temperature(fit), iterations=iterations)


def commitment_profile(fit, source_labels, target_labels, iterations: int = 400) -> CommitmentReport:
    """Order cell types by committor and locate where commitment happens."""
    q = np.asarray(developmental_coordinate(fit, source_labels, target_labels, iterations), dtype=float)
    labels = np.asarray(list(fit.labels))
    means = {}
    for lab in set(labels.tolist()):
        vals = q[labels == lab]
        if vals.size:
            means[str(lab)] = float(vals.mean())
    order = sorted(means, key=lambda k: means[k])
    commitment_label = next((lab for lab in order if means[lab] > 0.5), None)

    edges = _edges_from_fit(fit)
    flux = reactive_flux(edges, q.tolist(), fit.energies, temperature=_temperature(fit))
    reactive = [f for f in flux if f["flux"] > 0.0]
    peak_label = None
    if reactive:
        top = max(reactive, key=lambda f: f["flux"])
        peak_label = str(labels[top["source"]])
    return CommitmentReport(
        order=order,
        mean_committor=means,
        commitment_label=commitment_label,
        n_reactive_edges=len(reactive),
        peak_flux_label=peak_label,
    )


#: Resolution of the fixed-bandwidth potential of mean force, in committor units. Chosen so that the
#: barrier is stable from a few hundred cells upward (benchmarks/barrier_resolution.py).
RESOLVED_BANDWIDTH = 0.08


def committor_free_energy_profile(committor, temperature: float = 1.0, grid: int = 60, floor: float = 1e-9,
                                  bandwidth: float | None = None) -> FreeEnergyProfile:
    """Potential of mean force along the committor: F(q) = -kT ln P(q).

    Estimates the density of cells along the committor reaction coordinate with a
    1D Gaussian KDE and Boltzmann-inverts it. The maximum over the populated
    interior is the commitment barrier in kT (the transition state), the standard
    free-energy-profile method from reaction-rate theory.

    ``bandwidth=None`` (the default) uses Silverman's rule, whose bandwidth shrinks as
    n^(-1/5): the profile sharpens and the barrier rises with the number of cells, so
    barriers are comparable only at matched cell numbers. A float fixes the resolution
    in committor units (``RESOLVED_BANDWIDTH`` is the recommended value) and reflects the
    kernel at q = 0 and q = 1, where cells pile up; the barrier then estimates the same
    quantity whatever the number of cells, and differently sized datasets can be compared.
    """
    q = np.asarray(committor, dtype=float)
    q = q[np.isfinite(q)]
    if q.size < 3:
        raise ValueError("need at least 3 committor values")
    xs = np.linspace(0.0, 1.0, grid)
    norm = q.size * np.sqrt(2.0 * np.pi)
    if bandwidth is None:
        bw = 1.06 * (float(np.std(q)) + 1e-9) * q.size ** (-1.0 / 5.0)
        d = (xs[:, None] - q[None, :]) / bw
        dens = np.exp(-0.5 * d * d).sum(axis=1) / (norm * bw)
    else:
        bw = float(bandwidth)
        if not bw > 0:
            raise ValueError("bandwidth must be positive")
        dens = np.zeros_like(xs)
        for mirror in (q, -q, 2.0 - q):
            d = (xs[:, None] - mirror[None, :]) / bw
            dens += np.exp(-0.5 * d * d).sum(axis=1)
        dens /= norm * bw
    dens = np.clip(dens, floor, None)
    f = -temperature * np.log(dens)
    f = f - f.min()
    lo, hi = float(np.quantile(q, 0.05)), float(np.quantile(q, 0.95))
    interior = np.where((xs >= lo) & (xs <= hi))[0]
    top = interior[int(np.argmax(f[interior]))] if interior.size else int(np.argmax(f))
    return FreeEnergyProfile(xs.tolist(), f.tolist(), float(f[top]), float(xs[top]))


@dataclass(frozen=True)
class ResolvedBarrier:
    barrier_kt: float        # fixed-resolution barrier, bootstrap bias-corrected
    raw_kt: float            # fixed-resolution barrier before the correction
    bias_kt: float           # estimated small-sample bias (raw minus corrected)
    se_kt: float             # bootstrap standard error over resampled cells
    barrier_q: float         # transition state of the raw profile
    bandwidth: float
    n_cells: int

    def to_dict(self) -> dict:
        return asdict(self)


def resolved_barrier(committor, bandwidth: float = RESOLVED_BANDWIDTH, bootstraps: int = 60,
                     temperature: float = 1.0, seed: int = 0) -> ResolvedBarrier:
    """The commitment barrier at a fixed resolution, corrected for small-sample bias.

    The barrier is the maximum of -ln of an estimated density over the transition region, and
    sampling noise in a sparse region pushes that maximum upward. The bootstrap estimates the
    push: resample the cells with replacement, recompute the barrier, and subtract the mean excess
    (2 * raw - mean(bootstrap)). On a double well with an exact answer this takes the estimate from
    +0.24 kT to +0.01 kT at 400 cells (benchmarks/barrier_resolution.py). The standard error
    covers resampling of these cells only; independent fits of the landscape vary more.
    """
    q = np.asarray(committor, dtype=float)
    q = q[np.isfinite(q)]
    raw = committor_free_energy_profile(q, temperature=temperature, bandwidth=bandwidth)
    rng = np.random.default_rng(seed)
    boot = np.asarray([committor_free_energy_profile(q[rng.integers(0, q.size, q.size)], temperature=temperature,
                                                     bandwidth=bandwidth).barrier_kt for _ in range(bootstraps)])
    corrected = 2.0 * raw.barrier_kt - float(boot.mean())
    return ResolvedBarrier(float(corrected), float(raw.barrier_kt), float(raw.barrier_kt - corrected),
                           float(boot.std(ddof=1)), float(raw.barrier_q), float(bandwidth), int(q.size))


@dataclass(frozen=True)
class BarrierConvergence:
    sizes: list[int]             # numbers of cells the barrier was estimated from
    mean_kt: list[float]         # mean barrier over subsamples at each size
    low_kt: list[float]          # 2.5th percentile
    high_kt: list[float]         # 97.5th percentile
    bandwidth: float | None      # None: Silverman's rule; float: fixed resolution
    drift_kt: float              # barrier at the largest size minus barrier at the smallest

    def to_dict(self) -> dict:
        return asdict(self)


def barrier_convergence(committor, sizes=(200, 400, 800, 1600), bandwidth: float | None = RESOLVED_BANDWIDTH,
                        replicates: int = 40, temperature: float = 1.0, seed: int = 0) -> BarrierConvergence:
    """How the commitment barrier depends on the number of cells it is estimated from.

    Subsamples the committor values without replacement at each size (sizes larger than
    the data are dropped; the full data is always included) and reports the barrier's
    mean and 95% range. A resolution-consistent estimator has a small ``drift_kt``; the
    default Silverman estimator drifts upward. This measures the estimator, not the
    variability across independent fits of the landscape.
    """
    q = np.asarray(committor, dtype=float)
    q = q[np.isfinite(q)]
    rng = np.random.default_rng(seed)
    use = sorted({int(n) for n in sizes if 3 <= int(n) < q.size} | {int(q.size)})
    means, lows, highs = [], [], []
    for n in use:
        reps = replicates if n < q.size else 1
        b = [committor_free_energy_profile(q if n == q.size else q[rng.choice(q.size, n, replace=False)],
                                           temperature=temperature, bandwidth=bandwidth).barrier_kt
             for _ in range(reps)]
        means.append(float(np.mean(b)))
        lows.append(float(np.percentile(b, 2.5)))
        highs.append(float(np.percentile(b, 97.5)))
    return BarrierConvergence(use, means, lows, highs, bandwidth, means[-1] - means[0])
