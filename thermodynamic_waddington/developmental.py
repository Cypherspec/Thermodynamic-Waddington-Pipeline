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


def committor_free_energy_profile(committor, temperature: float = 1.0, grid: int = 60, floor: float = 1e-9) -> FreeEnergyProfile:
    """Potential of mean force along the committor: F(q) = -kT ln P(q).

    Estimates the density of cells along the committor reaction coordinate with a
    1D Gaussian KDE and Boltzmann-inverts it. The maximum over the populated
    interior is the commitment barrier in kT (the transition state), the standard
    free-energy-profile method from reaction-rate theory.
    """
    q = np.asarray(committor, dtype=float)
    q = q[np.isfinite(q)]
    if q.size < 3:
        raise ValueError("need at least 3 committor values")
    xs = np.linspace(0.0, 1.0, grid)
    bw = 1.06 * (float(np.std(q)) + 1e-9) * q.size ** (-1.0 / 5.0)
    d = (xs[:, None] - q[None, :]) / bw
    dens = np.exp(-0.5 * d * d).sum(axis=1) / (q.size * bw * np.sqrt(2.0 * np.pi))
    dens = np.clip(dens, floor, None)
    f = -temperature * np.log(dens)
    f = f - f.min()
    lo, hi = float(np.quantile(q, 0.05)), float(np.quantile(q, 0.95))
    interior = np.where((xs >= lo) & (xs <= hi))[0]
    top = interior[int(np.argmax(f[interior]))] if interior.size else int(np.argmax(f))
    return FreeEnergyProfile(xs.tolist(), f.tolist(), float(f[top]), float(xs[top]))
