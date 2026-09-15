"""One-line plots for the pipeline outputs.

Needs matplotlib (the `viz` extra). Each function draws on a coordinate set you
supply (a UMAP, a PCA embedding, or `fit.embedding`), so it fits any workflow.

    from thermodynamic_waddington import plots
    plots.free_energy_profile(committor)                 # the commitment barrier
    plots.committor(committor, coords=adata.obsm["X_umap"])
    plots.landscape(fit.energies, coords=adata.obsm["X_umap"])
"""
from __future__ import annotations

import numpy as np

from .developmental import committor_free_energy_profile


def _ax(ax):
    if ax is not None:
        return ax
    import matplotlib.pyplot as plt
    _, ax = plt.subplots(figsize=(6, 5))
    return ax


def committor(values, coords, ax=None, cmap="viridis", size=14):
    """Scatter of `coords` colored by the committor (0 progenitor, 1 terminal)."""
    ax = _ax(ax)
    xy = np.asarray(coords, dtype=float)
    sc = ax.scatter(xy[:, 0], xy[:, 1], c=np.asarray(values, dtype=float),
                    cmap=cmap, s=size, edgecolors="none")
    ax.figure.colorbar(sc, ax=ax, label="committor q")
    ax.set_title("Commitment coordinate")
    ax.set_xticks([]); ax.set_yticks([])
    return ax


def landscape(energies, coords, ax=None, cmap="magma_r", size=14):
    """Scatter of `coords` colored by the effective free energy."""
    ax = _ax(ax)
    xy = np.asarray(coords, dtype=float)
    sc = ax.scatter(xy[:, 0], xy[:, 1], c=np.asarray(energies, dtype=float),
                    cmap=cmap, s=size, edgecolors="none")
    ax.figure.colorbar(sc, ax=ax, label="free energy")
    ax.set_title("Free-energy landscape")
    ax.set_xticks([]); ax.set_yticks([])
    return ax


def free_energy_profile(values, temperature=1.0, ax=None):
    """Potential of mean force along the committor, with the barrier marked."""
    ax = _ax(ax)
    prof = committor_free_energy_profile(values, temperature=temperature)
    grid = np.asarray(prof.q_grid)
    f = np.asarray(prof.free_energy_kt)
    ax.plot(grid, f, color="#0f7d99", lw=2.4)
    ax.axvline(prof.barrier_q, color="#b0413e", ls="--", lw=1)
    ax.annotate(f"barrier {prof.barrier_kt:.1f} kT",
                xy=(prof.barrier_q, prof.barrier_kt),
                xytext=(prof.barrier_q + 0.03, prof.barrier_kt * 0.85),
                color="#b0413e")
    ax.set_xlabel("committor q  (0 = progenitor, 1 = terminal fate)")
    ax.set_ylabel("free energy (kT)")
    ax.set_title("Commitment barrier")
    ax.grid(True, alpha=0.3)
    return ax
