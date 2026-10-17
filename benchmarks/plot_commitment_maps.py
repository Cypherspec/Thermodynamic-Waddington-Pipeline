"""Commitment maps from the example exports (figures/commitment_maps.png).

For each lineage: cells on their UMAP coloured by the committor, with commitment
flow (each cell toward its neighbours further along the committor) as arrows, and
the potential of mean force along the committor with the commitment phases shaded.

    python benchmarks/export_examples.py      # writes examples/exports/
    python benchmarks/plot_commitment_maps.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

EXPORTS = Path("examples/exports")
TONES = ["#006aff", "#7c4dff", "#d17a00", "#e11d6a", "#0a9a4e"]   # one tone per commitment phase
PHASES = [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]
INK, MUTED, GRID = "#171717", "#5f636e", "#dfe2ea"
LINEAGES = [("pancreas", "Pancreatic endocrinogenesis", "Ductal", "Beta"),
            ("fly_testis", "Drosophila spermatogenesis", "spermatogonia", "elongating spermatids")]


def load(prefix):
    rows = list(csv.DictReader((EXPORTS / f"{prefix}_cells.csv").open()))
    summary = json.loads((EXPORTS / f"{prefix}_summary.json").read_text())
    return rows, summary


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8, length=0)


def main(out="figures/commitment_maps.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list("commitment", TONES)
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "bold", "axes.titlesize": 11})
    available = [lin for lin in LINEAGES if (EXPORTS / f"{lin[0]}_cells.csv").exists()]
    fig, axes = plt.subplots(len(available), 2, figsize=(12, 4.6 * len(available)),
                             gridspec_kw={"width_ratios": [1.25, 1]}, squeeze=False)
    for row, (prefix, title, src, tgt) in enumerate(available):
        rows, s = load(prefix)
        xy = np.array([[float(r["tw_x"]), float(r["tw_y"])] for r in rows])
        q = np.array([float(r["tw_committor"]) for r in rows])
        index = {r["tw_id"]: i for i, r in enumerate(rows)}
        ax = axes[row, 0]
        ax.scatter(xy[:, 0], xy[:, 1], c=q, cmap=cmap, vmin=0, vmax=1, s=5, lw=0, alpha=0.9, rasterized=True)
        rng = np.random.default_rng(0)
        span = np.ptp(xy, axis=0).max()
        for i in rng.choice(len(rows), size=min(260, len(rows)), replace=False):
            ds = [index[d] for d in rows[i]["tw_downstream"].split("|") if d in index]
            if not ds:
                continue
            v = xy[ds].mean(axis=0) - xy[i]
            n = np.linalg.norm(v)
            if n == 0:
                continue
            v = v / n * span * 0.035
            ax.annotate("", xy=xy[i] + v, xytext=xy[i],
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.6, alpha=0.55, mutation_scale=6))
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_visible(False)
        ax.set_title(f"{title}: committor and commitment flow", loc="left", color=INK)
        ax.text(0.0, -0.04, f"0 = {src}  ·  1 = {tgt}  ·  {s['n_cells']:,} cells  ·  commits at "
                f"{s['commitment']['commitment_stage']}", transform=ax.transAxes, color=MUTED, fontsize=8.5, va="top")

        ax = axes[row, 1]
        style(ax)
        qg, f = np.array(s["commitment"]["pmf_q"]), np.array(s["commitment"]["pmf_kt"])
        for (lo, hi), tone in zip(PHASES, TONES):
            ax.axvspan(lo, hi, color=tone, alpha=0.07, lw=0)
        ax.fill_between(qg, f, color="#006aff", alpha=0.12, lw=0)
        ax.plot(qg, f, color="#006aff", lw=2)
        bq, bk = s["commitment"]["transition_state_q"], s["commitment"]["barrier_kt"]
        ax.axvline(bq, color="#d17a00", lw=1.2, ls=(0, (3, 3)))
        ax.scatter([bq], [bk], color="#d17a00", s=28, zorder=3)
        ax.annotate(f"{bk:.2f} kT", (bq, bk), xytext=(8, 6), textcoords="offset points", color=INK, fontsize=9, weight="bold")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, max(f) * 1.18)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_xlabel("committor q", color=MUTED, fontsize=9)
        ax.set_ylabel("free energy along q (kT)", color=MUTED, fontsize=9)
        irr = s["irreversibility"]
        ax.set_title(f"Commitment barrier · irreversibility {irr['noise_matched_excess']:+.3f} (z = {irr['z']:+.1f})",
                     loc="left", color=INK)
    names = ["progenitor", "early", "transition state", "late", "committed"]
    handles = [plt.Line2D([0], [0], marker="s", ls="", color=t, markersize=8) for t in TONES]
    fig.legend(handles, names, loc="lower center", ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
