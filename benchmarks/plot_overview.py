"""Pipeline overview diagram (figures/tw_pipeline_overview.png).

Each result card carries a small view of the real readout, drawn from the pancreas example export
when it is present: the cyclic fraction against its noise-matched null, the potential of mean force
along the committor, and the density landscape over the embedding coloured by the committor.

    python benchmarks/plot_overview.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

INK, MUTED, LINE, CARD = "#1f2328", "#5f6672", "#9aa1ab", "#ffffff"
TEAL, PURPLE, AMBER, RED = "#0f7d99", "#6a4c93", "#b7791f", "#c0392b"
EXPORT = Path("examples/exports/pancreas")


def load_export():
    cells, summary = Path(f"{EXPORT}_cells.csv"), Path(f"{EXPORT}_summary.json")
    if not (cells.exists() and summary.exists()):
        return None, None, None
    rows = list(csv.DictReader(cells.open()))
    q = np.array([float(r["tw_committor"]) for r in rows])
    xy = np.array([[float(r["tw_x"]), float(r["tw_y"])] for r in rows])
    return q, xy, json.loads(summary.read_text())


def main(out="figures/tw_pipeline_overview.png"):
    family = "Inter" if any("Inter" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
    plt.rcParams.update({"font.family": family})
    q, xy, summary = load_export()

    fig = plt.figure(figsize=(14, 6.6), dpi=200, facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 6.6)
    ax.axis("off")

    def card(x, y, w, h, edge=LINE, lw=1.0, fill=CARD):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.1", fc=fill, ec=edge,
                                    lw=lw, zorder=2))

    def arrow(x0, y0, x1, y1, color=LINE):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=3,
                    arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.22", color=color, lw=1.2,
                                    shrinkA=2, shrinkB=2))

    ax.text(0.35, 6.35, "Thermodynamic Waddington Pipeline", fontsize=16, weight="bold", color=INK, va="top")
    ax.text(0.35, 5.88, "A calibrated thermodynamic layer on top of any RNA-velocity workflow: where a lineage "
            "commits, how high the barrier is, and whether the flow is irreversible beyond its own noise.",
            fontsize=9.5, color=MUTED, va="top")

    def step(x, y, w, h, tag, title, body):
        card(x, y, w, h, fill="#f7f8fa", edge="#c4c9d1")
        ax.text(x + 0.18, y + h - 0.18, tag, fontsize=8, color=MUTED, weight="bold", va="top", zorder=3)
        ax.text(x + 0.18, y + h - 0.48, title, fontsize=11.5, weight="bold", color=INK, va="top", zorder=3)
        ax.text(x + 0.18, y + h - 0.9, body, fontsize=8.8, color=INK, va="top", linespacing=1.5, zorder=3)

    step(0.35, 1.9, 2.45, 2.75, "1  INPUT", "Cells + velocity",
         "scRNA-seq expression\n+ RNA velocity\n(scVelo, or a spliced/\nunspliced proxy)")
    step(3.25, 1.9, 2.45, 2.75, "2  GRAPH", "Cell-state graph",
         "PCA embedding\nexact KD-tree kNN\n(50k cells in seconds)\n\nvelocity = flow\non every edge")
    arrow(2.8, 3.27, 3.25, 3.27)

    rows = [
        (3.75, TEAL, "3  HODGE", "Hodge decomposition", "flow = gradient (reversible)\n+ cyclic (irreversible)",
         "Calibrated irreversibility", "noise-matched null:\nfalse positives 100% → 0%"),
        (2.15, PURPLE, "4  TPT", "Transition-path theory", "committor q: 0 = progenitor,\n1 = terminal fate",
         "Commitment barrier", "F(q) = −kT ln P(q)\npancreas 1.6 kT [1.4, 1.9]"),
        (0.55, AMBER, "5  JARZYNSKI", "Path work", "free-energy landscape with a\ndensity-anchored kT scale",
         "kT landscape", "pseudopotential in kT\n≈ 5 kT deep (pancreas)"),
    ]
    for y, color, tag, t1, b1, t2, b2 in rows:
        card(6.2, y, 3.0, 1.4, edge="#c4c9d1")
        ax.text(6.38, y + 1.24, tag, fontsize=7.5, color=MUTED, weight="bold", va="top", zorder=3)
        ax.text(6.38, y + 0.98, t1, fontsize=11, weight="bold", color=color, va="top", zorder=3)
        ax.text(6.38, y + 0.62, b1, fontsize=8.5, color=INK, va="top", linespacing=1.45, zorder=3)
        card(9.75, y, 3.9, 1.4, edge=color, lw=1.6)
        ax.text(9.93, y + 1.22, t2, fontsize=11, weight="bold", color=color, va="top", zorder=3)
        ax.text(9.93, y + 0.84, b2, fontsize=8.5, color=INK, va="top", linespacing=1.45, zorder=3)
        arrow(5.7, 3.27, 6.2, y + 0.7)
        arrow(9.2, y + 0.7, 9.75, y + 0.7, color)

    def inset(y):
        a = fig.add_axes([12.0 / 14, (y + 0.17) / 6.6, 1.5 / 14, 1.06 / 6.6])
        a.axis("off")
        return a

    a = inset(3.75)
    vals = [0.100, 0.085] if summary is None else [summary["irreversibility"]["noise_matched_null"],
                                                  summary["irreversibility"]["cyclic_fraction"]]
    a.barh([0, 1], vals, color=["#c4c9d1", TEAL], height=0.55)
    a.set_xlim(0, max(vals) * 1.45)
    for i, (name, v) in enumerate(zip(["noise-matched null", "observed"], vals)):
        a.text(v, i, f" {v:.3f}", va="center", fontsize=7, color=INK)
        a.text(0, i + 0.38, name, fontsize=6.5, color=MUTED)

    a = inset(2.15)
    if q is not None:
        xs = np.linspace(0, 1, 60)
        bw = 1.06 * q.std() * len(q) ** -0.2
        dens = np.exp(-0.5 * ((xs[:, None] - q[None, :]) / bw) ** 2).sum(1)
        f = -np.log(dens / dens.max())
        a.plot(xs, f, color=PURPLE, lw=1.6)
        a.fill_between(xs, f, 0, color=PURPLE, alpha=0.08)
        inner = (xs > np.quantile(q, 0.05)) & (xs < np.quantile(q, 0.95))
        top = int(np.argmax(np.where(inner, f, -1)))
        a.axvline(xs[top], color=RED, lw=0.8, ls="--")
        a.text(0.5, -0.02, "committor q", transform=a.transAxes, ha="center", va="top", fontsize=6.5, color=MUTED)
        a.set_ylim(0, f.max() * 1.1)

    a = inset(0.55)
    if xy is not None:
        n = 90
        gx, gy = np.meshgrid(*(np.linspace(xy[:, k].min(), xy[:, k].max(), n) for k in range(2)))
        h = 0.05 * float(np.ptp(xy, axis=0).max())
        d2 = (gx.ravel()[:, None] - xy[None, :, 0]) ** 2 + (gy.ravel()[:, None] - xy[None, :, 1]) ** 2
        w = np.exp(-0.5 * d2 / h / h)
        dens = w.sum(1).reshape(gx.shape)
        qf = (w @ q / np.maximum(w.sum(1), 1e-300)).reshape(gx.shape)
        u = -np.log(dens / dens.max() + 0.04)
        rgba = plt.get_cmap("viridis")(qf)
        rgba[..., 3] = np.clip((u.max() - u) / (u.max() - u.min()), 0, 1) ** 0.7
        a.imshow(rgba, origin="lower", aspect="auto", extent=(0, 1, 0, 1))
        a.contour(np.linspace(0, 1, n), np.linspace(0, 1, n), u, levels=6, colors="#5f6672", linewidths=0.4)

    ax.text(0.35, 1.38, "Validated on 4 lineages in 2 phyla", fontsize=9.5, color=INK, weight="bold", va="top")
    ax.text(0.35, 1.06, "pancreas, gastrulation, bone marrow,\nDrosophila testis (Fly Cell Atlas)", fontsize=9,
            color=MUTED, va="top", linespacing=1.45)
    ax.text(0.35, 0.36, "233 tests · lint-clean CI · typed · MIT", fontsize=8.5, color=MUTED, va="top")
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
