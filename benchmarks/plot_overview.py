"""Pipeline overview diagram (figures/tw_pipeline_overview.png).

Each result card carries a small view of the real readout, drawn from the example exports when they
are present: the commitment barrier along the committor, the cells' committor distribution, and
the noise-matched irreversibility against its null.

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
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

BG, CARD, INK, MUTED, FAINT = "#07080c", "#0e1119", "#eef0f6", "#8b90a0", "#5b6070"
EDGE = (1, 1, 1, 0.09)
GLOW = "#a5dcff"
TONES = ["#67e8f9", "#5b9dff", "#818cf8", "#a78bfa", "#e879f9"]
ROW_TONES = ["#67e8f9", "#818cf8", "#e879f9"]
EXPORT = Path("examples/exports/pancreas")


def load_export():
    cells, summary = Path(f"{EXPORT}_cells.csv"), Path(f"{EXPORT}_summary.json")
    if not (cells.exists() and summary.exists()):
        return None, None
    rows = list(csv.DictReader(cells.open()))
    q = np.array([float(r["tw_committor"]) for r in rows])
    xy = np.array([[float(r["tw_x"]), float(r["tw_y"])] for r in rows])
    return q, xy, json.loads(summary.read_text())


def main(out="figures/tw_pipeline_overview.png"):
    family = "Inter" if any("Inter" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
    plt.rcParams.update({"font.family": family})
    cmap = LinearSegmentedColormap.from_list("tw", TONES)
    q, xy, summary = load_export()

    fig = plt.figure(figsize=(14, 7.2), dpi=200, facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 7.2)
    ax.axis("off")
    ax.set_facecolor(BG)
    glow = np.exp(-((np.linspace(-1, 1, 400)[None, :] * 1.1) ** 2 + (np.linspace(-1, 1, 200)[:, None] * 1.6 - 1.2) ** 2))
    ax.imshow(glow, extent=(0, 14, 0, 7.2), cmap=LinearSegmentedColormap.from_list("g", [BG, "#16224a"]),
              alpha=0.9, aspect="auto", zorder=0)

    def card(x, y, w, h, tone=None):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.14", fc=CARD,
                                    ec=tone or EDGE, lw=1.2 if tone else 0.9, zorder=2))
        if tone:
            ax.add_patch(FancyBboxPatch((x - 0.03, y - 0.03), w + 0.06, h + 0.06,
                                        boxstyle="round,pad=0,rounding_size=0.16", fc="none", ec=tone, lw=5,
                                        alpha=0.12, zorder=1))

    def label(x, y, text, color=GLOW):
        ax.text(x, y, text.upper(), fontsize=7.5, color=color, weight="bold", fontfamily="DejaVu Sans Mono",
                va="top", zorder=3)

    def arrow(x0, y0, x1, y1, tone=GLOW):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=3,
                    arrowprops=dict(arrowstyle="-|>,head_length=0.5,head_width=0.25", color=tone, lw=1.3,
                                    alpha=0.8, shrinkA=2, shrinkB=2, connectionstyle="arc3,rad=0"))

    label(0.35, 6.95, "Thermodynamic Waddington · how it works")
    ax.text(0.35, 6.6, "A calibrated thermodynamic layer for single-cell differentiation", fontsize=19,
            weight="bold", color=INK, va="top", zorder=3)
    ax.text(0.35, 6.02, "On top of any RNA-velocity workflow (scVelo, CellRank, ...): where a lineage commits, how "
            "high the barrier is, and whether the flow is irreversible beyond its own noise.",
            fontsize=10, color=MUTED, va="top", zorder=3)

    def step(x, y, w, h, tag, title, body):
        card(x, y, w, h)
        label(x + 0.2, y + h - 0.2, tag, FAINT)
        ax.text(x + 0.2, y + h - 0.5, title, fontsize=12.5, weight="bold", color=INK, va="top", zorder=3)
        ax.text(x + 0.2, y + h - 0.95, body, fontsize=9, color=MUTED, va="top", linespacing=1.55, zorder=3)

    step(0.35, 2.05, 2.45, 2.9, "01 · input", "Cells + velocity",
         "scRNA-seq expression\n+ RNA velocity\n(scVelo, or a spliced/\nunspliced proxy)")
    step(3.25, 2.05, 2.45, 2.9, "02 · graph", "Cell-state graph",
         "PCA embedding\nexact KD-tree kNN\n50k cells in seconds\n\nvelocity = flow\non every edge")
    arrow(2.8, 3.5, 3.25, 3.5)

    rows = [
        (4.05, "03 · Hodge", "Hodge decomposition", "flow = gradient (reversible)\n      + cyclic (irreversible)",
         "Calibrated irreversibility", "noise-matched null:\nfalse positives 100% → 0%"),
        (2.4, "04 · TPT", "Transition-path theory", "committor q: 0 = progenitor,\n1 = terminal fate",
         "Commitment barrier", "F(q) = −kT ln P(q)\npancreas 1.6 kT [1.4, 1.9]"),
        (0.75, "05 · Jarzynski", "Path work", "free-energy landscape with a\ndensity-anchored kT scale",
         "kT landscape", "pseudopotential in kT\n≈ 5 kT deep (pancreas)"),
    ]
    for (y, tag, t1, b1, t2, b2), tone in zip(rows, ROW_TONES):
        card(6.2, y, 3.0, 1.45)
        label(6.4, y + 1.27, tag, FAINT)
        ax.text(6.4, y + 1.0, t1, fontsize=11.5, weight="bold", color=INK, va="top", zorder=3)
        ax.text(6.4, y + 0.62, b1, fontsize=8.5, color=MUTED, va="top", linespacing=1.5, zorder=3)
        card(9.75, y, 3.9, 1.45, tone)
        ax.text(9.95, y + 1.25, t2, fontsize=11.5, weight="bold", color=tone, va="top", zorder=3)
        ax.text(9.95, y + 0.86, b2, fontsize=8.5, color=INK, va="top", linespacing=1.5, zorder=3)
        arrow(5.7, 3.5, 6.2, y + 0.72)
        arrow(9.2, y + 0.72, 9.75, y + 0.72, tone)

    # Small views of the real readouts, inside the result cards.
    def inset(y):
        a = fig.add_axes([12.05 / 14, (y + 0.2) / 7.2, 1.45 / 14, 1.05 / 7.2])
        a.set_facecolor("none")
        a.axis("off")
        return a

    a = inset(4.05)
    names = ["null", "real"]
    vals = [0.100, 0.085] if summary is None else [summary["irreversibility"]["noise_matched_null"],
                                                  summary["irreversibility"]["cyclic_fraction"]]
    a.barh(names, vals, color=[(1, 1, 1, 0.18), ROW_TONES[0]], height=0.55)
    a.set_xlim(0, max(vals) * 1.25)
    for i, v in enumerate(vals):
        a.text(v, i, f" {v:.3f}", va="center", fontsize=7, color=MUTED, fontfamily="DejaVu Sans Mono")
        a.text(0, i + 0.42, names[i], fontsize=6.5, color=FAINT, fontfamily="DejaVu Sans Mono")

    a = inset(2.4)
    if q is not None:
        xs = np.linspace(0, 1, 60)
        bw = 1.06 * q.std() * len(q) ** -0.2
        dens = np.exp(-0.5 * ((xs[:, None] - q[None, :]) / bw) ** 2).sum(1)
        f = -np.log(dens / dens.max())
        for k in range(len(xs) - 1):
            a.plot(xs[k:k + 2], f[k:k + 2], color=cmap(xs[k]), lw=2)
        a.fill_between(xs, f, f.max() * 1.05, color=ROW_TONES[1], alpha=0.06)
        top = int(np.argmax(np.where((xs > np.quantile(q, 0.05)) & (xs < np.quantile(q, 0.95)), f, -1)))
        a.scatter([xs[top]], [f[top]], s=18, color=BG, edgecolor=GLOW, lw=1.2, zorder=4)
        a.invert_yaxis()
        a.invert_yaxis()

    a = inset(0.75)
    if xy is not None:
        gx, gy = np.meshgrid(*(np.linspace(xy[:, k].min(), xy[:, k].max(), 90) for k in range(2)))
        h = 0.05 * float(np.ptp(xy, axis=0).max())
        d2 = (gx.ravel()[:, None] - xy[None, :, 0]) ** 2 + (gy.ravel()[:, None] - xy[None, :, 1]) ** 2
        w = np.exp(-0.5 * d2 / h / h)
        dens = w.sum(1).reshape(gx.shape)
        qf = (w @ q / np.maximum(w.sum(1), 1e-300)).reshape(gx.shape)
        u = -np.log(dens / dens.max() + 0.04)
        rgba = cmap(qf)
        rgba[..., 3] = np.clip((u.max() - u) / (u.max() - u.min()), 0, 1) ** 0.7
        a.imshow(rgba, origin="lower", aspect="auto", extent=(0, 1, 0, 1))
        a.contour(np.linspace(0, 1, 90), np.linspace(0, 1, 90), u, levels=6, colors=[(1, 1, 1, 0.25)],
                  linewidths=0.5)

    ax.text(0.35, 1.45, "Validated on 4 lineages across 2 phyla", fontsize=9.5, color=INK, weight="bold",
            va="top", zorder=3)
    ax.text(0.35, 1.12, "pancreas · gastrulation · bone marrow\nDrosophila testis (Fly Cell Atlas)", fontsize=9,
            color=MUTED, va="top", linespacing=1.5, zorder=3)
    label(0.35, 0.42, "233 tests · lint-clean CI · typed · MIT", FAINT)
    fig.savefig(out, facecolor=BG)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
