"""Pipeline overview diagram (figures/tw_pipeline_overview.png).

    python benchmarks/plot_overview.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

INK, MUTED, LINE = "#1f2933", "#52606d", "#9aa5b1"
TEAL, PURPLE, AMBER = "#0f7d99", "#6a4c93", "#b7791f"


def main(out="figures/tw_pipeline_overview.png"):
    fig, ax = plt.subplots(figsize=(13, 6.4), dpi=200)
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 6.4)
    ax.axis("off")

    def box(x, y, w, h, title, body, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc="#ffffff", ec=color, lw=1.8))
        ax.text(x + w / 2, y + h - 0.3, title, ha="center", va="top", fontsize=11.5, weight="bold", color=color)
        ax.text(x + w / 2, y + h - 0.72, body, ha="center", va="top", fontsize=9, color=INK, linespacing=1.4)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=LINE, lw=1.6, shrinkA=2, shrinkB=2))

    ax.text(0.2, 6.15, "Thermodynamic Waddington Pipeline: how it works", fontsize=15, weight="bold",
            color=INK, va="top")
    ax.text(0.2, 5.7, "A thermodynamic layer on top of any RNA-velocity workflow (scVelo, CellRank, ...)",
            fontsize=10, color=MUTED, va="top")
    box(0.2, 1.9, 2.3, 2.6, "Input", "scRNA-seq\nexpression\n+ RNA velocity\n(scVelo, or a\nspliced/unspliced\nproxy)", MUTED)
    box(3.1, 1.9, 2.3, 2.6, "Cell-state graph",
        "PCA embedding\nexact KD-tree kNN\n(50k cells\nin seconds)\n\nvelocity = flow\non each edge", MUTED)
    arrow(2.5, 3.2, 3.1, 3.2)
    rows = [
        (3.85, "Hodge decomposition", "flow = gradient (reversible)\n+ cyclic (irreversible)", TEAL,
         "Calibrated irreversibility", "noise-matched null: false\npositives 100% -> 0%"),
        (2.2, "Transition-path theory", "committor q: 0 = progenitor,\n1 = terminal fate", PURPLE,
         "Commitment barrier", "PMF along q, in kT\npancreas: 1.6 kT [1.4, 1.9]"),
        (0.55, "Jarzynski path work", "free-energy landscape +\ndensity-anchored kT scale", AMBER,
         "kT landscape", "pseudopotential in kT\n~5 kT deep (pancreas)"),
    ]
    for y, t1, b1, c, t2, b2 in rows:
        box(6.0, y, 3.0, 1.4, t1, b1, c)
        box(9.6, y, 3.2, 1.4, t2, b2, c)
        arrow(5.4, 3.2, 6.0, y + 0.7)
        arrow(9.0, y + 0.7, 9.6, y + 0.7)
    ax.text(0.2, 1.25, "Validated on 4 lineages, 2 phyla:\npancreas, gastrulation, bone marrow,\nDrosophila testis "
            "(Fly Cell Atlas)", fontsize=8.5, color=MUTED, va="top", linespacing=1.4)
    ax.text(0.2, 0.35, "230+ tests · lint-clean CI · typed · MIT", fontsize=8.5, color=MUTED, va="top")
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
