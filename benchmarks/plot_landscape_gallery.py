"""Landscape gallery from the example exports.

Three figures, one row per lineage (pancreatic endocrinogenesis, Drosophila spermatogenesis):

    figures/waddington_landscape_3d.png   the Waddington surface: a density pseudo-potential
                                          U(x) = -ln(p(x)/p_max + 0.04) over the UMAP, coloured by the
                                          committor, with the q = 0.5 separatrix traced on it
    figures/commitment_flow_field.png     the committor as a field over the UMAP (kernel regression),
                                          streamlines of its gradient (the direction of commitment),
                                          isocommittor contours and the q = 0.5 separatrix
    figures/barrier_bootstrap.png         the potential of mean force along the committor with a 95%
                                          band from resampling cells, the bootstrap distribution of
                                          the barrier and transition state, and the barrier against
                                          the number of cells it is estimated from

The surface and the field are descriptive views of the exported per-cell readouts; the barrier
band measures sampling variability of one fit (cells resampled), not variability across refits.
The figures follow the manuscript's style: white ground, panel letters, the committor in viridis
(as in the Drosophila figure) and the separatrix or transition state in red.

    python benchmarks/export_examples.py      # writes examples/exports/
    python benchmarks/plot_landscape_gallery.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

EXPORTS = Path("examples/exports")
LINEAGES = [("pancreas", "Pancreatic endocrinogenesis", "Ductal", "Beta"),
            ("fly_testis", "Drosophila spermatogenesis", "spermatogonia", "elongating spermatids")]
INK, MUTED, GRID = "#1f2328", "#5f6672", "#e3e6eb"
TEAL, PURPLE, RED = "#0f7d99", "#6a4c93", "#c0392b"
FLOOR = np.array([0.93, 0.94, 0.95])          # the landscape where there are no cells
GRID_N = 140
BOOTSTRAPS = 400


def load(prefix):
    rows = list(csv.DictReader((EXPORTS / f"{prefix}_cells.csv").open()))
    summary = json.loads((EXPORTS / f"{prefix}_summary.json").read_text())
    xy = np.array([[float(r["tw_x"]), float(r["tw_y"])] for r in rows])
    q = np.array([float(r["tw_committor"]) for r in rows])
    return xy, q, summary


def grid_over(xy, pad=0.08, n=GRID_N):
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    span = hi - lo
    gx = np.linspace(lo[0] - pad * span[0], hi[0] + pad * span[0], n)
    gy = np.linspace(lo[1] - pad * span[1], hi[1] + pad * span[1], n)
    return np.meshgrid(gx, gy)


def kernel_fields(xy, q, X, Y, scale=0.035):
    """Density and committor fields on the grid: Gaussian kernel density and Nadaraya-Watson regression."""
    h = scale * float(np.ptp(xy, axis=0).max())
    pts = np.stack([X.ravel(), Y.ravel()], axis=1)
    dens = np.zeros(len(pts))
    num = np.zeros(len(pts))
    for start in range(0, len(pts), 2000):
        chunk = pts[start:start + 2000]
        d2 = ((chunk[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2)
        w = np.exp(-0.5 * d2 / (h * h))
        dens[start:start + 2000] = w.sum(axis=1)
        num[start:start + 2000] = w @ q
    qf = num / np.maximum(dens, 1e-300)
    dens = dens / (len(xy) * 2 * np.pi * h * h)
    return dens.reshape(X.shape), qf.reshape(X.shape)


def setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    family = "Inter" if any("Inter" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
    plt.rcParams.update({
        "font.family": family, "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "semibold",
        "axes.titlelocation": "left", "axes.labelsize": 9, "axes.labelcolor": INK, "axes.edgecolor": "#9aa1ab",
        "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white",
        "savefig.facecolor": "white", "legend.frameon": False,
    })
    return plt


def letter(fig, x, y, text):
    fig.text(x, y, text, fontsize=13, weight="bold", color=INK, va="top")


def colorbar(fig, rect, label):
    import matplotlib as mpl
    ax = fig.add_axes(rect)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1), cmap="viridis"), cax=ax,
                      orientation="horizontal")
    cb.outline.set_linewidth(0.5)
    cb.set_ticks([0, 0.5, 1])
    cb.set_ticklabels(["0 (source)", "0.5", "1 (target)"])
    ax.tick_params(labelsize=8, length=2)
    ax.set_title(label, fontsize=8.5, loc="center", pad=3, weight="normal", color=INK)


def landscape_figure(plt):
    """The Waddington-surface figure, both lineages; returns the figure and its two 3D axes."""
    import contourpy
    from matplotlib.colors import LightSource
    from scipy.interpolate import RegularGridInterpolator

    cm = plt.get_cmap("viridis")
    fig = plt.figure(figsize=(12, 5.6), dpi=220)
    axes = []
    for col, (prefix, title, src, tgt) in enumerate(LINEAGES):
        xy, q, s = load(prefix)
        X, Y = grid_over(xy, pad=0.12)
        dens, qf = kernel_fields(xy, q, X, Y, scale=0.032)
        support = dens / dens.max()
        # A density floor keeps the walls soft: a continuous terrain whose valleys are the regions holding
        # cells, rising to a plateau of about -ln(0.04) = 3.2 where there are none.
        U = -np.log(support + 0.04)
        U = U - U.min()
        cap = float(U.max())
        ax = fig.add_axes([0.0 + 0.5 * col, 0.08, 0.5, 0.8], projection="3d")
        axes.append(ax)
        mix = np.clip(support / 0.08, 0, 1) ** 0.6
        rgb = cm(np.clip(qf, 0, 1))[..., :3] * mix[..., None] + FLOOR * (1 - mix[..., None])
        shaded = LightSource(azdeg=300, altdeg=42).shade_rgb(rgb, U, blend_mode="soft", vert_exag=3.0)
        ax.plot_surface(X, Y, U, facecolors=shaded, rstride=1, cstride=1, linewidth=0, antialiased=True, shade=False)
        height = RegularGridInterpolator((Y[:, 0], X[0]), U)
        for line in contourpy.contour_generator(X, Y, np.where(support > 0.05, qf, np.nan)).lines(0.5):
            if len(line) >= 4:
                ax.plot(line[:, 0], line[:, 1], height(line[:, ::-1]) - 0.05, color=RED, lw=1.6, zorder=10)
        ax.view_init(elev=56, azim=-64)
        ax.set_zlim(-0.1, cap + 0.2)
        ax.set_box_aspect((1, 1, 0.42))
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_pane_color((1, 1, 1, 0))
            axis.line.set_color("#c4c9d1")
            axis._axinfo["grid"]["color"] = (0.85, 0.87, 0.9, 1)
            axis._axinfo["grid"]["linewidth"] = 0.5
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticks([0, 1, 2, 3])
        ax.tick_params(labelsize=7, pad=0, colors=MUTED)
        ax.set_xlabel("UMAP 1", fontsize=8, labelpad=-10, color=MUTED)
        ax.set_ylabel("UMAP 2", fontsize=8, labelpad=-10, color=MUTED)
        ax.set_zlabel("U = −ln(p/p_max + 0.04)", fontsize=8, labelpad=2, color=MUTED)
        c = s["commitment"]
        letter(fig, 0.02 + 0.5 * col, 0.97, "ab"[col])
        fig.text(0.05 + 0.5 * col, 0.965, title, fontsize=10.5, weight="semibold", color=INK, va="top")
        fig.text(0.05 + 0.5 * col, 0.925, f"{src} → {tgt}; {len(q):,} cells; commits at {c['commitment_stage']}",
                 fontsize=8.5, color=MUTED, va="top")
    colorbar(fig, [0.4, 0.07, 0.2, 0.018], "committor q (surface colour); red: q = 0.5 separatrix")
    return fig, axes


def landscape_3d(plt, out):
    fig, _ = landscape_figure(plt)
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"wrote {out}")


def flow_field(plt, out):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), dpi=220, gridspec_kw={"wspace": 0.08})
    fig.subplots_adjust(left=0.03, right=0.97, top=0.86, bottom=0.17)
    cm = plt.get_cmap("viridis")
    for col, (prefix, title, src, tgt) in enumerate(LINEAGES):
        ax = axes[col]
        xy, q, s = load(prefix)
        X, Y = grid_over(xy)
        dens, qf = kernel_fields(xy, q, X, Y)
        support = dens / dens.max()
        inside = support > 0.02
        rgba = cm(np.clip(qf, 0, 1))
        rgba[..., 3] = np.clip(support ** 0.45, 0, 1) * 0.35
        ax.imshow(rgba, origin="lower", extent=(X.min(), X.max(), Y.min(), Y.max()), aspect="auto",
                  interpolation="bilinear")
        gy, gx = np.gradient(qf, Y[:, 0], X[0])
        mag = np.hypot(gx, gy)
        gx = np.where(inside, gx, np.nan)
        gy = np.where(inside, gy, np.nan)
        width = 0.3 + 1.3 * np.nan_to_num(mag / np.nanpercentile(mag[inside], 95)).clip(0, 1) * support.clip(0, 1) ** 0.3
        ax.streamplot(X, Y, gx, gy, color="#3a3f47", linewidth=width, density=1.4, arrowsize=0.6,
                      arrowstyle="-|>", broken_streamlines=True)
        masked = np.where(inside, qf, np.nan)
        cs = ax.contour(X, Y, masked, levels=[0.2, 0.4, 0.6, 0.8], colors="#8a919c", linewidths=0.5)
        ax.clabel(cs, fontsize=6.5, fmt="%.1f", inline_spacing=2)
        ax.contour(X, Y, masked, levels=[0.5], colors=RED, linewidths=1.6)
        ax.scatter(xy[:, 0], xy[:, 1], c=q, cmap=cm, vmin=0, vmax=1, s=1.6, lw=0, alpha=0.9, rasterized=True)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(True)
        for spine in ax.spines.values():
            spine.set_color("#c4c9d1")
        c = s["commitment"]
        ax.set_title(title, pad=18)
        ax.text(0, 1.012, f"{src} → {tgt}; commits at {c['commitment_stage']}", transform=ax.transAxes,
                fontsize=8.5, color=MUTED, va="bottom")
        ax.text(-0.025, 1.09, "ab"[col], transform=ax.transAxes, fontsize=13, weight="bold", color=INK, va="bottom")
        ax.set_xlabel("UMAP 1", color=MUTED, fontsize=8)
        ax.set_ylabel("UMAP 2", color=MUTED, fontsize=8)
    colorbar(fig, [0.4, 0.07, 0.2, 0.022], "committor q; grey: isocommittors; red: q = 0.5; lines: gradient of q")
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"wrote {out}")


def barrier_bootstrap(plt, out):
    from thermodynamic_waddington.developmental import committor_free_energy_profile

    fig = plt.figure(figsize=(12, 7.6), dpi=220)
    rng = np.random.default_rng(0)
    for row, (prefix, title, src, tgt) in enumerate(LINEAGES):
        xy, q, s = load(prefix)
        base = committor_free_energy_profile(q)
        xs = np.asarray(base.q_grid)
        curves, heights, places = [], [], []
        for _ in range(BOOTSTRAPS):
            p = committor_free_energy_profile(q[rng.integers(0, len(q), len(q))])
            curves.append(p.free_energy_kt)
            heights.append(p.barrier_kt)
            places.append(p.barrier_q)
        curves = np.asarray(curves)
        lo, hi = np.percentile(curves, [2.5, 97.5], axis=0)
        h_lo, h_hi = np.percentile(heights, [2.5, 97.5])
        top = 0.9 - row * 0.47
        height = 0.33
        letters = "abc" if row == 0 else "def"

        ax = fig.add_axes([0.07, top - height, 0.4, height])
        for k in range(5):
            if k % 2:
                ax.axvspan(0.2 * k, 0.2 * (k + 1), color="#f1f3f6", lw=0)
        ax.fill_between(xs, lo, hi, color=TEAL, alpha=0.18, lw=0, label="95% band (cells resampled)")
        ax.plot(xs, base.free_energy_kt, color=TEAL, lw=2, label="F(q), one fit")
        ax.axvline(base.barrier_q, color=RED, lw=0.9, ls="--")
        ax.plot([base.barrier_q], [base.barrier_kt], "o", ms=5, mfc="white", mec=RED, mew=1.4)
        ax.annotate(f"{base.barrier_kt:.2f} kT [{h_lo:.2f}, {h_hi:.2f}]", (base.barrier_q, base.barrier_kt),
                    xytext=(8, 6), textcoords="offset points", fontsize=8.5, color=INK)
        ax.set_xlim(0, 1)
        ylim = (0, max(hi.max(), base.barrier_kt) * 1.2)
        ax.set_ylim(ylim)
        ax.set_xlabel("committor q")
        ax.set_ylabel("free energy F(q) (kT)")
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.legend(loc="upper right", fontsize=7.5)
        ax.set_title(f"{title}: one fit, {len(q):,} cells", pad=6)
        ax.text(-0.13, 1.1, letters[0], transform=ax.transAxes, fontsize=13, weight="bold", color=INK, va="bottom")

        hx = fig.add_axes([0.51, top - height, 0.12, height])
        hx.hist(heights, bins=24, orientation="horizontal", color="#7fb3c4", lw=0)
        hx.axhline(base.barrier_kt, color=RED, lw=0.9)
        hx.set_ylim(ylim)
        hx.set_yticklabels([])
        hx.set_xlabel("resamples")
        hx.set_title("barrier, resampled", pad=6, weight="normal", fontsize=9)
        hx.text(0.95, 0.04, f"TS q = {np.median(places):.2f}\n[{np.percentile(places, 2.5):.2f}, "
                f"{np.percentile(places, 97.5):.2f}]", transform=hx.transAxes, ha="right", va="bottom",
                fontsize=7.5, color=MUTED)
        hx.text(-0.12, 1.1, letters[1], transform=hx.transAxes, fontsize=13, weight="bold", color=INK, va="bottom")

        rx = fig.add_axes([0.71, top - height, 0.26, height])
        sizes = [n for n in (150, 250, 400, 600, 900, 1500, 2200, 3000) if n <= len(q)]
        means, lows, highs = [], [], []
        for n in sizes:
            b = [committor_free_energy_profile(q[rng.choice(len(q), n, replace=False)]).barrier_kt
                 for _ in range(60 if n < len(q) else 1)]
            means.append(np.mean(b))
            lows.append(np.percentile(b, 2.5))
            highs.append(np.percentile(b, 97.5))
        rx.fill_between(sizes, lows, highs, color=PURPLE, alpha=0.15, lw=0)
        rx.plot(sizes, means, color=PURPLE, lw=1.8, marker="o", ms=3.5, mfc="white", mec=PURPLE)
        rx.set_xscale("log")
        rx.set_xticks(sizes)
        rx.set_xticklabels([f"{n:,}" for n in sizes], fontsize=7)
        rx.minorticks_off()
        rx.set_ylim(0, max(highs) * 1.15)
        rx.set_xlabel("cells in the estimate (log scale)")
        rx.set_ylabel("barrier (kT)")
        rx.grid(color=GRID, lw=0.6)
        rx.set_title("barrier vs cell number", pad=6, weight="normal", fontsize=9)
        rx.text(-0.2, 1.1, letters[2], transform=rx.transAxes, fontsize=13, weight="bold", color=INK, va="bottom")
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    plt = setup()
    landscape_3d(plt, "figures/waddington_landscape_3d.png")
    flow_field(plt, "figures/commitment_flow_field.png")
    barrier_bootstrap(plt, "figures/barrier_bootstrap.png")


if __name__ == "__main__":
    main()
