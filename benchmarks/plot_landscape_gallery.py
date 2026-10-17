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
BG, PANEL, INK, MUTED, FAINT = "#07080c", "#0c0e14", "#eef0f6", "#8b90a0", "#5b6070"
GRID = (1, 1, 1, 0.07)
TONES = ["#67e8f9", "#5b9dff", "#818cf8", "#a78bfa", "#e879f9"]   # one tone per commitment phase
GLOW = "#a5dcff"
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


def cmap():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("commitment", TONES)


def setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    family = "Inter" if any("Inter" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
    plt.rcParams.update({
        "font.family": family, "figure.facecolor": BG, "axes.facecolor": PANEL, "savefig.facecolor": BG,
        "text.color": INK, "axes.labelcolor": MUTED, "xtick.color": FAINT, "ytick.color": FAINT,
        "axes.edgecolor": GRID, "axes.titlecolor": INK, "axes.titleweight": "semibold", "axes.titlesize": 12,
        "axes.titlelocation": "left", "axes.titlepad": 10, "font.size": 9,
    })
    return plt


def eyebrow(fig, text, title, sub):
    fig.text(0.035, 0.965, text.upper(), color=GLOW, fontsize=8.5, weight="bold", va="top",
             fontfamily="DejaVu Sans Mono")
    fig.text(0.035, 0.935, title, color=INK, fontsize=19, weight="bold", va="top")
    fig.text(0.035, 0.885, sub, color=MUTED, fontsize=10, va="top")


def colorbar(fig, rect, label):
    import matplotlib as mpl
    ax = fig.add_axes(rect)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1), cmap=cmap()), cax=ax,
                      orientation="horizontal")
    cb.outline.set_visible(False)
    cb.set_ticks([0, 0.5, 1])
    cb.set_ticklabels(["0 · source", "0.5", "1 · target"])
    ax.tick_params(colors=MUTED, labelsize=8, length=0)
    ax.set_title(label, color=MUTED, fontsize=8, loc="left", pad=4)


def landscape_3d(plt, out):
    fig = plt.figure(figsize=(15, 7.6), dpi=200)
    eyebrow(fig, "Thermodynamic Waddington · landscape", "The Waddington surface, measured",
            "Density pseudo-potential U = −ln(p/p_max + 0.04) over the UMAP (deeper = more cells), coloured by the committor; the white line "
            "is the q = 0.5 separatrix, traced where cells are.")
    cm = cmap()
    for col, (prefix, title, src, tgt) in enumerate(LINEAGES):
        xy, q, s = load(prefix)
        X, Y = grid_over(xy, pad=0.12)
        dens, qf = kernel_fields(xy, q, X, Y, scale=0.032)
        support = dens / dens.max()
        # A density floor keeps the walls soft: the surface is a continuous terrain whose valleys are the
        # regions holding cells, rising to a plateau of about -ln(0.04) = 3.2 where there are none.
        U = -np.log(support + 0.04)
        U = U - U.min()
        cap = float(U.max())
        ax = fig.add_axes([0.02 + 0.49 * col, 0.04, 0.47, 0.8], projection="3d")
        ax.set_facecolor(BG)
        from matplotlib.colors import LightSource
        slate = np.array([0.10, 0.12, 0.18])
        mix = np.clip(support / 0.08, 0, 1) ** 0.6
        base_rgb = cm(np.clip(qf, 0, 1))[..., :3] * mix[..., None] + slate * (1 - mix[..., None])
        shaded = LightSource(azdeg=300, altdeg=40).shade_rgb(base_rgb, U, blend_mode="soft", vert_exag=3.0)
        ax.plot_surface(X, Y, U, facecolors=shaded, rstride=1, cstride=1, linewidth=0, antialiased=True, shade=False)
        # The separatrix q = 0.5, traced on the surface where cells are.
        import contourpy
        from scipy.interpolate import RegularGridInterpolator
        height = RegularGridInterpolator((Y[:, 0], X[0]), U)
        for line in contourpy.contour_generator(X, Y, np.where(support > 0.05, qf, np.nan)).lines(0.5):
            if len(line) < 4:
                continue
            z = height(line[:, ::-1]) - 0.06
            for lw, alpha in ((6, 0.14), (1.8, 1.0)):
                ax.plot(line[:, 0], line[:, 1], z, color=(0.88, 0.95, 1.0, alpha), lw=lw, zorder=10)
        ax.view_init(elev=56, azim=-64)
        ax.set_zlim(-0.1, cap + 0.2)
        ax.set_box_aspect((1, 1, 0.42))
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.set_pane_color((0, 0, 0, 0))
            axis.line.set_color((1, 1, 1, 0.08))
            axis._axinfo["grid"]["color"] = (1, 1, 1, 0.05)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticks([0, 1, 2, 3])
        ax.tick_params(colors=FAINT, labelsize=7, pad=0)
        ax.set_zlabel("U = −ln p  (valleys = where cells are)", color=FAINT, fontsize=8, labelpad=4)
        c = s["commitment"]
        fig.text(0.05 + 0.49 * col, 0.82, title, color=INK, fontsize=13, weight="bold")
        fig.text(0.05 + 0.49 * col, 0.795, f"{src} → {tgt} · {len(q):,} cells · commits at {c['commitment_stage']} · "
                 f"barrier {c['barrier_kt']:.2f} kT at q = {c['transition_state_q']:.2f}", color=MUTED, fontsize=9)
    colorbar(fig, [0.39, 0.06, 0.22, 0.014], "committor q")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


def flow_field(plt, out):
    fig = plt.figure(figsize=(15, 7.4), dpi=200)
    eyebrow(fig, "Thermodynamic Waddington · commitment field", "Where cells commit, and which way",
            "Committor as a field over the UMAP (kernel regression); streamlines follow its gradient, the direction of "
            "commitment; contours every 0.2; the bright line is the q = 0.5 separatrix.")
    cm = cmap()
    for col, (prefix, title, src, tgt) in enumerate(LINEAGES):
        xy, q, s = load(prefix)
        X, Y = grid_over(xy)
        dens, qf = kernel_fields(xy, q, X, Y)
        support = dens / dens.max()
        inside = support > 0.02
        ax = fig.add_axes([0.035 + 0.485 * col, 0.11, 0.45, 0.66])
        ax.set_facecolor(BG)
        rgba = cm(np.clip(qf, 0, 1))
        rgba[..., 3] = np.clip(support ** 0.45, 0, 1) * 0.55
        ax.imshow(rgba, origin="lower", extent=(X.min(), X.max(), Y.min(), Y.max()), aspect="auto",
                  interpolation="bilinear")
        gy, gx = np.gradient(qf, Y[:, 0], X[0])
        mag = np.hypot(gx, gy)
        gx = np.where(inside, gx, np.nan)
        gy = np.where(inside, gy, np.nan)
        width = 0.35 + 1.6 * np.nan_to_num(mag / np.nanpercentile(mag[inside], 95)).clip(0, 1) * support.clip(0, 1) ** 0.3
        ax.streamplot(X, Y, gx, gy, color=(0.85, 0.92, 1.0, 0.55), linewidth=width, density=1.5, arrowsize=0.7,
                      arrowstyle="-|>", broken_streamlines=True)
        masked = np.where(inside, qf, np.nan)
        ax.contour(X, Y, masked, levels=[0.2, 0.4, 0.6, 0.8], colors=[(1, 1, 1, 0.22)], linewidths=0.7)
        for lw, alpha in ((7, 0.06), (4, 0.12), (1.6, 0.95)):
            ax.contour(X, Y, masked, levels=[0.5], colors=[(0.78, 0.9, 1.0, alpha)], linewidths=lw)
        ax.scatter(xy[:, 0], xy[:, 1], c=q, cmap=cm, vmin=0, vmax=1, s=2.2, lw=0, alpha=0.85, rasterized=True)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(GRID)
        c = s["commitment"]
        ax.set_title(title, fontsize=13, pad=26)
        ax.text(0, 1.015, f"{src} → {tgt} · commits at {c['commitment_stage']}", transform=ax.transAxes,
                color=MUTED, fontsize=9, va="bottom")
    colorbar(fig, [0.39, 0.055, 0.22, 0.014], "committor q")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


def barrier_bootstrap(plt, out):
    from matplotlib.collections import LineCollection

    from thermodynamic_waddington.developmental import committor_free_energy_profile

    fig = plt.figure(figsize=(15, 10.2), dpi=200)
    eyebrow(fig, "Thermodynamic Waddington · commitment barrier", "How sure is the barrier?",
            f"Left: F(q) with a 95% band from {BOOTSTRAPS} cell resamples of one fit. Middle: the resampled barrier. "
            "Right: the barrier against the number of cells it is estimated from.\nThe KDE bandwidth shrinks as "
            "n^(-1/5), so the barrier rises with cell number: barriers are comparable only at matched cell numbers.")
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
        top = 0.72 - row * 0.41
        height = 0.26

        ax = fig.add_axes([0.05, top - height, 0.42, height])
        for k, (a, b) in enumerate([(0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]):
            ax.axvspan(a, b, color=TONES[k], alpha=0.06, lw=0)
        ax.fill_between(xs, lo, hi, color="#5b9dff", alpha=0.22, lw=0)
        pts = np.stack([xs, np.asarray(base.free_energy_kt)], axis=1)
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        for lw, alpha in ((8, 0.08), (4.5, 0.16), (2.2, 1.0)):
            lc = LineCollection(segs, cmap=cmap(), linewidths=lw, alpha=alpha)
            lc.set_array((xs[:-1] + xs[1:]) / 2)
            ax.add_collection(lc)
        ax.axvline(base.barrier_q, color=GLOW, lw=0.9, ls=(0, (2, 3)), alpha=0.8)
        ax.scatter([base.barrier_q], [base.barrier_kt], s=60, color=BG, edgecolor=GLOW, lw=1.8, zorder=5)
        ax.annotate(f"‡ {base.barrier_kt:.2f} kT  [{h_lo:.2f}, {h_hi:.2f}]", (base.barrier_q, base.barrier_kt),
                    xytext=(10, 8), textcoords="offset points", color=INK, fontsize=9, weight="bold")
        ax.set_xlim(0, 1)
        ylim = (0, max(hi.max(), base.barrier_kt) * 1.18)
        ax.set_ylim(ylim)
        ax.set_xlabel("committor q")
        ax.set_ylabel("F(q)  [kT]")
        ax.grid(color=GRID, lw=0.6)
        ax.set_title(title, fontsize=13, pad=22)
        ax.text(0, 1.015, f"{src} → {tgt} · one fit, {len(q):,} cells", transform=ax.transAxes, color=MUTED,
                fontsize=9, va="bottom")

        hx = fig.add_axes([0.5, top - height, 0.13, height])
        hx.hist(heights, bins=24, orientation="horizontal", color="#818cf8", alpha=0.85, lw=0)
        hx.axhline(base.barrier_kt, color=GLOW, lw=1)
        hx.set_ylim(ylim)
        hx.set_yticklabels([])
        hx.set_xlabel("resamples")
        hx.grid(color=GRID, lw=0.6)
        hx.set_title("barrier, resampled", fontsize=9.5, color=MUTED, weight="normal", pad=6)
        hx.text(0.96, 0.03, f"TS q = {np.median(places):.2f}\n[{np.percentile(places, 2.5):.2f}, "
                f"{np.percentile(places, 97.5):.2f}]", transform=hx.transAxes, ha="right", va="bottom",
                color=MUTED, fontsize=8, fontfamily="DejaVu Sans Mono")

        rx = fig.add_axes([0.69, top - height, 0.27, height])
        sizes = [n for n in (150, 250, 400, 600, 900, 1500, 2200, 3000) if n <= len(q)]
        means, lows, highs = [], [], []
        for n in sizes:
            b = [committor_free_energy_profile(q[rng.choice(len(q), n, replace=False)]).barrier_kt
                 for _ in range(60 if n < len(q) else 1)]
            means.append(np.mean(b))
            lows.append(np.percentile(b, 2.5))
            highs.append(np.percentile(b, 97.5))
        rx.fill_between(sizes, lows, highs, color="#a78bfa", alpha=0.22, lw=0)
        for lw, alpha in ((6, 0.12), (2, 1.0)):
            rx.plot(sizes, means, color="#a78bfa", lw=lw, alpha=alpha)
        rx.scatter(sizes, means, s=16, color=BG, edgecolor="#e879f9", lw=1.2, zorder=4)
        rx.set_xscale("log")
        rx.set_xticks(sizes)
        rx.set_xticklabels([f"{n:,}" for n in sizes], rotation=0, fontsize=7)
        rx.minorticks_off()
        rx.set_ylim(0, max(highs) * 1.15)
        rx.set_xlabel("cells in the estimate")
        rx.set_ylabel("barrier  [kT]")
        rx.grid(color=GRID, lw=0.6)
        rx.set_title("barrier vs cell number", fontsize=9.5, color=MUTED, weight="normal", pad=6)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    plt = setup()
    landscape_3d(plt, "figures/waddington_landscape_3d.png")
    flow_field(plt, "figures/commitment_flow_field.png")
    barrier_bootstrap(plt, "figures/barrier_bootstrap.png")


if __name__ == "__main__":
    main()
