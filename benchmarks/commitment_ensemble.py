"""How reproducible is the committor? An ensemble of independent fits on two real lineages.

One fit gives one committor per cell. Refit the landscape on random 80% subsets of the cells and
every cell gets a committor mean and spread; every fit gets its own commitment stage and its own
resolution-matched barrier. Two questions:

  1  Where is the committor uncertain? If the spread is concentrated in the transition state, the
     fits agree on who is committed and differ only on how far along the uncommitted cells are.
  2  Does a barrier's bootstrap error over one fit's cells describe how much it changes between
     fits? If not, single-fit barriers need the across-fit spread as their error bar.

Pancreatic endocrinogenesis (Ductal to Beta) and Drosophila spermatogenesis, 1,500 cells each.
Writes experiments/commitment_ensemble.json and figures/commitment_ensemble.png.

    python benchmarks/commitment_ensemble.py [--fits 10] [--cells 1500]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from thermodynamic_waddington import commitment_ensemble  # noqa: E402

PANCREAS = "data/real/endocrinogenesis_day15.h5ad"


def lineage(name, cells, seed=0):
    if name == "pancreas":
        import anndata as ad
        from pancreas_entropy_validation import BRANCH, load_branch
        expr, vel, labels = load_branch(PANCREAS, BRANCH, 10 ** 9, 50, 0)
        a = ad.read_h5ad(PANCREAS, backed="r")
        keep = np.isin(a.obs["clusters"].astype(str).to_numpy(), BRANCH)
        umap = np.asarray(a.obsm["X_umap"])[keep]
        expr, vel, labels, src, tgt, order = np.asarray(expr), np.asarray(vel), np.asarray(labels), ["Ductal"], ["Beta"], BRANCH
    else:
        import fly_testis as ft
        d = ft.load()
        scored = d["stage"] >= 0
        expr, vel, labels, umap = d["spliced"][scored], d["velocity"][scored], d["labels"][scored], d["umap"][scored]
        src, tgt = ft.SOURCE, ft.TARGET
        order = sorted(ft.STAGES, key=ft.STAGES.get)
    idx = np.sort(np.random.default_rng(seed).choice(len(labels), size=min(cells, len(labels)), replace=False))
    return expr[idx], vel[idx], labels[idx], umap[idx], src, tgt, order


def run(name, fits, cells):
    expr, vel, labels, umap, src, tgt, order = lineage(name, cells)
    ens = commitment_ensemble(expr, vel, labels, src, tgt, fits=fits, fraction=0.8, seed=11)
    q, sd = np.asarray(ens.q_mean), np.asarray(ens.q_sd)
    edges = np.linspace(0, 1, 11)
    mid = 0.5 * (edges[1:] + edges[:-1])
    binned = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (q >= lo) & ((q < hi) | (hi == 1)) & np.isfinite(sd)
        binned.append([round(float(np.median(sd[sel])), 4), round(float(np.percentile(sd[sel], 25)), 4),
                       round(float(np.percentile(sd[sel], 75)), 4), int(sel.sum())] if sel.sum() >= 5 else None)
    present = [s for s in order if s in set(labels.tolist())]
    stage_sd = {s: round(float(np.nanmedian(sd[labels == s])), 4) for s in present}
    stage_q = {s: round(float(np.nanmean(q[labels == s])), 4) for s in present}
    se = float(np.mean(ens.barrier_se_kt))
    out = ens.to_dict()
    out.update({
        "cells": int(len(labels)), "source": src, "target": tgt, "stage_order": present,
        "stage_median_sd": stage_sd, "stage_mean_q": stage_q,
        "sd_vs_q": {"q": mid.round(3).tolist(), "rows": binned},
        "within_fit_se_kt": round(se, 4),
        "across_to_within_ratio": round(ens.barrier_sd_kt / se, 3) if se > 0 else None,
        "transition_to_committed_sd_ratio": round(ens.sd_by_phase["transition state"] / max(ens.sd_by_phase["committed"], 1e-9), 2),
    })
    print(f"  {name}: rank agreement {ens.rank_agreement:.3f}; barrier {ens.barrier_mean_kt:.2f} +/- {ens.barrier_sd_kt:.2f} kT across fits, "
          f"within-fit se {se:.2f} kT; stage agreement {ens.commitment_agreement:.2f}; sd by phase "
          + ", ".join(f"{k} {v:.3f}" for k, v in ens.sd_by_phase.items()), flush=True)
    return out, {"umap": umap, "q": q, "sd": sd, "labels": labels}


# Figure -------------------------------------------------------------------------------------------

def plot(report, cells, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    family = "Inter" if any("Inter" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
    plt.rcParams.update({"font.family": family, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 10, "axes.titleweight": "semibold", "axes.titlelocation": "left",
                         "legend.frameon": False, "axes.edgecolor": "#9aa1ab", "xtick.color": "#5f6672", "ytick.color": "#5f6672"})
    COL = {"pancreas": "#0f7d99", "fly": "#c0392b"}
    NAME = {"pancreas": "Pancreas (Ductal → Beta)", "fly": "Fly testis (spermatogonia → spermatids)"}
    GRID = "#e3e6eb"
    fig = plt.figure(figsize=(15, 7.6))
    gs = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.3, left=0.04, right=0.965, top=0.93, bottom=0.08)
    letters = iter("abcdefgh")

    def letter(ax):
        ax.text(-0.08, 1.06, next(letters), transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")

    sd_max = max(np.nanpercentile(c["sd"], 99) for c in cells.values())
    for k, name in enumerate(("pancreas", "fly")):
        c = cells[name]
        order = np.argsort(c["q"])
        for j, (val, cmap, label, vmax) in enumerate(((c["q"], "viridis", "mean committor", 1.0),
                                                       (c["sd"], "magma_r", "committor sd across fits", sd_max))):
            ax = fig.add_subplot(gs[0, 2 * k + j])
            o = order if j == 0 else np.argsort(np.nan_to_num(c["sd"]))
            sc = ax.scatter(c["umap"][o, 0], c["umap"][o, 1], c=val[o], cmap=cmap, vmin=0, vmax=vmax, s=4, lw=0, rasterized=True)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.set_title(f"{NAME[name].split(' (')[0]}: {label}", pad=6)
            cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
            cb.outline.set_visible(False)
            cb.ax.tick_params(labelsize=7)
            letter(ax)

    ax = fig.add_subplot(gs[1, 0])
    for name in ("pancreas", "fly"):
        r = report[name]["sd_vs_q"]
        pts = [(x, row) for x, row in zip(r["q"], r["rows"]) if row]
        x = np.asarray([p[0] for p in pts])
        m, lo, hi = (np.asarray([p[1][i] for p in pts]) for i in range(3))
        ax.fill_between(x, lo, hi, color=COL[name], alpha=0.15, lw=0)
        ax.plot(x, m, "o-", color=COL[name], ms=3.5, lw=1.8, label=NAME[name].split(" (")[0])
    ax.set_xlabel("mean committor")
    ax.set_ylabel("committor sd across fits")
    ax.set_xlim(0, 1)
    ax.set_title("Uncommitted cells move ~0.2 between fits", pad=8)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_ylim(0, 0.34)
    ax.legend(fontsize=7.5, loc="upper right")
    letter(ax)

    ax = fig.add_subplot(gs[1, 1])
    rng = np.random.default_rng(1)
    for j, name in enumerate(("pancreas", "fly")):
        r = np.asarray(report[name]["pairwise_spearman"])
        ax.scatter(j + rng.uniform(-0.16, 0.16, len(r)), r, s=12, color=COL[name], alpha=0.55, lw=0)
        ax.boxplot([r], positions=[j], widths=0.5, showfliers=False, medianprops={"color": "#1f2328", "lw": 1.4},
                   boxprops={"color": "#5f6672", "lw": 0.8}, whiskerprops={"color": "#5f6672", "lw": 0.8},
                   capprops={"color": "#5f6672", "lw": 0.8})
        ax.text(j, 1.005, f"{np.mean(r):.2f}", ha="center", va="bottom", fontsize=8.5, color=COL[name], fontweight="semibold")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["pancreas", "fly testis"])
    ax.set_xlim(-0.6, 1.6)
    lo = min(min(report[n]["pairwise_spearman"]) for n in ("pancreas", "fly"))
    ax.set_ylim(max(0.0, lo - 0.08), 1.06)
    ax.set_ylabel("Spearman between two fits")
    ax.set_title("Do the fits order the cells alike?", pad=8)
    ax.text(0.5, 0.02, "uncommitted cells shared by both fits", transform=ax.transAxes, ha="center", fontsize=7, color="#5f6672")
    ax.grid(axis="y", color=GRID, lw=0.6)
    letter(ax)

    ax = fig.add_subplot(gs[1, 2])
    rng = np.random.default_rng(0)
    for j, name in enumerate(("pancreas", "fly")):
        r = report[name]
        b, se = np.asarray(r["barriers_kt"]), np.asarray(r["barrier_se_kt"])
        x = j + rng.uniform(-0.12, 0.12, len(b))
        ax.errorbar(x, b, yerr=se, fmt="o", color=COL[name], ms=4, elinewidth=0.9, capsize=2, alpha=0.8,
                    label="one fit ± bootstrap se" if j == 0 else None)
        ax.errorbar(j + 0.33, r["barrier_mean_kt"], yerr=r["barrier_sd_kt"], fmt="D", color="#1f2328", ms=5, capsize=4,
                    elinewidth=1.4, label="across fits: mean ± sd" if j == 0 else None)
        ax.text(j + 0.33, r["barrier_mean_kt"] + r["barrier_sd_kt"] + 0.12, f"{r['across_to_within_ratio']:.1f}×",
                ha="center", fontsize=8, color="#1f2328")
    ax.set_xticks([0.1, 1.1])
    ax.set_xticklabels(["pancreas", "fly testis"])
    ax.set_xlim(-0.45, 1.65)
    allb = [v for n in ("pancreas", "fly") for v in report[n]["barriers_kt"]]
    ax.set_ylim(min(allb) - 0.35, max(allb) + 0.75)
    ax.set_ylabel("resolution-matched barrier (kT)")
    ax.set_title("Between-fit vs within-fit error", pad=8)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.legend(fontsize=7, loc="upper left", ncol=2, handletextpad=0.3, columnspacing=0.8)
    letter(ax)

    ax = fig.add_subplot(gs[1, 3])
    ypos, ticks = 0, []
    for name in ("pancreas", "fly"):
        st = [s for s in report[name]["commitment_stages"]]
        counts = {s: st.count(s) for s in dict.fromkeys(st)}
        for s, cnt in counts.items():
            ax.barh(ypos, cnt / len(st), color=COL[name], alpha=0.85)
            ax.text(cnt / len(st) + 0.02, ypos, f"{cnt}/{len(st)}", va="center", fontsize=7.5)
            ticks.append((ypos, (s or "none").replace("elongation-stage", "elong.").replace("elongation stage", "elong.")))
            ypos += 1
        ypos += 0.6
    ax.set_yticks([t[0] for t in ticks])
    ax.set_yticklabels([t[1] for t in ticks], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.18)
    ax.set_xlabel("share of fits")
    ax.set_title("Commitment stage named by each fit", pad=8)
    letter(ax)

    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fits", type=int, default=10)
    ap.add_argument("--cells", type=int, default=1500)
    args = ap.parse_args()
    report, cells = {"fits": args.fits, "fraction": 0.8}, {}
    for name in ("pancreas", "fly"):
        report[name], cells[name] = run(name, args.fits, args.cells)
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/commitment_ensemble.json").write_text(json.dumps(report, indent=2) + "\n")
    plot(report, cells, "figures/commitment_ensemble.png")
    print("wrote experiments/commitment_ensemble.json and figures/commitment_ensemble.png")


if __name__ == "__main__":
    main()
