"""Noise-matched irreversibility on ground truth with realistic, uneven noise.

Same linear diffusion as synthetic_ground_truth.py (dx = -Bx dt + noise,
B = I + w J; zero entropy production at w = 0), now with measurement noise added
to the velocity: none, uniform, or heteroscedastic (one half of the cells 10x
noisier than the other, the way shallow-sequenced cells are). Three readouts:

  raw cyclic fraction        the Hodge measure alone
  naive excess               cyclic fraction minus the clean radial floor
  noise-matched excess       noise_matched_irreversibility (null built from the
                             data's own noise, bias-corrected on references)

Positives are w > 0, negatives w = 0, pooled across noise levels: a measure
that mistakes noise for circulation cannot separate noisy equilibrium from
clean driven fields. Also reported: the false-positive rate at equilibrium
(noise-matched z > 2) and the detection rate per drive and noise level.

    python benchmarks/noise_matched_ground_truth.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.irreversibility import (
    cyclic_fraction_per_cell,
    cyclic_irreversibility,
    noise_matched_irreversibility,
)

J = np.array([[0.0, -1.0], [1.0, 0.0]])
NOISE = {"none": (0.0, False), "uniform": (0.5, False), "uneven": (0.5, True), "uneven-high": (1.0, True)}
OMEGAS = (0.0, 0.25, 0.5, 1.0, 2.0)


def sample(omega, seed, n, noise, uneven):
    rng = np.random.default_rng(seed)
    pts = rng.normal(0.0, 1.0, size=(n, 2))
    vel = -(pts @ (np.eye(2) + omega * J).T)
    sd = noise * (np.where(pts[:, 0] > 0, 3.0, 0.3) if uneven else np.ones(n))
    return pts, vel + rng.normal(size=vel.shape) * sd[:, None]


def auroc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    wins = (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
    return float(wins / (len(pos) * len(neg)))


def run(seeds, n):
    rows = []
    for name, (noise, uneven) in NOISE.items():
        for w in OMEGAS:
            for s in range(seeds):
                pts, v = sample(w, s, n, noise, uneven)
                g = build_knn(pts.tolist(), 20)
                raw = cyclic_irreversibility(pts, v, g).cyclic_fraction
                floor = cyclic_irreversibility(pts, -pts, g).cyclic_fraction
                nm = noise_matched_irreversibility(pts, v, g, seed=s)
                rows.append({"noise": name, "omega": w, "seed": s, "raw": raw, "naive_excess": raw - floor,
                             "nm_excess": nm.excess, "nm_z": nm.z})
            sub = [r for r in rows if r["noise"] == name and r["omega"] == w]
            print(f"  noise={name:12s} w={w:<4}  raw={np.mean([r['raw'] for r in sub]):.3f}  "
                  f"naive={np.mean([r['naive_excess'] for r in sub]):+.3f}  "
                  f"noise-matched={np.mean([r['nm_excess'] for r in sub]):+.3f} (z={np.mean([r['nm_z'] for r in sub]):+.1f})",
                  flush=True)
    neg = [r for r in rows if r["omega"] == 0.0]
    pos = [r for r in rows if r["omega"] > 0.0]
    pooled = {k: round(auroc([r[k] for r in pos], [r[k] for r in neg]), 3) for k in ("raw", "naive_excess", "nm_excess")}
    # how often each readout ranks a noisy equilibrium above the cleanest driven field
    clean_driven = [r for r in pos if r["noise"] == "none" and r["omega"] == 0.25]
    noisy_eq = [r for r in neg if r["noise"] in ("uneven", "uneven-high")]
    confusion = {k: round(1.0 - auroc([r[k] for r in clean_driven], [r[k] for r in noisy_eq]), 3)
                 for k in ("raw", "naive_excess", "nm_excess")}
    fpr = {name: round(float(np.mean([r["nm_z"] > 2 for r in neg if r["noise"] == name])), 3) for name in NOISE}
    naive_fpr = {name: round(float(np.mean([r["naive_excess"] > 0.02 for r in neg if r["noise"] == name])), 3) for name in NOISE}
    power = {name: {str(w): round(float(np.mean([r["nm_z"] > 2 for r in pos if r["noise"] == name and r["omega"] == w])), 3)
                    for w in OMEGAS[1:]} for name in NOISE}
    per_cell = []
    for name, (noise, uneven) in NOISE.items():
        if not uneven:
            continue
        for w in (0.0, 1.0):
            pts, v = sample(w, 0, n, noise, uneven)
            g = build_knn(pts.tolist(), 20)
            naive_local = cyclic_fraction_per_cell(pts, v, g) - cyclic_fraction_per_cell(pts, -pts, g)
            nm = noise_matched_irreversibility(pts, v, g, per_cell=True)
            noisy = pts[:, 0] > 0
            per_cell.append({"noise": name, "omega": w,
                             "naive_local_excess_noisy_half": round(float(naive_local[noisy].mean()), 3),
                             "naive_local_excess_clean_half": round(float(naive_local[~noisy].mean()), 3),
                             "nm_local_excess_noisy_half": round(float(nm.per_cell_excess[noisy].mean()), 3),
                             "nm_local_excess_clean_half": round(float(nm.per_cell_excess[~noisy].mean()), 3)})
    summary = {}
    for name in NOISE:
        for w in OMEGAS:
            sub = [r for r in rows if r["noise"] == name and r["omega"] == w]
            summary[f"{name}|{w}"] = {k: round(float(np.mean([r[k] for r in sub])), 4) for k in ("raw", "naive_excess", "nm_excess", "nm_z")}
    return {
        "benchmark": "noise_matched_ground_truth",
        "model": "dx=-Bx dt, B=I+wJ, velocity + measurement noise (none / uniform / uneven 10x)",
        "n_cells": n, "seeds": seeds,
        "pooled_auroc_equilibrium_vs_driven": pooled,
        "p_noisy_equilibrium_ranked_above_clean_weak_drive": confusion,
        "false_positive_rate_nm_z_gt_2": fpr,
        "false_positive_rate_naive_excess_gt_0.02": naive_fpr,
        "detection_rate_nm_z_gt_2": power,
        "per_cell_uneven_noise": per_cell,
        "summary": summary,
    }


def plot(report, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    teal, grey, red = "#0f7d99", "#9aa5b1", "#b0413e"
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    s = report["summary"]
    for name, color, ls in (("none", teal, "-"), ("uniform", teal, "--"), ("uneven", red, "-"), ("uneven-high", red, "--")):
        ax[0].plot(OMEGAS, [s[f"{name}|{w}"]["naive_excess"] for w in OMEGAS], ls, color=color, marker="o", label=name)
        ax[1].plot(OMEGAS, [s[f"{name}|{w}"]["nm_excess"] for w in OMEGAS], ls, color=color, marker="o", label=name)
    for a, title in ((ax[0], "naive excess (vs clean radial floor)"), (ax[1], "noise-matched excess")):
        a.axhline(0, color=grey, lw=1)
        a.set_xlabel("rotational drive w (0 = equilibrium)")
        a.set_ylabel("excess cyclic fraction")
        a.set_title(title)
        a.grid(alpha=0.3)
        a.legend(title="velocity noise", frameon=False, fontsize=9)
    ax[0].set_ylim(ax[1].get_ylim()[0], max(ax[0].get_ylim()[1], ax[1].get_ylim()[1]))
    ax[1].set_ylim(ax[0].get_ylim())
    p = report["pooled_auroc_equilibrium_vs_driven"]
    names = ["raw cyclic\nfraction", "naive\nexcess", "noise-matched\nexcess"]
    vals = [p["raw"], p["naive_excess"], p["nm_excess"]]
    ax[2].bar(names, vals, color=[grey, grey, teal])
    for i, v in enumerate(vals):
        ax[2].text(i, v + 0.01, f"{v:.2f}", ha="center")
    ax[2].axhline(0.5, color="k", ls="--", lw=1)
    ax[2].set_ylim(0, 1.05)
    ax[2].set_ylabel("AUROC, equilibrium vs driven (pooled over noise)")
    fpr = report["false_positive_rate_nm_z_gt_2"]
    ax[2].set_title(f"Pooled discrimination; noise-matched FPR at w=0: {max(fpr.values()):.2f} max")
    fig.suptitle("Noise mistaken for circulation, and the noise-matched null that removes it", fontsize=13, weight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--cells", type=int, default=400)
    args = ap.parse_args()
    report = run(args.seeds, args.cells)
    print("pooled AUROC:", report["pooled_auroc_equilibrium_vs_driven"])
    print("P(noisy equilibrium ranked above clean weak drive):", report["p_noisy_equilibrium_ranked_above_clean_weak_drive"])
    print("FPR noise-matched z>2:", report["false_positive_rate_nm_z_gt_2"])
    print("FPR naive excess>0.02:", report["false_positive_rate_naive_excess_gt_0.02"])
    print("detection rate:", report["detection_rate_nm_z_gt_2"])
    for r in report["per_cell_uneven_noise"]:
        print("  per-cell", r)
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/noise_matched_ground_truth.json").write_text(json.dumps(report, indent=2))
    plot(report, "figures/noise_matched_ground_truth.png")
    print("wrote experiments/noise_matched_ground_truth.json and figures/noise_matched_ground_truth.png")


if __name__ == "__main__":
    main()
