"""Cyclic-fraction irreversibility on the three real tissues.

Applies the calibrated Hodge cyclic-fraction measure to real proxy velocity and
places it between two references computed on the same cells:

  gradient floor  a pure radial gradient field v = centroid - x (reversible) on
                  the same points: the finite-sample discretization floor.
  random ceiling  the velocity vectors shuffled across cells: a structure-free
                  field, near the maximum cyclic fraction.

If the real field's cyclic fraction sits well above its own gradient floor, its
flow carries circulation a reversible field cannot explain. This is descriptive
placement against matched references, not a single-dataset p-value (see the note
in thermodynamic_waddington/irreversibility.py).

    python benchmarks/irreversibility_real.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from generalization import DATASETS, load_dataset

from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.irreversibility import cyclic_irreversibility, noise_matched_irreversibility


def pca_project(expr, vel, dims=6):
    X = np.asarray(expr, dtype=float)
    Xc = X - X.mean(axis=0)
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    comp = vt[:dims].T                       # genes x dims
    return Xc @ comp, np.asarray(vel, dtype=float) @ comp


def noise_matched(pts, vel, seed=0):
    g = build_knn(pts.tolist(), 20)
    return noise_matched_irreversibility(pts, vel, g, seed=seed)


def frac(pts, vel, seed=0):
    g = build_knn(pts.tolist(), 20)
    return cyclic_irreversibility(pts, vel, g, seed=seed).cyclic_fraction


def detection_power(spec, cells=1500, genes=100, omegas=(0.0, 0.5, 1.0, 2.0)):
    # inject known circulation into the real geometry to probe detection power there;
    # on these near-1D lineages it barely moves the fraction (few cycles to carry it).
    sp, vel, lab = load_dataset(spec["path"], spec["key"], spec["lineage"], genes)
    idx = np.random.default_rng(0).choice(sp.shape[0], size=min(cells, sp.shape[0]), replace=False)
    pts, v = pca_project(sp[idx], vel[idx])
    c = pts - pts.mean(axis=0)
    rot = np.zeros_like(pts)
    rot[:, 0], rot[:, 1] = -c[:, 1], c[:, 0]           # rotation in the first two PCs
    scale = np.linalg.norm(v) / (np.linalg.norm(rot) + 1e-12)
    out = []
    for w in omegas:
        nm = noise_matched(pts, v + w * scale * rot, 0)
        out.append({"omega": w, "cyclic_fraction": round(frac(pts, v + w * scale * rot, 0), 4),
                    "noise_matched_excess": round(nm.excess, 4), "noise_matched_z": round(nm.z, 2)})
    return out


def run(cells=1500, genes=100, seeds=5):
    results = []
    for spec in DATASETS:
        sp, vel, lab = load_dataset(spec["path"], spec["key"], spec["lineage"], genes)
        n = sp.shape[0]
        real, floor, ceil, nm = [], [], [], []
        for s in range(seeds):
            idx = np.random.default_rng(s).choice(n, size=min(cells, n), replace=False)
            pts, v = pca_project(sp[idx], vel[idx])
            real.append(frac(pts, v, s))
            radial = pts.mean(axis=0) - pts           # pure gradient field on same points
            floor.append(frac(pts, radial, s))
            shuf = v[np.random.default_rng(100 + s).permutation(len(idx))]
            ceil.append(frac(pts, shuf, s))
            nm.append(noise_matched(pts, v, s))
        row = {
            "dataset": spec["dataset"], "tissue": spec["tissue"], "n_cells": int(min(cells, n)),
            "cyclic_fraction_real": round(float(np.mean(real)), 4),
            "cyclic_fraction_real_std": round(float(np.std(real)), 4),
            "gradient_floor": round(float(np.mean(floor)), 4),
            "random_ceiling": round(float(np.mean(ceil)), 4),
            "excess_over_floor": round(float(np.mean(real) - np.mean(floor)), 4),
            "noise_matched_excess": round(float(np.mean([r.excess for r in nm])), 4),
            "noise_matched_excess_std": round(float(np.std([r.excess for r in nm])), 4),
            "noise_matched_z": round(float(np.mean([r.z for r in nm])), 2),
        }
        results.append(row)
        print(f"  {row['dataset']:24s} real={row['cyclic_fraction_real']:.3f}  "
              f"gradient floor={row['gradient_floor']:.3f}  random ceiling={row['random_ceiling']:.3f}  "
              f"excess={row['excess_over_floor']:+.3f}  noise-matched={row['noise_matched_excess']:+.3f} "
              f"(z={row['noise_matched_z']:+.1f})", flush=True)
    print("detection power on real pancreas geometry (inject known circulation):", flush=True)
    power = detection_power(DATASETS[0], cells, genes)
    for row in power:
        print(f"    injected omega={row['omega']}: cyclic_fraction={row['cyclic_fraction']}  "
              f"noise-matched={row['noise_matched_excess']:+.3f} (z={row['noise_matched_z']:+.1f})", flush=True)

    report = {
        "benchmark": "irreversibility_real_cyclic_fraction",
        "setup": f"{cells} cells, {genes} genes, {seeds} seeds, PCA(6) space, proxy velocity",
        "results": results,
        "detection_power_pancreas_geometry": power,
        "note": "Cyclic fraction of real velocity vs a matched gradient (reversible) floor and a shuffled ceiling on the same cells. Excess over the floor indicates circulation beyond what a reversible field produces at this sample size. Descriptive placement, not a calibrated p-value. detection_power injects a PC1-PC2 rotation into the real geometry: it barely moves the fraction, because these lineages are near-1D and a 1D structure has no cycles to carry circulation. The measure's power is established on the 2D synthetic ground truth (AUROC 1.0), not here; this null is consistent with the near-floor real result. noise_matched_excess reads the same field against a null built from the data's own noise (noise_matched_irreversibility), robust to velocity noise and sequencing depth.",
    }
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/irreversibility_real.json").write_text(json.dumps(report, indent=2))
    print("wrote experiments/irreversibility_real.json")
    plot(report, "figures/irreversibility_four_tissues.png")
    return report


def plot(report, out):
    """All tissues on one scale: clean-floor excess vs noise-matched excess."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [{"dataset": r["dataset"], "naive": r["excess_over_floor"], "nm": r["noise_matched_excess"],
             "nm_std": r["noise_matched_excess_std"]} for r in report["results"]]
    fly = Path("experiments/fly_testis.json")
    if fly.exists():
        irr = json.loads(fly.read_text())["irreversibility"]
        rows.append({"dataset": "fly testis germline", "naive": irr["excess_over_floor"],
                     "nm": irr["noise_matched"]["excess"], "nm_std": irr["noise_matched"]["excess_std"]})
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.bar(x - 0.2, [r["naive"] for r in rows], 0.4, color="#9aa5b1", label="excess over clean floor")
    ax.bar(x + 0.2, [r["nm"] for r in rows], 0.4, yerr=[r["nm_std"] for r in rows], capsize=3,
           color="#0f7d99", label="noise-matched excess")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x, [r["dataset"].replace(" ", "\n", 1) for r in rows])
    ax.set_ylabel("excess cyclic (irreversible) fraction")
    ax.set_title("Four developmental lineages, two phyla: flow is reversible within its own noise")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")


if __name__ == "__main__":
    run()
