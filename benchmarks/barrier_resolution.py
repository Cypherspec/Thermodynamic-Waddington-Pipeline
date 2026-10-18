"""Is the commitment barrier comparable across datasets of different sizes?

The potential of mean force F(q) = -kT ln P(q) along the committor is estimated with a Gaussian
KDE. With Silverman's rule the bandwidth shrinks as n^(-1/5), and the committor distribution is
far from Gaussian: every cell in the source or target set sits at exactly q = 0 or q = 1, and the
transition region is sparse. The profile keeps sharpening as cells are added, so the barrier has
no fixed target and rises with n. A fixed resolution in committor units, with the kernel reflected
at q = 0 and q = 1, estimates a defined quantity at every n.

Three tests, written to experiments/barrier_resolution.json and figures/barrier_resolution.png:

  A  Ground truth. Overdamped diffusion in a double well U(x) = h (x^2 - 1)^2, kT = 1; cells drawn
     from the Boltzmann density, the exact committor between the wells, and the exact h-smoothed
     PMF by quadrature. Both estimators against the truth from 100 to 6,400 cells.
  B  Independent refits. The full pipeline (landscape fit, committor, PMF) on fresh subsamples of
     pancreatic endocrinogenesis and Drosophila spermatogenesis at three cell numbers each.
  C  Resolution. The fixed-resolution barrier against the bandwidth, on the full refits.

    python benchmarks/barrier_resolution.py [--seeds 5]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from thermodynamic_waddington import (  # noqa: E402
    RESOLVED_BANDWIDTH,
    FitConfig,
    committor_free_energy_profile,
    developmental_coordinate,
    fit_landscape,
    resolved_barrier,
)

H = RESOLVED_BANDWIDTH
WELL_HEIGHT = 2.5
SIZES_GT = [100, 200, 400, 800, 1600, 3200, 6400]


# A: ground truth ---------------------------------------------------------------------------------

def double_well(n_grid=40001):
    x = np.linspace(-2.2, 2.2, n_grid)
    u = WELL_HEIGHT * (x * x - 1.0) ** 2
    p = np.exp(-u)
    p /= np.trapezoid(p, x)
    # Committor of overdamped diffusion with constant D between the wells (A: x <= -1, B: x >= 1).
    # q(x) = int_{-1}^{x} e^{U} / int_{-1}^{1} e^{U}, integrated between the wells only (the tails'
    # e^{U} reaches 1e16 and would swamp the interior in floating point).
    a, b = np.searchsorted(x, -1.0), np.searchsorted(x, 1.0)
    eu = np.exp(u[a:b + 1])
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (eu[1:] + eu[:-1]) * np.diff(x[a:b + 1]))])
    q = np.where(x >= 1, 1.0, 0.0)
    q[a:b + 1] = cum / cum[-1]
    return x, p, q


def true_smoothed_barrier(x, p, q, bandwidth, grid=60):
    """The PMF the fixed-resolution estimator targets: E[reflected kernel(xs - q(X))], X ~ p."""
    xs = np.linspace(0.0, 1.0, grid)
    w = p * np.gradient(x)
    dens = np.zeros_like(xs)
    for mirror in (q, -q, 2.0 - q):
        d = (xs[:, None] - mirror[None, :]) / bandwidth
        dens += (np.exp(-0.5 * d * d) * w[None, :]).sum(axis=1)
    dens /= bandwidth * np.sqrt(2 * np.pi)
    f = -np.log(np.clip(dens, 1e-12, None))
    f -= f.min()
    # interior as the estimator defines it, from the population quantiles of q
    cdf = np.cumsum(w) / w.sum()
    lo, hi = q[np.searchsorted(cdf, 0.05)], q[np.searchsorted(cdf, 0.95)]
    inner = (xs >= lo) & (xs <= hi)
    return float(f[inner].max()) if inner.any() else float(f.max())


def ground_truth(replicates, seed=0):
    x, p, q = double_well()
    rng = np.random.default_rng(seed)
    cdf = np.cumsum(p)
    cdf /= cdf[-1]
    out = {"well_height_kt": WELL_HEIGHT, "sizes": SIZES_GT, "silverman": [], "resolved": [], "corrected": [],
           "truth_resolved_kt": round(true_smoothed_barrier(x, p, q, H), 4),
           "atoms": {"q0": round(float(p[q == 0].sum() * (x[1] - x[0])), 3), "q1": round(float(p[q == 1].sum() * (x[1] - x[0])), 3)}}
    for n in SIZES_GT:
        rows = {"silverman": [], "resolved": [], "corrected": []}
        for r in range(replicates):
            qs = np.interp(rng.random(n), cdf, q)
            rows["silverman"].append(committor_free_energy_profile(qs).barrier_kt)
            rb = resolved_barrier(qs, bootstraps=30, seed=r)
            rows["resolved"].append(rb.raw_kt)
            rows["corrected"].append(rb.barrier_kt)
        for k, v in rows.items():
            out[k].append([round(float(np.mean(v)), 4), round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)])
        print(f"  ground truth n={n}: Silverman {np.mean(rows['silverman']):.2f}  resolved {np.mean(rows['resolved']):.2f}"
              f"  corrected {np.mean(rows['corrected']):.2f}  (truth {out['truth_resolved_kt']:.2f})", flush=True)
    return out


# B, C: independent refits on real lineages --------------------------------------------------------

def _fit_q(expr, vel, labels, source, target, seed):
    cfg = FitConfig(neighbors=20, dimensions=6, bootstrap_replicates=2, enable_cycle_decomposition=False,
                    entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2, seed=seed)
    fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=list(labels), metadata={"velocity_status": "derived_proxy"})
    return np.asarray(developmental_coordinate(fit, source, target))


def lineage_data(name):
    if name == "pancreas":
        from pancreas_entropy_validation import BRANCH, load_branch
        expr, vel, labels = load_branch("data/real/endocrinogenesis_day15.h5ad", BRANCH, 10 ** 9, 50, 0)
        return np.asarray(expr), np.asarray(vel), np.asarray(labels), ["Ductal"], ["Beta"], [300, 600, 1200]
    import fly_testis as ft
    d = ft.load()
    return d["spliced"], d["velocity"], d["labels"], ft.SOURCE, ft.TARGET, [500, 1000, 2000]


def refits(name, seeds):
    expr, vel, labels, src, tgt, sizes = lineage_data(name)
    out = {"sizes": sizes, "silverman": [], "resolved": [], "corrected": [], "atoms": [], "sweep": None}
    bandwidths = [0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.14]
    for n in sizes:
        sil, res, cor, atoms, sweeps = [], [], [], [], []
        for s in range(seeds):
            idx = np.random.default_rng(5000 + 97 * n + s).choice(len(labels), size=min(n, len(labels)), replace=False)
            q = _fit_q(expr[idx], vel[idx], labels[idx], src, tgt, s)
            sil.append(committor_free_energy_profile(q).barrier_kt)
            rb = resolved_barrier(q, seed=s)
            res.append(rb.raw_kt)
            cor.append(rb.barrier_kt)
            atoms.append(float(np.mean((q < 1e-6) | (q > 1 - 1e-6))))
            if n == sizes[-1]:
                sweeps.append([committor_free_energy_profile(q, bandwidth=h).barrier_kt for h in bandwidths])
            print(f"  {name} n={n} seed {s}: Silverman {sil[-1]:.2f}  resolved {res[-1]:.2f}  corrected {cor[-1]:.2f}", flush=True)
        for k, v in (("silverman", sil), ("resolved", res), ("corrected", cor)):
            out[k].append({"mean": round(float(np.mean(v)), 4), "sd": round(float(np.std(v, ddof=1)), 4), "values": [round(x, 4) for x in v]})
        out["atoms"].append(round(float(np.mean(atoms)), 3))
    sw = np.asarray(sweeps)
    out["sweep"] = {"bandwidths": bandwidths, "mean": sw.mean(axis=0).round(4).tolist(), "sd": sw.std(axis=0, ddof=1).round(4).tolist()}
    return out


def drift(series):
    return round(series[-1]["mean"] - series[0]["mean"], 3)


# Figure -------------------------------------------------------------------------------------------

def plot(report, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    family = "Inter" if any("Inter" in f.name for f in font_manager.fontManager.ttflist) else "DejaVu Sans"
    plt.rcParams.update({"font.family": family, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 10, "axes.titleweight": "semibold", "axes.titlelocation": "left",
                         "legend.frameon": False, "axes.edgecolor": "#9aa1ab", "xtick.color": "#5f6672", "ytick.color": "#5f6672"})
    SIL, RES, TRUTH, GRID = "#c0392b", "#0f7d99", "#1f2328", "#e3e6eb"
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9), gridspec_kw={"width_ratios": [1.15, 1, 1, 1], "wspace": 0.38})
    fig.subplots_adjust(left=0.05, right=0.99, top=0.82, bottom=0.17)

    ax = axes[0]
    gt = report["ground_truth"]
    n = np.asarray(gt["sizes"])
    for key, color, label in (("silverman", SIL, "Silverman bandwidth"), ("resolved", "#7fb3c4", f"fixed resolution h = {H}"),
                              ("corrected", RES, "fixed resolution, bias-corrected")):
        v = np.asarray(gt[key])
        ax.fill_between(n, v[:, 1], v[:, 2], color=color, alpha=0.15, lw=0)
        ax.plot(n, v[:, 0], "o-", color=color, ms=3.5, lw=1.8, label=label)
    ax.axhline(gt["truth_resolved_kt"], color=TRUTH, ls="--", lw=1, label="exact (h-smoothed)")
    ax.set_xscale("log")
    ax.set_xticks(n)
    ax.set_xticklabels([f"{k:,}" for k in n], fontsize=7)
    ax.minorticks_off()
    ax.set_xlabel("cells")
    ax.set_ylabel("barrier (kT)")
    ax.set_title("Ground truth: double well", pad=8)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.legend(fontsize=7.5, loc="upper left")

    for k, name in enumerate(("pancreas", "fly")):
        ax = axes[1 + k]
        r = report[name]
        n = np.asarray(r["sizes"])
        for key, color in (("silverman", SIL), ("corrected", RES)):
            m = np.asarray([x["mean"] for x in r[key]])
            sd = np.asarray([x["sd"] for x in r[key]])
            ax.errorbar(n, m, yerr=sd, fmt="o-", color=color, ms=4, lw=1.8, capsize=3,
                        label={"silverman": "Silverman", "corrected": "resolved, bias-corrected"}[key])
            for j, x in enumerate(r[key]):
                ax.scatter(np.full(len(x["values"]), n[j]) * (1 + 0.04 * (np.arange(len(x["values"])) - 2) / 2),
                           x["values"], s=8, color=color, alpha=0.35, lw=0)
        ax.set_xscale("log")
        ax.set_xticks(n)
        ax.set_xticklabels([f"{k:,}" for k in n], fontsize=8)
        ax.minorticks_off()
        ax.set_xlabel("cells per independent fit")
        ax.set_title(("Pancreatic endocrinogenesis" if name == "pancreas" else "Drosophila spermatogenesis"), pad=8)
        ax.text(0.98, 0.03, f"drift {n[0]:,} \u2192 {n[-1]:,} cells\nSilverman {r['drift_silverman_kt']:+.2f} kT"
                f"\ncorrected {r['drift_corrected_kt']:+.2f} kT\nsd across fits {np.mean([x['sd'] for x in r['corrected']]):.2f} kT",
                transform=ax.transAxes, va="bottom", ha="right", fontsize=7.2, color="#5f6672",
                bbox={"boxstyle": "round,pad=0.35", "fc": "white", "ec": GRID, "lw": 0.6})
        ax.grid(axis="y", color=GRID, lw=0.6)
        top = max(max(v) for key in ("silverman", "corrected") for x in r[key] for v in [x["values"]])
        ax.set_ylim(0, top * 1.28)
        ax.legend(fontsize=7.5, loc="upper left", ncol=2)

    ax = axes[3]
    for name, color in (("pancreas", "#6a4c93"), ("fly", "#0f7d99")):
        sw = report[name]["sweep"]
        b, m, sd = np.asarray(sw["bandwidths"]), np.asarray(sw["mean"]), np.asarray(sw["sd"])
        ax.fill_between(b, m - sd, m + sd, color=color, alpha=0.15, lw=0)
        ax.plot(b, m, "o-", color=color, ms=3.5, lw=1.8, label=f"{name} ({report[name]['sizes'][-1]:,} cells)")
    ax.axvline(H, color="#1f2328", ls=":", lw=1)
    ax.set_xlabel("resolution h (committor units)")
    ax.set_ylabel("barrier (kT)")
    ax.set_title("Barrier against resolution", pad=8)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.legend(fontsize=7.5)
    for k, ax in enumerate(axes):
        ax.text(-0.14 if k else -0.16, 1.1, "abcd"[k], transform=ax.transAxes, fontsize=13, weight="bold", va="bottom")
    fig.savefig(out, dpi=220)
    plt.close(fig)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--replicates", type=int, default=60)
    ap.add_argument("--out", default="experiments/barrier_resolution.json")
    ap.add_argument("--figure", default="figures/barrier_resolution.png")
    ap.add_argument("--ground-truth-only", action="store_true",
                    help="recompute test A only and reuse the refits already in --out")
    args = ap.parse_args()
    if args.ground_truth_only:
        report = json.loads(Path(args.out).read_text())
    else:
        report = {"benchmark": "barrier_resolution", "resolved_bandwidth": H}
    print("A: ground truth")
    report["ground_truth"] = ground_truth(args.replicates)
    for name in ("pancreas", "fly") if not args.ground_truth_only else ():
        print(f"B: independent refits, {name}")
        r = refits(name, args.seeds)
        r["drift_silverman_kt"] = drift(r["silverman"])
        r["drift_resolved_kt"] = drift(r["resolved"])
        r["drift_corrected_kt"] = drift(r["corrected"])
        report[name] = r
    gt = report["ground_truth"]
    at = gt["sizes"].index(400)
    report["summary"] = {
        "ground_truth_truth_kt": gt["truth_resolved_kt"],
        "ground_truth_error_at_400_kt": {k: round(gt[k][at][0] - gt["truth_resolved_kt"], 3) for k in ("resolved", "corrected")},
        "ground_truth_error_at_largest_kt": {k: round(gt[k][-1][0] - gt["truth_resolved_kt"], 3) for k in ("resolved", "corrected")},
        "ground_truth_silverman_400_to_largest_kt": round(gt["silverman"][-1][0] - gt["silverman"][at][0], 3),
        "pancreas_drift_kt": {k: report["pancreas"][f"drift_{k}_kt"] for k in ("silverman", "resolved", "corrected")},
        "fly_drift_kt": {k: report["fly"][f"drift_{k}_kt"] for k in ("silverman", "resolved", "corrected")},
        "refit_sd_kt": {name: round(float(np.mean([x["sd"] for x in report[name]["corrected"]])), 3) for name in ("pancreas", "fly")},
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(json.dumps(report["summary"], indent=2))
    plot(report, args.figure)


if __name__ == "__main__":
    main()
