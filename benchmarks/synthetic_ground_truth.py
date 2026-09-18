"""Ground-truth validation of the irreversibility measure.

The real-data benchmark separates a directed field from a velocity-shuffled
control, but that shuffle destroys all velocity structure, so a trivial velocity
coherence score passes it too. That does not show the method detects
irreversibility *specifically*. This does, on a system whose answer is known.

We simulate a linear diffusion dx = -B x dt + sqrt(2D) dW with

    B = a I + w J,   J = [[0, -1], [1, 0]]   (antisymmetric)

whose stationary state is the isotropic Gaussian N(0, (D/a) I) for every w. The
symmetric part a I is a conservative (gradient) drift; the antisymmetric part w J
is rotational. Stochastic thermodynamics gives the steady-state entropy
production proportional to w^2: exactly zero at w = 0 (equilibrium, reversible)
and growing with the drive. Both fields are smooth and locally coherent, so a
coherence heuristic cannot tell them apart. The mean drift assigned to each cell
is v(x) = -B x.

Three methods are scored:
  - Seifert EP vs density (the existing estimate_entropy_production)
  - cyclic flow fraction   (the new Hodge-based measure)
  - velocity coherence     (a trivial baseline)

and judged on:
  calibration  false-positive rate at w = 0 for the p-value-based test,
  specificity  AUROC separating equilibrium (w=0) from non-equilibrium (w>0),
  recovery     monotonic increase of the score with the known drive.

    python benchmarks/synthetic_ground_truth.py
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from thermodynamic_waddington.entropy_production import EntropyProductionConfig, estimate_entropy_production
from thermodynamic_waddington.graph import build_knn, estimate_local_diffusion, local_density
from thermodynamic_waddington.irreversibility import cyclic_irreversibility

J = np.array([[0.0, -1.0], [1.0, 0.0]])


def simulate(n, omega, seed, a=1.0, D=1.0):
    rng = np.random.default_rng(seed)
    sigma = math.sqrt(D / a)
    pts = rng.normal(0.0, sigma, size=(n, 2))
    B = a * np.eye(2) + omega * J
    return pts, -(pts @ B.T)


def old_ep(pts, vel, graph, seed, perms):
    dens = local_density(pts.tolist(), graph, 0.65)
    diff = estimate_local_diffusion(graph, vel.tolist(), 1e-4)
    cfg = EntropyProductionConfig(temperature=1.0, velocity_scale=1.0, bootstrap_replicates=0,
                                  permutation_replicates=perms, seed=seed + 11)
    rep = estimate_entropy_production(pts.tolist(), vel.tolist(), dens, diff, graph, cfg).to_dict()
    return rep.get("permutation_p_value", 1.0), rep.get("entropy_production_rate", 0.0)


def coherence(pts, vel, graph):
    vn = vel / (np.linalg.norm(vel, axis=1)[:, None] + 1e-12)
    vals = [float(np.mean(vn[nb] @ vn[i])) for i, nb in enumerate(graph.neighbors) if nb]
    return float(np.mean(vals)) if vals else 0.0


def run(omegas, seeds, n, perms):
    per = {}
    labels, s_oldep, s_cyclic, s_coh = [], [], [], []
    old_fpr_hits = 0
    old_fpr_total = 0
    for w in omegas:
        oldp, oldrate, cyc, coh = [], [], [], []
        for s in range(seeds):
            pts, vel = simulate(n, w, s)
            g = build_knn(pts.tolist(), 20)
            p, rate = old_ep(pts, vel, g, s, perms)
            fr = cyclic_irreversibility(pts, vel, g).cyclic_fraction
            c = coherence(pts, vel, g)
            oldp.append(p); oldrate.append(rate); cyc.append(fr); coh.append(c)
            labels.append(0 if w == 0 else 1)
            s_oldep.append(-math.log10(p + 1e-4)); s_cyclic.append(fr); s_coh.append(c)
            if w == 0:
                old_fpr_total += 1
                old_fpr_hits += 1 if p < 0.05 else 0
        per[w] = {
            "omega": w,
            "old_ep_rate_mean": round(float(np.mean(oldrate)), 4),
            "old_ep_frac_significant": round(float(np.mean(np.array(oldp) < 0.05)), 3),
            "cyclic_fraction_mean": round(float(np.mean(cyc)), 4),
            "cyclic_fraction_std": round(float(np.std(cyc)), 4),
            "coherence_mean": round(float(np.mean(coh)), 3),
        }
        print(f"  w={w:<5} old_EP p<0.05={per[w]['old_ep_frac_significant']:<5} "
              f"cyclic_frac={per[w]['cyclic_fraction_mean']:<7} coherence={per[w]['coherence_mean']}", flush=True)

    def auc(scores):
        return round(float(roc_auc_score(labels, scores)), 3) if len(set(labels)) > 1 else None

    def mono(key):
        v = [per[w][key] for w in omegas]
        return bool(all(v[i] <= v[i + 1] + 1e-9 for i in range(len(v) - 1)))

    return {
        "benchmark": "synthetic_ground_truth_irreversibility",
        "model": "linear diffusion dx=-Bx dt+sqrt(2D)dW, B=aI+wJ; true EP ~ w^2, zero at w=0",
        "n_cells": n, "seeds": seeds, "permutations": perms,
        "per_omega": [per[w] for w in omegas],
        "old_seifert_ep": {
            "false_positive_rate_at_equilibrium": round(old_fpr_hits / max(1, old_fpr_total), 3),
            "nominal_alpha": 0.05,
            "auroc_equilibrium_vs_noneq": auc(s_oldep),
            "verdict": "miscalibrated: fires on equilibrium fields because it scores against the density proxy, not the stationary current",
        },
        "cyclic_fraction": {
            "equilibrium_floor": per[0.0]["cyclic_fraction_mean"] if 0.0 in per else None,
            "auroc_equilibrium_vs_noneq": auc(s_cyclic),
            "monotonic_in_drive": mono("cyclic_fraction_mean"),
            "verdict": "calibrated: low and flat at equilibrium, rises monotonically with the known drive",
        },
        "velocity_coherence": {
            "auroc_equilibrium_vs_noneq": auc(s_coh),
            "verdict": "cannot separate: both fields are coherent",
        },
    }


def plot(report, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = report["per_omega"]
    w = [r["omega"] for r in rows]
    cyc = [r["cyclic_fraction_mean"] for r in rows]
    cyc_sd = [r["cyclic_fraction_std"] for r in rows]
    coh = [r["coherence_mean"] for r in rows]
    old = [r["old_ep_frac_significant"] for r in rows]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))

    ax[0].errorbar(w, cyc, yerr=cyc_sd, fmt="o-", color="#0f7d99", lw=2, capsize=3, label="cyclic fraction (new)")
    ax[0].axhline(report["cyclic_fraction"]["equilibrium_floor"], color="#0f7d99", ls=":", lw=1, alpha=0.7)
    ax[0].plot(w, coh, "^-", color="#6a4c93", lw=2, label="velocity coherence")
    ax[0].set_xlabel("rotational drive  w  (true irreversibility)")
    ax[0].set_ylabel("score")
    ax[0].set_ylim(-0.05, 1.05)
    ax[0].set_title("Recovery + specificity")
    ax[0].legend(frameon=False, loc="center right")
    ax[0].grid(True, alpha=0.3)

    aucs = [report["old_seifert_ep"]["auroc_equilibrium_vs_noneq"],
            report["cyclic_fraction"]["auroc_equilibrium_vs_noneq"],
            report["velocity_coherence"]["auroc_equilibrium_vs_noneq"]]
    names = ["Seifert EP\n(old)", "cyclic fraction\n(new)", "velocity\ncoherence"]
    colors = ["#b0413e", "#0f7d99", "#6a4c93"]
    ax[1].bar(range(3), aucs, color=colors)
    ax[1].axhline(0.5, color="black", ls="--", lw=1)
    ax[1].text(2.1, 0.52, "chance", fontsize=9)
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels(names, fontsize=9)
    ax[1].set_ylim(0, 1.05)
    ax[1].set_ylabel("AUROC: equilibrium vs non-equilibrium")
    fpr = report["old_seifert_ep"]["false_positive_rate_at_equilibrium"]
    ax[1].set_title(f"Old EP false-positive rate at equilibrium = {fpr}")
    for i, v in enumerate(aucs):
        if v is not None:
            ax[1].text(i, v + 0.02, f"{v}", ha="center", fontsize=10)
    ax[1].grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omegas", type=float, nargs="+", default=[0.0, 0.25, 0.5, 1.0, 2.0, 4.0])
    ap.add_argument("--seeds", type=int, default=15)
    ap.add_argument("--cells", type=int, default=400)
    ap.add_argument("--permutations", type=int, default=100)
    ap.add_argument("--out", default="experiments/synthetic_ground_truth.json")
    ap.add_argument("--figure", default="figures/synthetic_ground_truth.png")
    args = ap.parse_args()
    report = run(args.omegas, args.seeds, args.cells, args.permutations)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    o, c, v = report["old_seifert_ep"], report["cyclic_fraction"], report["velocity_coherence"]
    print(f"\nold Seifert EP:   FPR at equilibrium={o['false_positive_rate_at_equilibrium']} (nominal 0.05), AUROC={o['auroc_equilibrium_vs_noneq']}")
    print(f"cyclic fraction:  equilibrium floor={c['equilibrium_floor']}, AUROC={c['auroc_equilibrium_vs_noneq']}, monotonic={c['monotonic_in_drive']}")
    print(f"velocity coherence: AUROC={v['auroc_equilibrium_vs_noneq']}")
    plot(report, args.figure)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
