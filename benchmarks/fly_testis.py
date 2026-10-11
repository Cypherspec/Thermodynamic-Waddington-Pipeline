"""Drosophila spermatogenesis: the pipeline on a fourth tissue and a second phylum.

Runs the calibrated thermodynamic readouts on the male germline of the Fly Cell
Atlas testis (Li et al., Science 2022), with spliced/unspliced counts built
from the raw reads by benchmarks/fly_testis_prepare.py. The lineage runs from
spermatogonia through the spermatocyte program to elongating spermatids, with
the Fly Cell Atlas annotation as stage labels.

The setup matches the three vertebrate tissues (irreversibility_real.py,
commitment_barrier.py) so the numbers are directly comparable:

  1. velocity check   unspliced fraction by stage. Transcription largely stops
                      after meiosis in Drosophila (spermatid mRNAs are made in
                      spermatocytes and stored), so a faithful quantification
                      must show the nascent fraction collapse in spermatids.
  2. irreversibility  Hodge cyclic fraction vs a matched gradient floor and a
                      shuffled ceiling, 1500 cells x 100 genes, 5 seeds.
  3. where            per-stage cyclic fraction (magnitude-invariant) against a
                      per-stage gradient floor on the same cells, plus the
                      cyclic-energy share (1 = uniform), which also tracks speed.
  4. commitment       committor spermatogonia -> late elongating spermatids:
                      where it crosses 0.5, the PMF barrier in kT with a 95%
                      CI, and stage-ordering Spearman against PC1.
  5. drivers          FlyBase-linked genes that switch hardest across the
                      commitment crossing (q = 0.5), plus the ranks of a
                      regulator list fixed before the analysis was run.

    python benchmarks/fly_testis.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from thermodynamic_waddington import (
    FitConfig,
    commitment_profile,
    committor_free_energy_profile,
    developmental_coordinate,
    fit_landscape,
)
from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.irreversibility import (
    cyclic_flow_per_cell,
    cyclic_fraction_per_cell,
    cyclic_irreversibility,
)

PATH = "data/fly/testis_velocity.h5ad"

# Fly Cell Atlas germline annotations in developmental order (stage index).
# The unnumbered "spermatocyte" cluster stays in the graph but is not scored,
# since its position inside the spermatocyte program is not resolved.
STAGES = {
    "spermatogonium": 0,
    "mid-late proliferating spermatogonia": 1,
    "spermatogonium-spermatocyte transition": 2,
    "spermatocyte 0": 3, "spermatocyte 1": 4, "spermatocyte 2": 5, "spermatocyte 3": 6,
    "spermatocyte 4": 7, "spermatocyte 5": 8, "spermatocyte 6": 9, "spermatocyte 7a": 10,
    "late primary spermatocyte": 11,
    "spermatid": 12,
    "early elongation stage spermatid": 13,
    "early-mid elongation-stage spermatid": 14,
    "mid-late elongation-stage spermatid": 15,
}
UNSCORED = ["spermatocyte"]
SOURCE = ["spermatogonium", "mid-late proliferating spermatogonia"]
TARGET = ["mid-late elongation-stage spermatid"]
PHASES = [("spermatogonia", 0, 2), ("spermatocytes", 3, 11), ("spermatids", 12, 15)]


def _dense(m):
    return np.asarray(m.toarray() if hasattr(m, "toarray") else m, dtype=float)


def load(path=PATH, n_genes=100):
    import anndata as ad

    a = ad.read_h5ad(path)
    lab = a.obs["annotation"].astype(str).to_numpy()
    keep = np.isin(lab, list(STAGES) + UNSCORED)
    a, lab = a[keep].copy(), lab[keep]
    s_all, u_all = a.layers["spliced"], a.layers["unspliced"]
    s_tot = np.asarray(s_all.sum(axis=1)).ravel()
    u_tot = np.asarray(u_all.sum(axis=1)).ravel()
    # top genes by spliced variance, the same rule as the vertebrate benchmarks
    mean = np.asarray(s_all.mean(axis=0)).ravel()
    mean_sq = np.asarray(s_all.multiply(s_all).mean(axis=0)).ravel() if hasattr(s_all, "multiply") else (s_all ** 2).mean(axis=0)
    top = np.argsort(mean_sq - mean * mean)[::-1][:n_genes]
    sp, un = _dense(s_all[:, top]), _dense(u_all[:, top])
    stage = np.array([STAGES.get(x, -1) for x in lab])
    return {"rows": np.where(keep)[0], "spliced": sp, "velocity": un - sp, "labels": lab, "stage": stage,
            "unspliced_fraction": u_tot / np.maximum(s_tot + u_tot, 1), "umi": s_tot + u_tot,
            "umap": np.asarray(a.obsm["X_umap_fca"]), "genes": a.var_names[top].tolist()}


def pca_project(expr, vel, dims=6):
    Xc = expr - expr.mean(axis=0)
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    comp = vt[:dims].T
    return Xc @ comp, vel @ comp


def velocity_check(d):
    rows = []
    for name, lo, hi in PHASES:
        m = (d["stage"] >= lo) & (d["stage"] <= hi)
        rows.append({"phase": name, "n_cells": int(m.sum()),
                     "unspliced_fraction": round(float(np.median(d["unspliced_fraction"][m])), 4),
                     "median_umi": int(np.median(d["umi"][m]))})
    by_stage = []
    for lab, st in sorted(STAGES.items(), key=lambda kv: kv[1]):
        m = d["labels"] == lab
        if m.sum():
            by_stage.append({"stage": st, "label": lab, "n_cells": int(m.sum()),
                             "unspliced_fraction": round(float(np.median(d["unspliced_fraction"][m])), 4)})
    sc = next(r["unspliced_fraction"] for r in rows if r["phase"] == "spermatocytes")
    st = next(r["unspliced_fraction"] for r in rows if r["phase"] == "spermatids")
    return {"by_phase": rows, "by_stage": by_stage,
            "spermatid_over_spermatocyte": round(st / max(sc, 1e-9), 3)}


def irreversibility(d, cells, seeds):
    n = len(d["labels"])
    real, floor, ceil = [], [], []
    stage_density = {lab: [] for lab in STAGES}
    stage_frac = {lab: [] for lab in STAGES}
    stage_floor = {lab: [] for lab in STAGES}
    for s in range(seeds):
        idx = np.random.default_rng(s).choice(n, size=min(cells, n), replace=False)
        pts, v = pca_project(d["spliced"][idx], d["velocity"][idx])
        g = build_knn(pts.tolist(), 20)
        real.append(cyclic_irreversibility(pts, v, g, seed=s).cyclic_fraction)
        floor.append(cyclic_irreversibility(pts, pts.mean(axis=0) - pts, g, seed=s).cyclic_fraction)
        shuf = v[np.random.default_rng(100 + s).permutation(len(idx))]
        ceil.append(cyclic_irreversibility(pts, shuf, g, seed=s).cyclic_fraction)
        share = np.asarray(cyclic_flow_per_cell(pts, v, g)) * len(idx)   # 1 = uniform
        local = cyclic_fraction_per_cell(pts, v, g)
        local_floor = cyclic_fraction_per_cell(pts, pts.mean(axis=0) - pts, g)
        lab = d["labels"][idx]
        for name in STAGES:
            m = lab == name
            if m.sum() >= 5:
                stage_density[name].append(float(share[m].mean()))
                stage_frac[name].append(float(local[m].mean()))
                stage_floor[name].append(float(local_floor[m].mean()))
    where = [{"label": k, "stage": STAGES[k],
              "cyclic_fraction": round(float(np.mean(stage_frac[k])), 3),
              "cyclic_fraction_std": round(float(np.std(stage_frac[k])), 3),
              "gradient_floor": round(float(np.mean(stage_floor[k])), 3),
              "excess_over_floor": round(float(np.mean(stage_frac[k]) - np.mean(stage_floor[k])), 3),
              "relative_cyclic_energy": round(float(np.mean(v)), 3),
              "relative_cyclic_energy_std": round(float(np.std(v)), 3), "seeds_present": len(v)}
             for k, v in stage_density.items() if v]
    return {
        "n_cells": int(min(cells, n)), "seeds": seeds,
        "cyclic_fraction_real": round(float(np.mean(real)), 4),
        "cyclic_fraction_real_std": round(float(np.std(real)), 4),
        "gradient_floor": round(float(np.mean(floor)), 4),
        "random_ceiling": round(float(np.mean(ceil)), 4),
        "excess_over_floor": round(float(np.mean(real) - np.mean(floor)), 4),
        "per_stage": sorted(where, key=lambda r: r["stage"]),
    }


def _spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


def commitment(d, cells, seeds, grid=60):
    n = len(d["labels"])
    barriers, bqs, labels_at, rho_q, rho_pc1, curves = [], [], [], [], [], []
    for s in range(seeds):
        idx = np.random.default_rng(1000 + s).choice(n, size=min(cells, n), replace=False)
        expr, vel, lab, stage = d["spliced"][idx], d["velocity"][idx], d["labels"][idx].tolist(), d["stage"][idx]
        cfg = FitConfig(neighbors=20, dimensions=6, bootstrap_replicates=2, enable_cycle_decomposition=False,
                        entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2, seed=s)
        fit = fit_landscape(expr.tolist(), vel.tolist(), config=cfg, labels=lab,
                            metadata={"velocity_status": "derived_proxy"})
        q = np.asarray(developmental_coordinate(fit, SOURCE, TARGET))
        prof = committor_free_energy_profile(q.tolist(), temperature=1.0, grid=grid)
        rep = commitment_profile(fit, SOURCE, TARGET)
        scored = stage >= 0
        pc1 = pca_project(expr, vel, dims=1)[0][:, 0]
        barriers.append(prof.barrier_kt)
        bqs.append(prof.barrier_q)
        curves.append(list(prof.free_energy_kt))
        labels_at.append(rep.commitment_label)
        rho_q.append(_spearman(q[scored], stage[scored]))
        rho_pc1.append(abs(_spearman(pc1[scored], stage[scored])))
        print(f"  seed {s}: barrier {prof.barrier_kt:.2f} kT at q={prof.barrier_q:.2f}, commits at "
              f"'{rep.commitment_label}', Spearman committor {rho_q[-1]:.3f} vs PC1 {rho_pc1[-1]:.3f}", flush=True)
    b = np.asarray(barriers)
    sem = float(b.std(ddof=1) / np.sqrt(len(b))) if len(b) > 1 else 0.0
    commit_counts = {lab: labels_at.count(lab) for lab in sorted(set(labels_at))}
    return {
        "n_cells": int(min(cells, n)), "seeds": seeds,
        "source": SOURCE, "target": TARGET,
        "barrier_kt_mean": round(float(b.mean()), 3),
        "barrier_kt_ci95": [round(float(b.mean() - 1.96 * sem), 3), round(float(b.mean() + 1.96 * sem), 3)],
        "barrier_q_mean": round(float(np.mean(bqs)), 3),
        "commitment_label_counts": commit_counts,
        "spearman_committor_mean": round(float(np.mean(rho_q)), 3),
        "spearman_committor_std": round(float(np.std(rho_q)), 3),
        "spearman_pc1_mean": round(float(np.mean(rho_pc1)), 3),
        "committor_wins": int(sum(a > b for a, b in zip(rho_q, rho_pc1))),
        "grid": np.linspace(0.0, 1.0, grid).tolist(),
        "curves": curves,
    }


# Fixed before the analysis was run (not chosen from its output): canonical
# regulators of the spermatocyte program and its exits, from the literature.
KNOWN_REGULATORS = {
    "bam": "spermatogonia -> spermatocyte switch", "bgcn": "spermatogonia -> spermatocyte switch",
    "aly": "meiotic arrest (tMAC)", "comr": "meiotic arrest (tMAC)", "tomb": "meiotic arrest (tMAC)",
    "topi": "meiotic arrest (tMAC)", "achi": "meiotic arrest (tMAC)", "vis": "meiotic arrest (tMAC)",
    "can": "meiotic arrest (tTAF)", "mia": "meiotic arrest (tTAF)", "sa": "meiotic arrest (tTAF)",
    "nht": "meiotic arrest (tTAF)", "rye": "meiotic arrest (tTAF)",
    "fzo": "spermatid mitochondria", "Mst87F": "spermatid tail", "ProtA": "protamine", "ProtB": "protamine",
    "dj": "spermatid tail",
}


def _gene_symbols(t2g="data/fly/ref/t2g.txt"):
    names = {}
    for line in Path(t2g).read_text().splitlines():
        f = line.split("\t")
        if len(f) > 2 and f[1] not in names:
            names[f[1]] = f[2] or f[1]
    return names


def _switch_scores(x, q, window):
    before = (q >= 0.5 - window) & (q < 0.5)
    after = (q >= 0.5) & (q <= 0.5 + window)
    sd = np.sqrt(np.maximum(np.asarray(x.multiply(x).mean(axis=0)).ravel() - np.asarray(x.mean(axis=0)).ravel() ** 2, 1e-12))
    mu_b = np.asarray(x[before].mean(axis=0)).ravel()
    mu_a = np.asarray(x[after].mean(axis=0)).ravel()
    return (mu_a - mu_b) / sd, int(before.sum()), int(after.sum())


def drivers(d, maps, path=PATH, top=15, window=0.2, min_frac=0.05, top_k=50):
    """Genes that switch reproducibly as cells cross the commitment point (q = 0.5).

    Each independent committor fit (a random cell subsample) gives a switch score
    per gene: the change in mean log expression across q = 0.5, in SD units. A
    gene is reported only if it switches the same way in every fit; the overlap
    of top-k lists between fits measures how reproducible the ranking is.
    """
    import anndata as ad

    a = ad.read_h5ad(path)
    s_all = a.layers["spliced"][d["rows"]]
    detected = np.asarray((s_all > 0).mean(axis=0)).ravel()
    keep = np.where(detected >= min_frac)[0]
    tot = np.asarray(s_all.sum(axis=1)).ravel()
    x_all = s_all[:, keep].multiply(1e4 / np.maximum(tot, 1)[:, None]).tocsr().log1p()
    scores, counts = [], []
    for idx, q in maps:
        sc, nb, na = _switch_scores(x_all[idx], q, window)
        scores.append(sc)
        counts.append((nb, na))
    S = np.vstack(scores)
    med = np.median(S, axis=0)
    consistent = (np.sign(S) == np.sign(med)).all(axis=0) & (med != 0)
    tops = [set(np.argsort(-np.abs(row))[:top_k]) for row in S]
    overlaps = [len(tops[i] & tops[j]) / top_k for i in range(len(tops)) for j in range(i + 1, len(tops))]
    names = _gene_symbols()
    fbgn = a.var_names[keep].to_numpy()
    order = np.argsort(-np.abs(med))
    rank = {fbgn[i]: r + 1 for r, i in enumerate(order)}

    def row(i):
        return {"fbgn": str(fbgn[i]), "symbol": names.get(fbgn[i], fbgn[i]),
                "switch_score": round(float(med[i]), 3),
                "score_range": [round(float(S[:, i].min()), 3), round(float(S[:, i].max()), 3)],
                "flybase": f"https://flybase.org/reports/{fbgn[i]}"}

    up = [row(i) for i in np.argsort(-med) if consistent[i] and med[i] > 0][:top]
    down = [row(i) for i in np.argsort(med) if consistent[i] and med[i] < 0][:top]
    by_symbol = {names.get(g, g): k for k, g in enumerate(fbgn)}
    known = []
    for sym, role in KNOWN_REGULATORS.items():
        k = by_symbol.get(sym)
        known.append({"symbol": sym, "role": role, "fbgn": str(fbgn[k]) if k is not None else None,
                      "rank": rank.get(fbgn[k]) if k is not None else None,
                      "switch_score": round(float(med[k]), 3) if k is not None else None})
    ranked = [k["rank"] for k in known if k["rank"]]
    hits = int(sum(r <= 0.05 * len(keep) for r in ranked))
    from math import comb
    p_enrich = sum(comb(len(ranked), k) * 0.05 ** k * 0.95 ** (len(ranked) - k) for k in range(hits, len(ranked) + 1))
    return {
        "n_genes_tested": int(len(keep)), "n_fits": len(maps), "window": window,
        "cells_before_after_per_fit": counts,
        "top_k_overlap_between_fits_mean": round(float(np.mean(overlaps)), 3) if overlaps else None,
        "top_k_overlap_chance": round(top_k / len(keep), 4), "top_k": top_k,
        "n_sign_consistent_genes": int(consistent.sum()),
        "top_up": up, "top_down": down, "known_regulators": known,
        "known_in_top_5pct": hits, "known_detected": len(ranked),
        "known_expected_by_chance": round(0.05 * len(ranked), 2), "known_enrichment_binomial_p": round(float(p_enrich), 4),
        "note": "switch score = median over independent committor fits of the change in mean log expression "
                "across q = 0.5 (SD units); reported genes switch the same way in every fit. A ranking for "
                "follow-up, not a causal test.",
    }


def committor_map(d, cells=3000, seed=0):
    """One larger fit for the figure: committor on the FCA UMAP."""
    n = len(d["labels"])
    idx = np.random.default_rng(seed).choice(n, size=min(cells, n), replace=False)
    cfg = FitConfig(neighbors=20, dimensions=6, bootstrap_replicates=2, enable_cycle_decomposition=False,
                    entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2, seed=seed)
    fit = fit_landscape(d["spliced"][idx].tolist(), d["velocity"][idx].tolist(), config=cfg,
                        labels=d["labels"][idx].tolist(), metadata={"velocity_status": "derived_proxy"})
    return idx, np.asarray(developmental_coordinate(fit, SOURCE, TARGET))


def plot(report, d, idx, q, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    teal, grey, red, ink = "#0f7d99", "#9aa5b1", "#b0413e", "#1f2933"
    fig, axs = plt.subplots(2, 3, figsize=(19, 9.5))
    ax = axs[:, :2]

    a = ax[0, 0]
    um = d["umap"]
    a.scatter(um[:, 0], um[:, 1], s=2, c="#e4e7eb", lw=0)
    sc = a.scatter(um[idx, 0], um[idx, 1], s=5, c=q, cmap="viridis", lw=0)
    fig.colorbar(sc, ax=a, label="committor q")
    a.set_title("A  Committor on the Fly Cell Atlas UMAP\n(0 = spermatogonia, 1 = late elongating spermatids)")
    a.set_xticks([]); a.set_yticks([])

    a = ax[0, 1]
    c = report["commitment"]
    grid, curves = np.asarray(c["grid"]), np.asarray(c["curves"])
    for cv in curves:
        a.plot(grid, cv, color=teal, alpha=0.18, lw=1)
    a.plot(grid, curves.mean(axis=0), color=teal, lw=2.6)
    a.axvline(c["barrier_q_mean"], color=red, ls="--", lw=1)
    ci = c["barrier_kt_ci95"]
    a.set_title(f"B  Commitment barrier {c['barrier_kt_mean']:.1f} kT [{ci[0]:.1f}, {ci[1]:.1f}] ({c['seeds']} seeds)")
    a.set_xlabel("committor q"); a.set_ylabel("free energy (kT)"); a.grid(alpha=0.3)

    a = ax[1, 0]
    rows = report["vertebrate_reference"] + [{"dataset": "fly testis germline", **{k: report["irreversibility"][k] for k in ("cyclic_fraction_real", "gradient_floor", "random_ceiling")}}]
    x = np.arange(len(rows))
    a.bar(x - 0.27, [r["cyclic_fraction_real"] for r in rows], 0.27, color=teal, label="real velocity")
    a.bar(x, [r["gradient_floor"] for r in rows], 0.27, color=grey, label="reversible floor")
    a.bar(x + 0.27, [r["random_ceiling"] for r in rows], 0.27, color="#d9dde2", label="shuffled ceiling")
    a.set_xticks(x, [r["dataset"].replace(" ", "\n", 1) for r in rows], fontsize=9)
    a.set_ylabel("cyclic (irreversible) fraction")
    a.set_title("C  Calibrated irreversibility: fly vs three vertebrate tissues")
    a.legend(frameon=False, fontsize=9); a.grid(axis="y", alpha=0.3)

    a = ax[1, 1]
    vs = report["velocity_check"]["by_stage"]
    st = [r["stage"] for r in vs]
    a.plot(st, [r["unspliced_fraction"] for r in vs], "o-", color=ink, lw=2, label="unspliced fraction (median)")
    a.axvspan(11.5, 15.5, color=red, alpha=0.08)
    a.text(13.5, a.get_ylim()[1] * 0.92 if a.get_ylim()[1] else 0.1, "post-meiotic", ha="center", color=red, fontsize=9)
    a.set_xlabel("germline stage (0 = spermatogonium ... 15 = mid-late elongating spermatid)")
    a.set_ylabel("unspliced / total UMIs")
    b2 = a.twinx()
    w = report["irreversibility"]["per_stage"]
    b2.errorbar([r["stage"] for r in w], [r["cyclic_fraction"] for r in w], yerr=[r["cyclic_fraction_std"] for r in w],
                fmt="s-", color=teal, lw=1.5, capsize=2, label="local cyclic fraction")
    b2.plot([r["stage"] for r in w], [r["gradient_floor"] for r in w], ":", color=grey, lw=1.5,
            label="reversible floor (same cells)")
    b2.set_ylabel("local cyclic fraction", color=teal)
    a.set_title("D  Velocity check (post-meiotic shutdown) and local\ncyclic fraction vs its reversible floor, by stage")
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3)

    for a, key, color, title in ((axs[0, 2], "top_up", teal, "E  Switch on at commitment"),
                                 (axs[1, 2], "top_down", red, "F  Switch off at commitment")):
        genes = report["drivers"][key][:12][::-1]
        a.barh([g["symbol"] for g in genes], [g["switch_score"] for g in genes], color=color)
        a.set_xlabel("switch score across q = 0.5 (SD units)")
        a.set_title(f"{title}\n(FlyBase genes, same sign in all {report['drivers']['n_fits']} fits; candidates, not causal)")
        a.grid(axis="x", alpha=0.3)
    fig.suptitle("Thermodynamic Waddington on Drosophila spermatogenesis (Fly Cell Atlas)", fontsize=14, weight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--path", default=PATH)
    ap.add_argument("--cells", type=int, default=1500)
    ap.add_argument("--genes", type=int, default=100)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--commit-cells", type=int, default=600)
    ap.add_argument("--commit-seeds", type=int, default=8)
    ap.add_argument("--driver-fits", type=int, default=5)
    ap.add_argument("--driver-cells", type=int, default=3000)
    args = ap.parse_args()

    d = load(args.path, args.genes)
    print(f"germline cells: {len(d['labels'])}")
    print("1. velocity check (unspliced fraction by phase)")
    vc = velocity_check(d)
    for r in vc["by_phase"]:
        print(f"  {r['phase']:14s} n={r['n_cells']:5d}  unspliced={r['unspliced_fraction']:.3f}  median UMI={r['median_umi']}")
    print("2-3. calibrated irreversibility")
    irr = irreversibility(d, args.cells, args.seeds)
    print(f"  real={irr['cyclic_fraction_real']:.3f}  floor={irr['gradient_floor']:.3f}  "
          f"ceiling={irr['random_ceiling']:.3f}  excess={irr['excess_over_floor']:+.3f}")
    depth = {r["label"]: r for r in vc["by_stage"]}
    umi = {lab: float(np.median(d["umi"][d["labels"] == lab])) for lab in depth}
    ps = [r for r in irr["per_stage"] if r["label"] in umi]
    irr["excess_vs_depth_spearman"] = round(_spearman([r["excess_over_floor"] for r in ps], [umi[r["label"]] for r in ps]), 3)
    irr["excess_vs_unspliced_spearman"] = round(_spearman([r["excess_over_floor"] for r in ps],
                                                          [depth[r["label"]]["unspliced_fraction"] for r in ps]), 3)
    irr["energy_vs_depth_spearman"] = round(_spearman([r["relative_cyclic_energy"] for r in ps], [umi[r["label"]] for r in ps]), 3)
    pre = [r["excess_over_floor"] for r in ps if r["stage"] <= 11]
    post = [r["excess_over_floor"] for r in ps if r["stage"] >= 12]
    irr["excess_pre_meiotic_mean"], irr["excess_post_meiotic_mean"] = round(float(np.mean(pre)), 3), round(float(np.mean(post)), 3)
    print(f"  per stage: excess pre-meiotic {irr['excess_pre_meiotic_mean']:+.3f}, post-meiotic {irr['excess_post_meiotic_mean']:+.3f}")
    print(f"  Spearman over stages: excess vs depth {irr['excess_vs_depth_spearman']}, excess vs unspliced "
          f"{irr['excess_vs_unspliced_spearman']}; cyclic energy vs depth {irr['energy_vs_depth_spearman']}")
    print("4. commitment")
    com = commitment(d, args.commit_cells, args.commit_seeds)
    ref_path = Path("experiments/irreversibility_real.json")
    vertebrate = json.loads(ref_path.read_text())["results"] if ref_path.exists() else []
    report = {
        "benchmark": "fly_testis_germline",
        "dataset": "Fly Cell Atlas testis, library FCA59 (E-MTAB-10519), spliced/unspliced via kallisto|bustools nac",
        "n_germline_cells": int(len(d["labels"])),
        "velocity": "unspliced-minus-spliced proxy, top 100 genes by spliced variance (as for the vertebrate tissues)",
        "velocity_check": vc,
        "irreversibility": irr,
        "commitment": com,
        "vertebrate_reference": [{k: r[k] for k in ("dataset", "cyclic_fraction_real", "gradient_floor", "random_ceiling")} for r in vertebrate],
    }
    print("5. commitment drivers (FlyBase), over independent committor fits")
    maps = [committor_map(d, cells=args.driver_cells, seed=k) for k in range(args.driver_fits)]
    idx, q = maps[0]
    report["drivers"] = drivers(d, maps)
    dr = report["drivers"]
    print("  up:  ", ", ".join(g["symbol"] for g in dr["top_up"][:10]))
    print("  down:", ", ".join(g["symbol"] for g in dr["top_down"][:10]))
    print(f"  top-{dr['top_k']} overlap between fits: {dr['top_k_overlap_between_fits_mean']} "
          f"(chance {dr['top_k_overlap_chance']}); sign-consistent genes: {dr['n_sign_consistent_genes']}")
    print(f"  known regulators in top 5%: {dr['known_in_top_5pct']}/{dr['known_detected']} detected "
          f"(chance {dr['known_expected_by_chance']}, binomial p={dr['known_enrichment_binomial_p']}; "
          f"{dr['n_genes_tested']} genes tested)")
    for k in dr["known_regulators"]:
        print(f"    {k['symbol']:7s} rank {k['rank']}  score {k['switch_score']}  ({k['role']})")
    Path("experiments").mkdir(exist_ok=True)
    Path("experiments/fly_testis.json").write_text(json.dumps(report, indent=2))
    Path("figures").mkdir(exist_ok=True)
    plot(report, d, idx, q, "figures/fly_testis.png")
    print("wrote experiments/fly_testis.json and figures/fly_testis.png")


if __name__ == "__main__":
    main()
