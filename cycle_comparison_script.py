"""Compare Schnakenberg cycle EP between Beta and Alpha branches.

Which fundamental cycles in the kNN graph are specifically irreversible
for one fate vs the other? Cycles shared between branches reflect the
common progenitor path; branch-specific cycles may reflect fate-determining
dynamics.
"""
import sys, json, time
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad

from thermodynamic_waddington.preprocessing import velocity_in_transformed_space
from thermodynamic_waddington.graph import build_knn, local_density, estimate_local_diffusion, annotate_edges
from thermodynamic_waddington.entropy_production import EntropyProductionConfig
from thermodynamic_waddington.cycle_decomposition import schnakenberg_decomposition
from thermodynamic_waddington.divergence import estimate_divergence
from run_entropy_replication_v3_dynamical import load_branch_idx, DATA_PATH

adata = ad.read_h5ad(DATA_PATH)
conv  = adata.var_names[~adata.var["fit_alpha"].isna()]
cfg   = EntropyProductionConfig(permutation_replicates=0, bootstrap_replicates=0)

results = {}
all_cycle_eps = {}

for fate, seed in [("Beta", 7), ("Alpha", 1013)]:
    print(f"processing {fate}...", flush=True)
    idx = load_branch_idx(adata, fate, seed=seed)
    sub = adata[idx][:, conv]
    expr = np.asarray(sub.layers["Ms"], dtype=float)
    vel  = np.asarray(sub.layers["velocity"], dtype=float)
    ok   = ~np.isnan(vel).any(axis=0)
    pts, v, _ = velocity_in_transformed_space(expr[:, ok].tolist(), vel[:, ok].tolist())
    graph  = build_knn(pts, 24)
    dens   = local_density(pts, graph, 0.65)
    diff   = estimate_local_diffusion(graph, v, 1e-4)

    t0 = time.time()
    cyc = schnakenberg_decomposition(pts, v, dens, diff, graph, cfg)
    div = estimate_divergence(pts, v, graph, dens, bandwidth=0.65)
    print(f"  {fate}: {len(cyc.cycles)} cycles, global_EP={cyc.global_ep:.4f}, "
          f"cycle_sum={cyc.total_ep_from_cycles:.4f} ({time.time()-t0:.1f}s)")
    print(f"  divergence: mean={div.mean_divergence:.4f}, "
          f"converging={div.fraction_converging:.0%}, "
          f"sign_consistency={div.sign_consistency:.3f}")

    all_cycle_eps[fate] = sorted([c.ep_contribution for c in cyc.cycles], reverse=True)
    results[fate] = {
        "n_cycles": cyc.n_fundamental_cycles,
        "global_ep": cyc.global_ep,
        "cycle_sum": cyc.total_ep_from_cycles,
        "rel_error": cyc.relative_error,
        "n_irreversible": sum(1 for c in cyc.cycles if c.is_irreversible()),
        "top10_ep": [c.ep_contribution for c in cyc.top_cycles[:10]],
        "top10_affinity": [c.affinity for c in cyc.top_cycles[:10]],
        "top10_current": [c.current for c in cyc.top_cycles[:10]],
        "divergence": {
            "mean": div.mean_divergence,
            "weighted_mean": div.weighted_mean,
            "converging_fraction": div.fraction_converging,
            "sign_consistency": div.sign_consistency,
        },
    }

json.dump(results, open("experiments/cycle_comparison_beta_alpha.json", "w"), indent=2)

# ── figure: side by side cycle EP distributions + cumulative ─────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

colors = {"Beta": "#2a5f8f", "Alpha": "#c0622a"}

# top-left: log-scale cycle EP distributions
ax = axes[0, 0]
for fate in ["Beta", "Alpha"]:
    eps = [e for e in all_cycle_eps[fate] if abs(e) > 1e-8]
    ax.hist(eps, bins=50, alpha=0.55, color=colors[fate], label=fate,
            edgecolor="white", linewidth=0.2)
ax.set_xlabel("cycle EP contribution (J_c * A_c)")
ax.set_ylabel("count")
ax.set_title("Cycle EP distribution\nBeta vs Alpha branch")
ax.legend()

# top-right: cumulative EP by cycle rank (how many cycles carry the EP)
ax = axes[0, 1]
for fate in ["Beta", "Alpha"]:
    eps = sorted([abs(e) for e in all_cycle_eps[fate]], reverse=True)
    cumsum = np.cumsum(eps) / (sum(eps) or 1)
    ax.plot(range(1, len(cumsum)+1), cumsum, color=colors[fate], label=fate, lw=2)
ax.axhline(0.5, color="gray", ls="--", lw=1, label="50% of EP")
ax.axhline(0.9, color="gray", ls=":", lw=1, label="90% of EP")
ax.set_xlabel("cycle rank (by |EP contribution|)")
ax.set_ylabel("cumulative fraction of total |cycle EP|")
ax.set_title("How many cycles carry the entropy production?\n(cumulative by rank)")
ax.legend(fontsize=8)
ax.set_xscale("log")

# bottom-left: top-20 cycle EP values
ax = axes[1, 0]
x = np.arange(20)
w = 0.35
for i, fate in enumerate(["Beta", "Alpha"]):
    top20 = [results[fate]["top10_ep"][j] if j < 10 else 0
             for j in range(20)][:len(all_cycle_eps[fate])]
    top20 = [abs(all_cycle_eps[fate][j]) if j < len(all_cycle_eps[fate]) else 0
             for j in range(20)]
    ax.bar(x + i*w, top20, w, color=colors[fate], alpha=0.85, label=fate)
ax.set_xlabel("cycle rank")
ax.set_ylabel("|EP contribution|")
ax.set_title("Top 20 cycles by |EP contribution|")
ax.legend()

# bottom-right: summary comparison table as text
ax = axes[1, 1]
ax.axis("off")
rows = [
    ["Metric", "Beta branch", "Alpha branch"],
    ["Global EP (permutation)", f"{results['Beta']['global_ep']:.4f}", f"{results['Alpha']['global_ep']:.4f}"],
    ["Cycle sum EP", f"{results['Beta']['cycle_sum']:.4f}", f"{results['Alpha']['cycle_sum']:.4f}"],
    ["Relative error", f"{results['Beta']['rel_error']:.1%}", f"{results['Alpha']['rel_error']:.1%}"],
    ["Fundamental cycles", str(results['Beta']['n_cycles']), str(results['Alpha']['n_cycles'])],
    ["Irreversible cycles", str(results['Beta']['n_irreversible']), str(results['Alpha']['n_irreversible'])],
    ["Top-1 cycle EP", f"{results['Beta']['top10_ep'][0]:.5f}", f"{results['Alpha']['top10_ep'][0]:.5f}"],
    ["div(v) mean", f"{results['Beta']['divergence']['mean']:.4f}", f"{results['Alpha']['divergence']['mean']:.4f}"],
    ["% converging cells", f"{results['Beta']['divergence']['converging_fraction']:.0%}", f"{results['Alpha']['divergence']['converging_fraction']:.0%}"],
]
tbl = ax.table(cellText=rows[1:], colLabels=rows[0],
               loc="center", cellLoc="center")
tbl.auto_set_font_size(False)
tbl.set_fontsize(9)
tbl.scale(1, 1.6)
ax.set_title("Beta vs Alpha: Schnakenberg + divergence", pad=12)

plt.suptitle("Schnakenberg cycle decomposition — Beta vs Alpha differentiation branches\n"
             "real pancreatic endocrinogenesis data, 450 cells each, k=24, scVelo dynamical velocity",
             fontsize=10)
plt.tight_layout()
plt.savefig("figures/cycle_comparison_beta_alpha.png", dpi=150)
plt.close()
print("saved figure and results")
