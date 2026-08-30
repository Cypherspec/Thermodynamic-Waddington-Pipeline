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
from thermodynamic_waddington.multi_source_propagation import propagate_multi_source
from thermodynamic_waddington.jarzynski import attractor_candidates
from thermodynamic_waddington.mfpt import estimate_mfpt
from thermodynamic_waddington.gradient_alignment import (
    estimate_gradient_alignment, permutation_test_alignment)
from run_entropy_replication_v3_dynamical import load_branch_idx, DATA_PATH

adata = ad.read_h5ad(DATA_PATH)
conv  = adata.var_names[~adata.var["fit_alpha"].isna()]
idx   = load_branch_idx(adata, "Beta", seed=7)
sub   = adata[idx][:, conv]

expr = np.asarray(sub.layers["Ms"],      dtype=float)
vel  = np.asarray(sub.layers["velocity"], dtype=float)
ok   = ~np.isnan(vel).any(axis=0)
expr, vel = expr[:, ok].tolist(), vel[:, ok].tolist()

pts, v, _ = velocity_in_transformed_space(expr, vel)
graph = build_knn(pts, 24)
dens  = local_density(pts, graph, 0.65)
diff  = estimate_local_diffusion(graph, v, 1e-4)
ann   = annotate_edges(pts, v, graph, dens, diff, 1.0, 1.0, 0.0)

ms = propagate_multi_source(ann, len(pts), temperature=1.0, max_paths=64)
energies = ms.energies
attractors = [a for a in attractor_candidates(ann, energies, 0.9) if ms.coverage[a]]
print(f"n attractors (real cells only): {len(attractors)}")

# ── MFPT ──────────────────────────────────────────────────────────────────
# use top-5 attractors by lowest energy to keep runtime manageable
top_atts = sorted(attractors, key=lambda i: energies[i])[:5]
print(f"running MFPT for {len(top_atts)} attractors...", flush=True)
t0 = time.time()
mfpt_result = estimate_mfpt(ann, energies, top_atts, temperature=1.0, max_iter=300)
print(f"  done in {time.time()-t0:.1f}s")
print(f"  converged: {all(mfpt_result.converged)}")
print(f"  mean commitment time: {sum(mfpt_result.commitment_time)/len(mfpt_result.commitment_time):.1f} steps")
commit_times = sorted(t for t in mfpt_result.commitment_time if t < 1e8)
print(f"  commitment time range: [{min(commit_times):.1f}, {max(commit_times):.1f}] steps")

# ── GRADIENT ALIGNMENT ────────────────────────────────────────────────────
print(f"\nrunning gradient alignment + permutation test...", flush=True)
t0 = time.time()
ga = estimate_gradient_alignment(pts, v, energies, graph, bandwidth=0.65)
perm = permutation_test_alignment(pts, v, energies, graph, bandwidth=0.65,
                                   n_permutations=200, seed=42)
print(f"  done in {time.time()-t0:.1f}s")
print(f"  mean alignment: {ga.mean_alignment:.4f}  (>0 = downhill)")
print(f"  fraction downhill: {ga.fraction_downhill:.1%}")
print(f"  weighted mean:  {ga.weighted_mean:.4f}")
print(f"  permutation p-value: {perm['p_value']:.4f}")
print(f"  significant: {perm['significant']}")

# ── FIGURES ───────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# panel A: commitment time distribution
ax = axes[0, 0]
ct = [t for t in mfpt_result.commitment_time if t < 1e8]
ax.hist(ct, bins=35, color="#5a8abf", alpha=0.85, edgecolor="white", lw=0.3)
ax.axvline(sum(ct)/len(ct), color="#c0392b", lw=2,
           label=f"mean={sum(ct)/len(ct):.0f} steps")
ax.set_xlabel("steps to nearest attractor (MFPT)")
ax.set_ylabel("count")
ax.set_title(f"Fate commitment timescale\n"
             f"{len(top_atts)} attractors, 450 cells, Beta branch")
ax.legend(fontsize=8)

# panel B: MFPT landscape (cells colored by commitment time)
ax = axes[0, 1]
emb = np.array(pts)
ct_arr = np.array([t if t < 1e8 else np.nan for t in mfpt_result.commitment_time])
finite_mask = np.isfinite(ct_arr)
sc = ax.scatter(emb[finite_mask, 0], emb[finite_mask, 1],
                c=ct_arr[finite_mask], cmap="YlOrRd", s=22,
                edgecolor="white", lw=0.2, vmax=np.percentile(ct_arr[finite_mask], 95))
for a in top_atts:
    ax.scatter(emb[a, 0], emb[a, 1], marker="*", s=220,
               facecolor="cyan", edgecolor="black", lw=0.8, zorder=5)
plt.colorbar(sc, ax=ax, label="MFPT to nearest attractor (steps)")
ax.set_title("Cells colored by fate commitment time\n(cyan stars = attractors)")
ax.set_xlabel("PCA dim 1"); ax.set_ylabel("PCA dim 2")

# panel C: gradient alignment distribution
ax = axes[1, 0]
ax.hist(ga.per_cell, bins=35, color="#8abf5a", alpha=0.85, edgecolor="white", lw=0.3)
ax.axvline(0, color="black", lw=1, ls="--", label="no alignment")
ax.axvline(ga.mean_alignment, color="#c0392b", lw=2,
           label=f"mean={ga.mean_alignment:.3f}")
ax.set_xlabel("alignment: cos(velocity, -grad F)")
ax.set_ylabel("count")
ax.set_title(f"Velocity-landscape gradient alignment\n"
             f"{ga.fraction_downhill:.0%} of cells move downhill  "
             f"(p={perm['p_value']:.3f})")
ax.legend(fontsize=8)

# panel D: cells colored by alignment
ax = axes[1, 1]
al_arr = np.array(ga.per_cell)
sc = ax.scatter(emb[:, 0], emb[:, 1], c=al_arr, cmap="RdYlGn",
                s=22, edgecolor="white", lw=0.2, vmin=-1, vmax=1)
plt.colorbar(sc, ax=ax, label="alignment (green=downhill, red=uphill)")
ax.set_title("Per-cell gradient alignment\n(green = velocity points downhill on landscape)")
ax.set_xlabel("PCA dim 1"); ax.set_ylabel("PCA dim 2")

plt.suptitle("MFPT + Gradient Alignment — real pancreatic Beta branch\n"
             "scVelo dynamical velocity, multi-source free-energy landscape", fontsize=10)
plt.tight_layout()
plt.savefig("figures/mfpt_gradient_alignment.png", dpi=150)
plt.close()
print("\nfigure saved")

results = {
    "mfpt": {
        "n_attractors_used": len(top_atts),
        "attractor_indices": top_atts,
        "converged_all": all(mfpt_result.converged),
        "mean_commitment_time": sum(mfpt_result.commitment_time)/len(mfpt_result.commitment_time),
        "min_commitment_time": min(commit_times),
        "max_commitment_time": max(commit_times),
        **mfpt_result.to_dict(),
    },
    "gradient_alignment": {
        **ga.to_dict(),
        "permutation_test": perm,
    },
}
json.dump(results, open("experiments/mfpt_gradient_pancreas_beta.json","w"), indent=2)
print("results saved")
