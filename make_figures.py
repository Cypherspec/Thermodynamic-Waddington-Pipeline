"""Real figures generated only from this session's own verified experiment
outputs (v2/v3 pancreas replication, kNN sensitivity, dentate gyrus
independent-dataset tests). No numbers invented, nothing borrowed from
the older experiments/ files this session did not produce or verify."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def load(name):
    return json.load(open(f"experiments/{name}"))

# ---------- Figure 1: forest plot of every significance test run this session ----------
tests = [
    ("Pancreas Beta (proxy velocity, n=450)", 0.0025, True),
    ("Pancreas Alpha (proxy velocity, n=450)", 0.0050, True),
    ("Pancreas Beta (dynamical velocity, n=450)", 0.0025, True),
    ("Pancreas Alpha (dynamical velocity, n=450)", 0.0050, True),
    ("Dentate: 5-stage, n=75 (nIPC-bottlenecked)", 0.9526, False),
    ("Dentate: 4-stage, n=180 (nIPC dropped)", 0.9975, False),
    ("Dentate: P12 only, n=120", 0.6977, False),
    ("Dentate: P35 only, n=56", 0.4718, False),
    ("Dentate: RadialGlia->GABA, n=90", 0.9900, False),
    ("Dentate: RadialGlia->Granule (exploratory), n=90", 0.0465, False),
    ("Dentate: RadialGlia->Granule (confirmatory, pre-specified), n=90", 0.0270, False),
]
labels = [t[0] for t in tests]
pvals = [t[1] for t in tests]
passed = [t[2] for t in tests]
colors = ["#2a7f3f" if p else "#b3402f" for p in passed]

fig, ax = plt.subplots(figsize=(10, 6))
y = np.arange(len(tests))[::-1]
ax.barh(y, [-np.log10(max(p, 1e-4)) for p in pvals], color=colors, height=0.6)
ax.axvline(-np.log10(0.05), color="black", linestyle="--", linewidth=1, label="p = 0.05")
ax.axvline(-np.log10(0.01), color="gray", linestyle=":", linewidth=1, label="p = 0.01 (confirmatory bar)")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel("-log10(permutation p-value)")
ax.set_title("Every entropy-production significance test run this session\n(green = significant at p<0.05, red = not)")
ax.legend(loc="lower right", fontsize=8)
plt.tight_layout()
plt.savefig("figures/all_tests_forest_plot.png", dpi=150)
plt.close()
print("figure 1 done")

# ---------- Figure 2: kNN sensitivity sweep (real data) ----------
sens = load("knn_sensitivity_dynamical.json")
ks = [r["k"] for r in sens["Beta"]]
beta_ep = [r["entropy_production_rate"] for r in sens["Beta"]]
beta_p = [r["permutation_p_value"] for r in sens["Beta"]]
alpha_ep = [r["entropy_production_rate"] for r in sens["Alpha"]]
alpha_p = [r["permutation_p_value"] for r in sens["Alpha"]]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
ax1.plot(ks, beta_ep, "o-", color="#2a5f8f", label="Beta")
ax1.plot(ks, alpha_ep, "s-", color="#c0622a", label="Alpha")
ax1.set_xlabel("k (kNN neighbors)")
ax1.set_ylabel("entropy production rate\n(raw, not comparable across k -- see note)")
ax1.set_title("Magnitude vs k\n(climbs mechanically with k, expected)")
ax1.legend()

ax2.plot(ks, beta_p, "o-", color="#2a5f8f", label="Beta")
ax2.plot(ks, alpha_p, "s-", color="#c0622a", label="Alpha")
ax2.axhline(0.05, color="black", linestyle="--", linewidth=1, label="p = 0.05")
ax2.set_xlabel("k (kNN neighbors)")
ax2.set_ylabel("permutation p-value")
ax2.set_title("Significance vs k\nBeta robust throughout; Alpha only in k=18-24")
ax2.legend()
plt.tight_layout()
plt.savefig("figures/knn_sensitivity.png", dpi=150)
plt.close()
print("figure 2 done")

# ---------- Figure 3: stationarity residual vs k, showing best-fit region ----------
beta_res = [r["stationarity_residual"] for r in sens["Beta"]]
alpha_res = [r["stationarity_residual"] for r in sens["Alpha"]]
fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.plot(ks, beta_res, "o-", color="#2a5f8f", label="Beta")
ax.plot(ks, alpha_res, "s-", color="#c0622a", label="Alpha")
ax.axvspan(18, 24, color="gray", alpha=0.15, label="region of best model fit")
ax.set_xlabel("k (kNN neighbors)")
ax.set_ylabel("stationarity residual (lower = better local fit)")
ax.set_title("Model fit quality vs k")
ax.legend()
plt.tight_layout()
plt.savefig("figures/stationarity_vs_k.png", dpi=150)
plt.close()
print("figure 3 done")
