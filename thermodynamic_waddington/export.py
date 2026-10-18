"""Export a per-cell thermodynamic table for data platforms.

`export_cells` runs the calibrated readouts on one lineage and writes two files
any tabular or map-based platform can ingest:

    <prefix>_cells.csv    one row per cell: identifier, a one-sentence text
                          description (for platforms that embed text), stage,
                          2-D coordinates, committor, phase, free energy, local
                          cyclic fraction, noise-matched local excess, and the
                          identifiers of each cell's downstream neighbours along
                          the committor (for drawing commitment-flow edges)
    <prefix>_summary.json the lineage-level results: commitment barrier and its
                          potential-of-mean-force curve, where commitment
                          happens, per-stage means, and the noise-matched
                          irreversibility of the whole lineage

Column names start with `tw_` so they cannot collide with a platform's own
fields. The table is self-describing: the summary lists every column with its
meaning and units.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

PHASES = ((0.0, 0.2, "progenitor"), (0.2, 0.4, "early"), (0.4, 0.6, "transition state"),
          (0.6, 0.8, "late"), (0.8, 1.0001, "committed"))

COLUMNS = {
    "tw_id": "cell identifier",
    "tw_text": "one-sentence description of the cell's thermodynamic state",
    "tw_stage": "annotated cell type or stage",
    "tw_x": "2-D coordinate (the embedding passed in, or PCA)",
    "tw_y": "2-D coordinate (the embedding passed in, or PCA)",
    "tw_committor": "probability of reaching the target before the source (0..1)",
    "tw_phase": "committor band: progenitor / early / transition state / late / committed",
    "tw_energy": "effective free energy (path-work landscape units)",
    "tw_cyclic_fraction": "local cyclic (irreversible) fraction of the velocity flow (0..1)",
    "tw_cyclic_excess": "local cyclic fraction above the noise-matched null (~0 = reversible within noise)",
    "tw_downstream": "identifiers of kNN neighbours further along the committor, separated by '|'",
}
# Written only when the export runs a commitment ensemble (ensemble_fits > 0).
ENSEMBLE_COLUMNS = {
    "tw_committor_sd": "spread of the committor across independent refits on 80% subsets (blank if in fewer than two)",
}


@dataclass
class ExportResult:
    cells_path: Path
    summary_path: Path
    summary: dict


def phase_of(q: float) -> str:
    for lo, hi, name in PHASES:
        if lo <= q < hi:
            return name
    return "committed"


def _pca(x: np.ndarray, dims: int) -> np.ndarray:
    xc = x - x.mean(axis=0)
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    return vt[:dims].T


def export_cells(
    expression: Sequence[Sequence[float]],
    velocity: Sequence[Sequence[float]],
    labels: Sequence[str],
    source_labels: Sequence[str],
    target_labels: Sequence[str],
    out_dir: str | Path,
    prefix: str = "tw",
    coords: Sequence[Sequence[float]] | None = None,
    cell_ids: Sequence[str] | None = None,
    stage_order: Sequence[str] | None = None,
    dims: int = 6,
    neighbors: int = 20,
    downstream: int = 3,
    seed: int = 0,
    dataset: str = "",
    ensemble_fits: int = 0,
) -> ExportResult:
    """Run the calibrated readouts on one lineage and write the cell table and summary.

    With `ensemble_fits` > 1 the landscape is also refitted on that many random 80% subsets of the
    cells (`commitment_ensemble`), adding a per-cell `tw_committor_sd` column and the barrier's
    spread across fits to the summary.
    """
    from .commitment_ensemble import commitment_ensemble
    from .config import FitConfig
    from .developmental import (
        RESOLVED_BANDWIDTH,
        barrier_convergence,
        commitment_profile,
        committor_free_energy_profile,
        developmental_coordinate,
        resolved_barrier,
    )
    from .graph import build_knn
    from .irreversibility import cyclic_fraction_per_cell, noise_matched_irreversibility
    from .model import fit_landscape

    X = np.asarray(expression, dtype=float)
    V = np.asarray(velocity, dtype=float)
    labels = [str(x) for x in labels]
    n = len(X)
    if V.shape != X.shape or len(labels) != n:
        raise ValueError("expression, velocity and labels must describe the same cells")
    ids = [str(c) for c in cell_ids] if cell_ids is not None else [f"cell_{i}" for i in range(n)]

    cfg = FitConfig(neighbors=neighbors, dimensions=dims, bootstrap_replicates=2, enable_cycle_decomposition=False,
                    entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2, seed=seed)
    fit = fit_landscape(X.tolist(), V.tolist(), config=cfg, labels=labels, metadata={"velocity_status": "provided"})
    q = np.asarray(developmental_coordinate(fit, source_labels, target_labels))
    energy = np.asarray(fit.energies, dtype=float)
    profile = commitment_profile(fit, source_labels, target_labels)
    pmf = committor_free_energy_profile(q.tolist())
    resolved = committor_free_energy_profile(q.tolist(), bandwidth=RESOLVED_BANDWIDTH)
    corrected = resolved_barrier(q, seed=seed)
    convergence = barrier_convergence(q, bandwidth=RESOLVED_BANDWIDTH, seed=seed)

    comp = _pca(X, dims)
    pts, vp = (X - X.mean(axis=0)) @ comp, V @ comp
    graph = build_knn(pts.tolist(), min(neighbors, n - 1))
    local = cyclic_fraction_per_cell(pts, vp, graph)
    nm = noise_matched_irreversibility(pts, vp, graph, seed=seed, per_cell=True)
    xy = np.asarray(coords, dtype=float)[:, :2] if coords is not None else pts[:, :2]
    ens = None
    if ensemble_fits > 1:
        ens_cfg = FitConfig(neighbors=neighbors, dimensions=dims, bootstrap_replicates=2, enable_cycle_decomposition=False,
                            entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2,
                            seed=seed + 1)
        ens = commitment_ensemble(X, V, labels, source_labels, target_labels, fits=ensemble_fits, config=ens_cfg,
                                  seed=seed + 1)
    columns = {**COLUMNS, **(ENSEMBLE_COLUMNS if ens is not None else {})}

    down = []
    for i, nb in enumerate(graph.neighbors):
        ahead = sorted((j for j in nb if q[j] > q[i] + 1e-4), key=lambda j: -q[j])[:downstream]
        down.append("|".join(ids[j] for j in ahead))

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cells_path = out / f"{prefix}_cells.csv"
    with cells_path.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(list(columns))
        for i in range(n):
            ph = phase_of(q[i])
            text = (f"{labels[i]} cell in the {ph} band of commitment (committor q = {q[i]:.2f}); "
                    f"free energy {energy[i]:.2f}; local irreversibility beyond noise {nm.per_cell_excess[i]:+.3f}.")
            row = [ids[i], text, labels[i], f"{xy[i, 0]:.5f}", f"{xy[i, 1]:.5f}", f"{q[i]:.5f}", ph,
                   f"{energy[i]:.5f}", f"{local[i]:.5f}", f"{nm.per_cell_excess[i]:.5f}", down[i]]
            if ens is not None:
                row.append(f"{ens.q_sd[i]:.5f}" if np.isfinite(ens.q_sd[i]) else "")
            w.writerow(row)

    order = list(stage_order) if stage_order is not None else list(profile.order)
    stages = []
    for lab in order:
        m = np.array([x == lab for x in labels])
        if m.any():
            stages.append({"stage": lab, "n_cells": int(m.sum()), "mean_committor": round(float(q[m].mean()), 4),
                           "mean_energy": round(float(energy[m].mean()), 4),
                           "mean_cyclic_excess": round(float(nm.per_cell_excess[m].mean()), 4)})
    summary = {
        "format": "thermodynamic-waddington-export", "format_version": 1, "dataset": dataset,
        "n_cells": n, "source": list(source_labels), "target": list(target_labels),
        "commitment": {"commitment_stage": profile.commitment_label, "barrier_kt": round(float(pmf.barrier_kt), 4),
                       "transition_state_q": round(float(pmf.barrier_q), 4),
                       "barrier_resolved_kt": round(float(corrected.barrier_kt), 4),
                       "barrier_resolved_se_kt": round(float(corrected.se_kt), 4),
                       "barrier_resolved_raw_kt": round(float(resolved.barrier_kt), 4),
                       "resolved_bandwidth": RESOLVED_BANDWIDTH,
                       "resolved_barrier_drift_kt": round(float(convergence.drift_kt), 4),
                       "pmf_resolved_kt": [round(float(x), 4) for x in resolved.free_energy_kt],
                       "pmf_q": [round(float(x), 4) for x in pmf.q_grid],
                       "pmf_kt": [round(float(x), 4) for x in pmf.free_energy_kt]},
        "irreversibility": {"cyclic_fraction": round(nm.cyclic_fraction, 4), "noise_matched_null": round(nm.null_mean, 4),
                            "noise_matched_excess": round(nm.excess, 4), "z": round(nm.z, 3)},
        "phases": [{"name": name, "q_min": lo, "q_max": min(hi, 1.0),
                    "n_cells": int(sum(phase_of(x) == name for x in q))} for lo, hi, name in PHASES],
        "stages": stages,
        "columns": columns,
        "files": {"cells": cells_path.name},
    }
    if ens is not None:
        summary["ensemble"] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in ens.to_dict().items()
                               if k not in ("barriers_kt", "barrier_se_kt", "sd_by_phase")}
        summary["ensemble"].update({
            "barriers_kt": [round(b, 4) for b in ens.barriers_kt],
            "barrier_se_kt": [round(b, 4) for b in ens.barrier_se_kt],
            "sd_by_phase": {k: round(v, 4) for k, v in ens.sd_by_phase.items()},
        })
    summary_path = out / f"{prefix}_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    return ExportResult(cells_path, summary_path, summary)


def main(argv: Sequence[str] | None = None) -> None:
    """`tw-export`: run the calibrated readouts on an AnnData lineage and write the cell table."""
    import argparse

    ap = argparse.ArgumentParser(prog="tw-export", description="Export a per-cell Thermodynamic Waddington table.")
    ap.add_argument("h5ad", help="AnnData with 'spliced' and 'unspliced' layers")
    ap.add_argument("--label-key", required=True, help="obs column with cell types or stages")
    ap.add_argument("--source", nargs="+", required=True, help="progenitor label(s)")
    ap.add_argument("--target", nargs="+", required=True, help="terminal label(s)")
    ap.add_argument("--lineage", nargs="+", help="labels to keep, in developmental order (default: all)")
    ap.add_argument("--out", default="tw_export", help="output directory")
    ap.add_argument("--prefix", default="tw")
    ap.add_argument("--cells", type=int, default=3000, help="random subsample size (default 3000)")
    ap.add_argument("--genes", type=int, default=100, help="top genes by spliced variance")
    ap.add_argument("--embedding", default="X_umap", help="obsm key for 2-D coordinates (PCA if absent)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ensemble-fits", type=int, default=0,
                    help="also refit on this many 80%% subsets for a per-cell committor spread (0 = off)")
    args = ap.parse_args(argv)

    import anndata as ad

    a = ad.read_h5ad(args.h5ad)
    labels = a.obs[args.label_key].astype(str).to_numpy()
    keep = np.arange(a.n_obs) if not args.lineage else np.where(np.isin(labels, args.lineage))[0]
    if len(keep) > args.cells:
        keep = np.sort(np.random.default_rng(args.seed).choice(keep, size=args.cells, replace=False))
    a = a[keep].copy()

    def dense(m):
        return np.asarray(m.toarray() if hasattr(m, "toarray") else m, dtype=float)

    sp_all = a.layers["spliced"]
    mean = np.asarray(sp_all.mean(axis=0)).ravel()
    sq = sp_all.multiply(sp_all) if hasattr(sp_all, "multiply") else np.asarray(sp_all) ** 2
    top = np.argsort(np.asarray(sq.mean(axis=0)).ravel() - mean * mean)[::-1][: args.genes]
    sp, un = dense(sp_all[:, top]), dense(a.layers["unspliced"][:, top])
    coords = np.asarray(a.obsm[args.embedding]) if args.embedding in a.obsm else None
    res = export_cells(sp, un - sp, a.obs[args.label_key].astype(str).to_numpy(), args.source, args.target, args.out,
                       prefix=args.prefix, coords=coords, cell_ids=list(a.obs_names), stage_order=args.lineage,
                       seed=args.seed, dataset=str(args.h5ad), ensemble_fits=args.ensemble_fits)
    s = res.summary
    print(f"wrote {res.cells_path} and {res.summary_path}: {s['n_cells']} cells, commits at "
          f"{s['commitment']['commitment_stage']}, barrier {s['commitment']['barrier_kt']} kT, "
          f"noise-matched excess {s['irreversibility']['noise_matched_excess']:+.4f}")


if __name__ == "__main__":
    main()
