"""End-to-end tutorial on synthetic data, no download needed.

Runs the whole pipeline and prints each result in order: the landscape, the
entropy-production test, the kT calibration, the developmental coordinate, and
where fate commits.

    python examples/tutorial.py
"""
from __future__ import annotations

from thermodynamic_waddington import (
    FitConfig,
    calibrate_fit,
    commitment_profile,
    developmental_coordinate,
    fit_landscape,
    make_synthetic_dataset,
)


def main() -> None:
    ds = make_synthetic_dataset(cells=200, genes=16, seed=7)
    cfg = FitConfig(neighbors=15, dimensions=6, seed=7)
    fit = fit_landscape(ds.expression, ds.velocity, config=cfg, labels=ds.labels)

    print("1. landscape")
    print(f"   cells={len(fit.energies)}  attractors={len(fit.attractors)}  edges={len(fit.edges)}")

    print("2. is it irreversible?")
    ep = fit.diagnostics["entropy_production"]
    print(f"   entropy production rate={ep.get('entropy_production_rate')}")
    print(f"   permutation p-value={ep.get('permutation_p_value')}")

    print("3. calibrate to kT")
    rep = calibrate_fit(fit)
    print(f"   landscape depth={rep.boltzmann_energy_range_kt:.2f} kT")
    print(f"   path-work vs density R^2={rep.r_squared:.3f}")

    print("4. developmental coordinate and commitment")
    labels = sorted(set(ds.labels))
    source, target = [labels[0]], [labels[-1]]
    q = developmental_coordinate(fit, source, target)
    profile = commitment_profile(fit, source, target)
    print(f"   committor range=[{min(q):.2f}, {max(q):.2f}]")
    print(f"   order={profile.order}")
    print(f"   commitment at={profile.commitment_label}")


if __name__ == "__main__":
    main()
