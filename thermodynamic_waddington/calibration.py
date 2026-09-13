"""Physical calibration of the effective free-energy landscape.

The path-work free energy is in arbitrary units because the temperature is a
free parameter of the fit. This module expresses the landscape in units of kT
using the observed cell-density distribution as an equilibrium reference, the
standard Waddington pseudopotential convention

    F_boltz(cell) = -ln P(cell)      [units of kT]

It then fits the path-work landscape to that reference, so path-work energies
can be reported on the same kT scale, and returns how well the two agree (R^2)
as a consistency check between the non-equilibrium path-work landscape and the
equilibrium-density landscape.

Assumptions and limits, stated plainly:
- -ln P is a quasi-steady-state pseudopotential, not a molecular free energy.
  Differentiation is not at equilibrium, so kT here is the thermal-energy scale
  convention, valid up to a steady-state-occupancy assumption.
- Rare or under-sampled states read as high energy. Adequate sampling matters.
- R^2 measures how consistent the path-work landscape is with the density
  reference. It is not a claim that either equals a true molecular free energy.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class CalibrationReport:
    n_cells: int
    kt_per_work_unit: float
    offset_kt: float
    r_squared: float
    rmse_kt: float
    boltzmann_energy_range_kt: float
    calibrated_energy_range_kt: float
    boltzmann_energies_kt: list[float]
    calibrated_energies_kt: list[float]
    assumptions: str = (
        "Density-anchored pseudopotential (F = -ln P) as an equilibrium reference; "
        "kT is the thermal-energy convention, valid up to steady-state occupancy. "
        "Not a molecular free-energy measurement."
    )

    def to_dict(self) -> dict:
        return asdict(self)

    def summary(self) -> dict:
        d = asdict(self)
        d.pop("boltzmann_energies_kt")
        d.pop("calibrated_energies_kt")
        return d


def _kde_density(coords, bandwidth: float | None = None) -> np.ndarray:
    """Gaussian kernel density estimate over an embedding, numpy only.

    Density at each cell is the summed Gaussian kernel from all cells, using
    Scott's-rule bandwidth by default. This is the occupancy density that the
    Boltzmann inversion needs, the same idea as gmx sham over a 2D embedding.
    """
    x = np.asarray(coords, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    n, d = x.shape
    if bandwidth is None:
        spread = float(np.mean(x.std(axis=0))) + 1e-12
        bandwidth = spread * n ** (-1.0 / (d + 4)) + 1e-12
    sq_norm = np.sum(x * x, axis=1)
    dist2 = sq_norm[:, None] + sq_norm[None, :] - 2.0 * (x @ x.T)
    np.clip(dist2, 0.0, None, out=dist2)
    kernel = np.exp(-dist2 / (2.0 * bandwidth * bandwidth))
    return kernel.sum(axis=1) / (n * bandwidth ** d)


def boltzmann_free_energy(coords, bandwidth: float | None = None, floor: float = 1e-12) -> list[float]:
    """Density-anchored free energy in kT: F = -ln(p) over the embedding,
    shifted so the minimum is 0. p is a Gaussian KDE occupancy density."""
    dens = _kde_density(coords, bandwidth)
    dens = np.clip(dens, floor, None)
    p = dens / dens.sum()
    f = -np.log(p)
    return (f - f.min()).tolist()


def calibrate(path_energies, coords, bandwidth: float | None = None, floor: float = 1e-12) -> CalibrationReport:
    """Fit the path-work landscape to the density reference and report kT scale.

    Least-squares maps path-work F to the Boltzmann reference F_boltz = -ln P
    computed over the embedding:
        F_cal = a * F_path + b
    where a is the kT-per-work-unit conversion and b the offset. R^2 is the
    agreement between the path-work landscape and the equilibrium-density
    landscape; a low R^2 means the path-work landscape carries non-equilibrium
    structure that the density landscape does not, which is expected when the
    entropy-production test is significant.
    """
    fp = np.asarray(path_energies, dtype=float)
    fb = np.asarray(boltzmann_free_energy(coords, bandwidth, floor), dtype=float)
    if fp.shape[0] != fb.shape[0]:
        raise ValueError("path_energies and coords must have the same length")
    n = fp.shape[0]
    var = float(np.var(fp))
    if n < 3 or var <= 0.0:
        raise ValueError("need at least 3 cells with variation in path energies to calibrate")
    cov = float(np.mean((fp - fp.mean()) * (fb - fb.mean())))
    a = cov / var
    b = float(fb.mean() - a * fp.mean())
    cal = a * fp + b
    ss_res = float(np.sum((cal - fb) ** 2))
    ss_tot = float(np.sum((fb - fb.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    rmse = float(np.sqrt(np.mean((cal - fb) ** 2)))
    cal_shift = cal - cal.min()
    return CalibrationReport(
        n_cells=n,
        kt_per_work_unit=float(a),
        offset_kt=b,
        r_squared=float(r2),
        rmse_kt=rmse,
        boltzmann_energy_range_kt=float(fb.max() - fb.min()),
        calibrated_energy_range_kt=float(cal_shift.max() - cal_shift.min()),
        boltzmann_energies_kt=fb.tolist(),
        calibrated_energies_kt=cal_shift.tolist(),
    )


def basin_barriers_kt(calibrated_energies, labels, floor: float = 1e-12) -> dict:
    """Report per-label basin depth in kT: mean energy of each label relative to
    the global minimum. A shallow number means the label sits near the floor of
    the landscape (a basin), a large one means it sits high (a ridge)."""
    energies = np.asarray(calibrated_energies, dtype=float)
    base = float(energies.min())
    out = {}
    for lab in sorted(set(labels)):
        mask = np.asarray([l == lab for l in labels], dtype=bool)
        vals = energies[mask]
        if vals.size:
            out[str(lab)] = {
                "n": int(vals.size),
                "mean_kt": float(vals.mean() - base),
                "min_kt": float(vals.min() - base),
            }
    return out


def calibrate_fit(fit, bandwidth: float | None = None) -> CalibrationReport:
    """Convenience: calibrate a LandscapeFit using its energies and the
    embedding as the occupancy space for the density reference."""
    return calibrate(fit.energies, fit.embedding, bandwidth)
