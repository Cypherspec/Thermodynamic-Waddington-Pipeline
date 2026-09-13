import unittest

import numpy as np

from thermodynamic_waddington.calibration import (
    basin_barriers_kt,
    boltzmann_free_energy,
    calibrate,
)


def _two_clusters(seed=0):
    rng = np.random.default_rng(seed)
    dense = rng.normal(0.0, 0.3, size=(120, 2))
    sparse = rng.normal(6.0, 1.2, size=(20, 2))
    return np.vstack([dense, sparse])


class TestCalibration(unittest.TestCase):
    def test_boltzmann_shifted_and_finite(self):
        coords = _two_clusters()
        f = boltzmann_free_energy(coords)
        self.assertEqual(len(f), len(coords))
        self.assertTrue(all(np.isfinite(f)))
        self.assertAlmostEqual(min(f), 0.0, places=9)

    def test_sparse_region_is_higher_energy(self):
        coords = _two_clusters()
        f = np.asarray(boltzmann_free_energy(coords))
        dense_mean = f[:120].mean()
        sparse_mean = f[120:].mean()
        self.assertGreater(sparse_mean, dense_mean)

    def test_linear_recovery(self):
        coords = _two_clusters()
        fb = np.asarray(boltzmann_free_energy(coords))
        path = 2.0 * fb + 5.0
        rep = calibrate(path, coords)
        self.assertAlmostEqual(rep.kt_per_work_unit, 0.5, places=6)
        self.assertAlmostEqual(rep.offset_kt, -2.5, places=6)
        self.assertGreater(rep.r_squared, 0.999)

    def test_r_squared_bounded(self):
        coords = _two_clusters()
        rng = np.random.default_rng(1)
        path = rng.normal(size=len(coords))
        rep = calibrate(path, coords)
        self.assertGreaterEqual(rep.r_squared, 0.0)
        self.assertLessEqual(rep.r_squared, 1.0 + 1e-9)

    def test_deterministic(self):
        coords = _two_clusters()
        path = np.asarray(boltzmann_free_energy(coords)) * 1.5
        a = calibrate(path, coords).to_dict()
        b = calibrate(path, coords).to_dict()
        self.assertEqual(a["kt_per_work_unit"], b["kt_per_work_unit"])
        self.assertEqual(a["r_squared"], b["r_squared"])

    def test_basin_barriers(self):
        coords = _two_clusters()
        f = boltzmann_free_energy(coords)
        labels = ["dense"] * 120 + ["sparse"] * 20
        bb = basin_barriers_kt(f, labels)
        self.assertIn("dense", bb)
        self.assertIn("sparse", bb)
        self.assertGreater(bb["sparse"]["mean_kt"], bb["dense"]["mean_kt"])

    def test_needs_variation(self):
        coords = _two_clusters()
        with self.assertRaises(ValueError):
            calibrate([1.0] * len(coords), coords)


if __name__ == "__main__":
    unittest.main()
