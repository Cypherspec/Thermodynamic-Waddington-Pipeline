import unittest

import numpy as np

from thermodynamic_waddington import (
    RESOLVED_BANDWIDTH,
    FitConfig,
    barrier_convergence,
    commitment_profile,
    committor_free_energy_profile,
    developmental_coordinate,
    fit_landscape,
    resolved_barrier,
)
from thermodynamic_waddington.synthetic import make_synthetic_dataset


def _fit():
    ds = make_synthetic_dataset(cells=120, genes=16, seed=5)
    cfg = FitConfig(
        neighbors=12, dimensions=6, bootstrap_replicates=2,
        enable_cycle_decomposition=False,
        entropy_production_bootstrap_replicates=2,
        entropy_production_permutation_replicates=2, seed=5,
    )
    fit = fit_landscape(ds.expression, ds.velocity, config=cfg, labels=ds.labels)
    return fit, sorted(set(ds.labels))


class TestDevelopmental(unittest.TestCase):
    def setUp(self):
        self.fit, self.labels = _fit()
        self.src = [self.labels[0]]
        self.tgt = [self.labels[-1]]

    def test_committor_bounded(self):
        q = developmental_coordinate(self.fit, self.src, self.tgt)
        self.assertEqual(len(q), len(self.fit.energies))
        self.assertTrue(all(0.0 <= v <= 1.0 for v in q))

    def test_boundary_conditions(self):
        q = developmental_coordinate(self.fit, self.src, self.tgt)
        labs = list(self.fit.labels)
        for i, lab in enumerate(labs):
            if lab == self.src[0]:
                self.assertAlmostEqual(q[i], 0.0, places=9)
            if lab == self.tgt[0]:
                self.assertAlmostEqual(q[i], 1.0, places=9)

    def test_source_below_target(self):
        rep = commitment_profile(self.fit, self.src, self.tgt)
        self.assertLess(rep.mean_committor[self.src[0]], rep.mean_committor[self.tgt[0]])
        self.assertEqual(rep.order[0], self.src[0])
        self.assertEqual(rep.order[-1], self.tgt[0])

    def test_deterministic(self):
        a = developmental_coordinate(self.fit, self.src, self.tgt)
        b = developmental_coordinate(self.fit, self.src, self.tgt)
        self.assertEqual(a, b)

    def test_free_energy_profile(self):
        q = developmental_coordinate(self.fit, self.src, self.tgt)
        prof = committor_free_energy_profile(q, temperature=1.0, grid=40)
        self.assertEqual(len(prof.q_grid), 40)
        self.assertEqual(len(prof.free_energy_kt), 40)
        self.assertGreaterEqual(prof.barrier_kt, 0.0)
        self.assertGreaterEqual(prof.barrier_q, 0.0)
        self.assertLessEqual(prof.barrier_q, 1.0)
        self.assertAlmostEqual(min(prof.free_energy_kt), 0.0, places=9)

    def test_resolved_profile_reflects_and_validates(self):
        q = developmental_coordinate(self.fit, self.src, self.tgt)
        prof = committor_free_energy_profile(q, bandwidth=RESOLVED_BANDWIDTH)
        self.assertEqual(len(prof.q_grid), 60)
        self.assertAlmostEqual(min(prof.free_energy_kt), 0.0, places=9)
        self.assertGreaterEqual(prof.barrier_kt, 0.0)
        with self.assertRaises(ValueError):
            committor_free_energy_profile(q, bandwidth=0.0)

    def test_resolved_barrier_does_not_drift_with_cell_number(self):
        # A committor distribution with atoms at 0 and 1 and a sparse transition region, as on real data.
        rng = np.random.default_rng(3)
        q = np.concatenate([np.zeros(1200), np.ones(800), rng.beta(0.7, 0.7, 2000)])
        silverman = barrier_convergence(q, sizes=(200, 800, 3200), bandwidth=None, replicates=30)
        resolved = barrier_convergence(q, sizes=(200, 800, 3200), bandwidth=RESOLVED_BANDWIDTH, replicates=30)
        self.assertEqual(resolved.sizes, [200, 800, 3200, 4000])
        self.assertGreater(silverman.drift_kt, 0.3)
        self.assertLess(abs(resolved.drift_kt), 0.15)
        for lo, m, hi in zip(resolved.low_kt, resolved.mean_kt, resolved.high_kt):
            self.assertLessEqual(lo, m)
            self.assertLessEqual(m, hi)

    def test_resolved_barrier_bias_correction(self):
        # Double-well-like committor: atoms at 0 and 1 and a sparse middle, small sample.
        rng = np.random.default_rng(8)
        q = np.concatenate([np.zeros(150), np.ones(120), rng.beta(0.6, 0.6, 130)])
        r = resolved_barrier(q, bootstraps=40, seed=1)
        raw = committor_free_energy_profile(q, bandwidth=RESOLVED_BANDWIDTH)
        self.assertAlmostEqual(r.raw_kt, raw.barrier_kt, places=12)
        self.assertAlmostEqual(r.raw_kt - r.bias_kt, r.barrier_kt, places=12)
        self.assertGreater(r.se_kt, 0.0)
        self.assertEqual(r.n_cells, 400)
        self.assertEqual(r, resolved_barrier(q, bootstraps=40, seed=1))

    def test_missing_labels_raise(self):
        with self.assertRaises(ValueError):
            developmental_coordinate(self.fit, ["nonexistent"], self.tgt)
        with self.assertRaises(ValueError):
            developmental_coordinate(self.fit, self.src, ["nonexistent"])


if __name__ == "__main__":
    unittest.main()
