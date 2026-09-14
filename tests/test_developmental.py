import unittest

from thermodynamic_waddington import (
    FitConfig,
    commitment_profile,
    committor_free_energy_profile,
    developmental_coordinate,
    fit_landscape,
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

    def test_missing_labels_raise(self):
        with self.assertRaises(ValueError):
            developmental_coordinate(self.fit, ["nonexistent"], self.tgt)
        with self.assertRaises(ValueError):
            developmental_coordinate(self.fit, self.src, ["nonexistent"])


if __name__ == "__main__":
    unittest.main()
