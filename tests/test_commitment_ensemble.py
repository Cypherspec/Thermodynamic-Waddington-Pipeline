import math
import unittest

import numpy as np

from thermodynamic_waddington import commitment_ensemble
from thermodynamic_waddington.synthetic import make_synthetic_dataset


class CommitmentEnsembleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ds = make_synthetic_dataset(cells=120, genes=16, seed=5)
        cls.labels = list(ds.labels)
        names = sorted(set(cls.labels))
        cls.src, cls.tgt = [names[0]], [names[-1]]
        cls.ens = commitment_ensemble(ds.expression, ds.velocity, cls.labels, cls.src, cls.tgt, fits=4, fraction=0.8, seed=2)
        cls.ds = ds

    def test_every_cell_is_pooled_consistently(self):
        e = self.ens
        n = len(self.labels)
        self.assertEqual(len(e.q_mean), n)
        self.assertEqual(sum(e.n_fits), 4 * round(0.8 * n))
        for m, s, k in zip(e.q_mean, e.q_sd, e.n_fits):
            if k == 0:
                self.assertTrue(math.isnan(m))
                continue
            self.assertGreaterEqual(m, -1e-12)
            self.assertLessEqual(m, 1 + 1e-12)
            if k > 1:
                self.assertGreaterEqual(s, 0.0)
            else:
                self.assertTrue(math.isnan(s))

    def test_boundary_cells_are_certain(self):
        e = self.ens
        lab = np.asarray(self.labels)
        sd = np.asarray(e.q_sd)
        for group in (self.src, self.tgt):
            sel = np.isin(lab, group) & np.isfinite(sd)
            self.assertTrue(np.allclose(sd[sel], 0.0))

    def test_fit_level_summaries(self):
        e = self.ens
        self.assertEqual(len(e.barriers_kt), 4)
        self.assertEqual(len(e.barrier_se_kt), 4)
        self.assertTrue(all(v >= 0.0 for v in e.barrier_se_kt))
        self.assertEqual(len(e.commitment_stages), 4)
        self.assertAlmostEqual(e.barrier_mean_kt, float(np.mean(e.barriers_kt)))
        self.assertGreaterEqual(e.barrier_sd_kt, 0.0)
        self.assertGreater(e.commitment_agreement, 0.0)
        self.assertLessEqual(e.commitment_agreement, 1.0)
        self.assertEqual(set(e.sd_by_phase), {"progenitor", "early", "transition state", "late", "committed"})
        self.assertNotIn("q_sd", e.to_dict())
        self.assertEqual(len(e.pairwise_spearman), 6)
        self.assertTrue(all(-1.0 <= r <= 1.0 for r in e.pairwise_spearman))
        self.assertAlmostEqual(e.rank_agreement, float(np.mean(e.pairwise_spearman)))

    def test_spearman_handles_ties(self):
        from thermodynamic_waddington.commitment_ensemble import _average_ranks, _spearman
        np.testing.assert_allclose(_average_ranks(np.array([0.0, 0.0, 1.0, 0.5])), [0.5, 0.5, 3.0, 2.0])
        x = np.array([0.1, 0.4, 0.4, 0.9, 0.2])
        self.assertAlmostEqual(_spearman(x, x ** 3), 1.0)
        self.assertAlmostEqual(_spearman(x, -x), -1.0)

    def test_deterministic_and_validated(self):
        again = commitment_ensemble(self.ds.expression, self.ds.velocity, self.labels, self.src, self.tgt, fits=4,
                                    fraction=0.8, seed=2)
        np.testing.assert_allclose(again.q_mean, self.ens.q_mean, equal_nan=True)
        with self.assertRaises(ValueError):
            commitment_ensemble(self.ds.expression, self.ds.velocity, self.labels, self.src, self.tgt, fits=1)
        with self.assertRaises(ValueError):
            commitment_ensemble(self.ds.expression, self.ds.velocity, self.labels, self.src, self.tgt, fraction=0)


if __name__ == "__main__":
    unittest.main()
