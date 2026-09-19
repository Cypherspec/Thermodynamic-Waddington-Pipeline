"""Tests for the calibrated cyclic-fraction irreversibility measure."""
import math
import unittest

import numpy as np

from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.irreversibility import (
    cycle_affinities,
    cyclic_flow_per_cell,
    cyclic_irreversibility,
)

J = np.array([[0.0, -1.0], [1.0, 0.0]])


def _sample(omega, seed, n=400, a=1.0, D=1.0):
    rng = np.random.default_rng(seed)
    pts = rng.normal(0.0, math.sqrt(D / a), size=(n, 2))
    B = a * np.eye(2) + omega * J
    return pts, -(pts @ B.T)


def _fraction(omega, seed):
    pts, vel = _sample(omega, seed)
    g = build_knn(pts.tolist(), 20)
    return cyclic_irreversibility(pts, vel, g, seed=seed).cyclic_fraction


class TestCyclicIrreversibility(unittest.TestCase):
    def test_bounded(self):
        for w in (0.0, 1.0, 4.0):
            f = _fraction(w, 0)
            self.assertGreaterEqual(f, 0.0)
            self.assertLessEqual(f, 1.0)

    def test_gradient_field_is_low(self):
        # a conservative (equilibrium) field sits near the discretization floor
        self.assertLess(_fraction(0.0, 1), 0.2)

    def test_rotational_field_is_high(self):
        # a strongly driven rotational field is clearly irreversible
        self.assertGreater(_fraction(4.0, 1), 0.5)

    def test_monotonic_in_drive(self):
        vals = [_fraction(w, 2) for w in (0.0, 0.5, 1.0, 2.0, 4.0)]
        for i in range(len(vals) - 1):
            self.assertLessEqual(vals[i], vals[i + 1] + 1e-6)

    def test_separates_equilibrium_from_noneq(self):
        eq = [_fraction(0.0, s) for s in range(4)]
        neq = [_fraction(2.0, s) for s in range(4)]
        self.assertLess(max(eq), min(neq))

    def test_bootstrap_ci(self):
        pts, vel = _sample(1.0, 3)
        g = build_knn(pts.tolist(), 20)
        rep = cyclic_irreversibility(pts, vel, g, bootstrap=50, seed=0)
        self.assertIsNotNone(rep.bootstrap_ci)
        lo, hi = rep.bootstrap_ci
        self.assertLessEqual(lo, rep.cyclic_fraction + 1e-9)
        self.assertLessEqual(rep.cyclic_fraction - 1e-9, hi + 0.2)

    def test_degenerate_inputs(self):
        pts = np.zeros((2, 2))
        g = build_knn([[0.0, 0.0], [1.0, 0.0]], 1)
        rep = cyclic_irreversibility(pts, np.zeros((2, 2)), g)
        self.assertEqual(rep.cyclic_fraction, 0.0)

    def test_sparse_matches_dense(self):
        # the scalable sparse projection must reproduce the dense pseudo-inverse
        from thermodynamic_waddington.irreversibility import (
            _edge_flow, _frac_dense, _frac_sparse, _incidence, _laplacian_pinv, _pairs,
        )
        pts, vel = _sample(1.0, 5, n=500)
        g = build_knn(pts.tolist(), 20)
        pairs = _pairs(g)
        f = _edge_flow(pts, vel, pairs)
        sparse = _frac_sparse(_incidence(pairs, 500), f)
        dense = _frac_dense(pairs, f, 500, _laplacian_pinv(pairs, 500))
        self.assertAlmostEqual(sparse, dense, places=6)


def _affinity(omega, seed):
    pts, vel = _sample(omega, seed)
    g = build_knn(pts.tolist(), 20)
    return cycle_affinities(pts, vel, g).rms_affinity_kt


class TestCycleAffinities(unittest.TestCase):
    def test_nonneg_and_has_cycles(self):
        pts, vel = _sample(1.0, 0)
        g = build_knn(pts.tolist(), 20)
        rep = cycle_affinities(pts, vel, g)
        self.assertGreater(rep.n_cycles, 0)
        self.assertGreaterEqual(rep.rms_affinity_kt, 0.0)
        self.assertGreaterEqual(rep.max_affinity_kt, rep.mean_abs_affinity_kt)

    def test_rotational_exceeds_gradient(self):
        self.assertGreater(_affinity(4.0, 1), 3.0 * _affinity(0.0, 1))

    def test_monotonic_in_drive(self):
        vals = [_affinity(w, 2) for w in (0.0, 0.5, 1.0, 2.0, 4.0)]
        for i in range(len(vals) - 1):
            self.assertLessEqual(vals[i], vals[i + 1] + 1e-6)

    def test_degenerate(self):
        g = build_knn([[0.0, 0.0], [1.0, 0.0]], 1)
        rep = cycle_affinities(np.zeros((2, 2)), np.zeros((2, 2)), g)
        self.assertEqual(rep.n_cycles, 0)
        self.assertEqual(rep.rms_affinity_kt, 0.0)


class TestPerCellMap(unittest.TestCase):
    def test_shape_nonneg_normalized(self):
        pts, vel = _sample(2.0, 0)
        g = build_knn(pts.tolist(), 20)
        per = cyclic_flow_per_cell(pts, vel, g)
        self.assertEqual(per.shape, (len(pts),))
        self.assertTrue((per >= 0).all())
        self.assertAlmostEqual(float(per.sum()), 1.0, places=6)

    def test_degenerate_returns_zeros(self):
        g = build_knn([[0.0, 0.0], [1.0, 0.0]], 1)
        per = cyclic_flow_per_cell(np.zeros((2, 2)), np.zeros((2, 2)), g)
        self.assertEqual(per.shape, (2,))
        self.assertEqual(float(per.sum()), 0.0)


if __name__ == "__main__":
    unittest.main()
