"""Tests for the calibrated cyclic-fraction irreversibility measure."""
import math
import unittest

import numpy as np

from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.irreversibility import cyclic_irreversibility

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


if __name__ == "__main__":
    unittest.main()
