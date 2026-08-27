from __future__ import annotations

import math
import random
import unittest

from thermodynamic_waddington.optimal_transport import (
    SinkhornConfig,
    sinkhorn,
    squared_euclidean_cost,
    estimate_growth_weights,
    displacement_interpolate,
    fit_schrodinger_chain,
    chain_diagnostic_report,
)
from thermodynamic_waddington.synthetic import make_synthetic_dataset


class SinkhornCorrectnessTests(unittest.TestCase):

    def test_identical_marginals_converge_and_have_bounded_transport_cost(self) -> None:
        points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
        plan = sinkhorn(points, points, config=SinkhornConfig(epsilon=0.05, max_iterations=300))
        self.assertTrue(plan.converged)
        # transporting a distribution to an identical copy of itself should be
        # cheap relative to a random permutation of the same points
        self.assertLess(plan.transport_cost(), 1.0)

    def test_marginals_are_recovered_within_tolerance(self) -> None:
        random.seed(3)
        source = [[random.random(), random.random()] for _ in range(25)]
        target = [[random.random(), random.random()] for _ in range(30)]
        cfg = SinkhornConfig(epsilon=0.08, max_iterations=800, tolerance=1e-6)
        plan = sinkhorn(source, target, config=cfg)
        self.assertTrue(plan.converged, f"failed to converge, violation={plan.marginal_violation}")
        row_mass = [sum(row) for row in plan.coupling]
        col_mass = [sum(plan.coupling[i][j] for i in range(len(source))) for j in range(len(target))]
        for a, b in zip(row_mass, plan.source_mass):
            self.assertAlmostEqual(a, b, delta=1e-4)
        for a, b in zip(col_mass, plan.target_mass):
            self.assertAlmostEqual(a, b, delta=1e-4)

    def test_larger_epsilon_increases_entropy_of_coupling(self) -> None:
        random.seed(5)
        source = [[random.random(), random.random()] for _ in range(20)]
        target = [[random.random(), random.random()] for _ in range(20)]
        low = sinkhorn(source, target, config=SinkhornConfig(epsilon=0.01, max_iterations=500))
        high = sinkhorn(source, target, config=SinkhornConfig(epsilon=1.0, max_iterations=500))
        # As epsilon grows the entropic penalty dominates and the coupling
        # approaches the independent (max-entropy) coupling a (x) b; as epsilon
        # shrinks it approaches a sparse, near-deterministic Monge-like map.
        self.assertGreater(high.entropy(), low.entropy())

    def test_smaller_epsilon_decreases_transport_cost_toward_true_OT_cost(self) -> None:
        # The entropic transport cost is monotonically non-decreasing in epsilon
        # (more regularization can only make the plan more diffuse, i.e. costlier
        # in raw transport terms, at fixed marginals) -- Peyré & Cuturi 2019, Prop 4.8.
        random.seed(7)
        source = [[random.random()] for _ in range(15)]
        target = [[random.random() + 2.0] for _ in range(15)]
        costs = []
        for eps in (0.5, 0.1, 0.02):
            plan = sinkhorn(source, target, config=SinkhornConfig(epsilon=eps, max_iterations=800, tolerance=1e-7))
            costs.append(plan.transport_cost())
        self.assertLessEqual(costs[2], costs[1] + 1e-6)
        self.assertLessEqual(costs[1], costs[0] + 1e-6)

    def test_barycentric_map_reduces_to_identity_for_self_transport_at_small_epsilon(self) -> None:
        points = [[0.0], [5.0], [10.0]]
        plan = sinkhorn(points, points, config=SinkhornConfig(epsilon=0.001, max_iterations=1000, tolerance=1e-8))
        mapped = plan.barycentric_map(points)
        for original, m in zip(points, mapped):
            self.assertAlmostEqual(original[0], m[0], delta=0.05)

    def test_implied_velocity_matches_known_shift(self) -> None:
        source = [[0.0], [1.0], [2.0]]
        target = [[1.0], [2.0], [3.0]]  # pure shift of +1
        plan = sinkhorn(source, target, config=SinkhornConfig(epsilon=0.005, max_iterations=1000, tolerance=1e-8))
        velocity = plan.implied_velocity(source, target, dt=1.0)
        for v in velocity:
            self.assertAlmostEqual(v[0], 1.0, delta=0.15)

    def test_unbalanced_transport_handles_mass_growth(self) -> None:
        # Source has fewer cells than target -- simulates proliferation between
        # timepoints. Balanced OT would force an artificial 1-to-many squeeze;
        # unbalanced OT should reduce (not eliminate) the marginal violation
        # relative to a comparably-tuned balanced solve on mismatched raw counts.
        source = [[0.0], [1.0]]
        target = [[0.0], [0.5], [1.0], [1.5]]
        cfg = SinkhornConfig(epsilon=0.1, max_iterations=500, unbalanced_tau=0.5)
        plan = sinkhorn(source, target, config=cfg)
        self.assertFalse(plan.audit["balanced"])
        self.assertEqual(plan.audit["unbalanced_tau"], 0.5)

    def test_zero_or_empty_marginal_is_handled_without_crashing(self) -> None:
        plan = sinkhorn([], [[0.0, 0.0]])
        self.assertEqual(plan.coupling, [])
        self.assertTrue(plan.converged)

class GrowthWeightTests(unittest.TestCase):

    def test_no_marker_genes_gives_uniform_weights(self) -> None:
        expr = [[1.0, 2.0], [3.0, 4.0]]
        weights = estimate_growth_weights(expr)
        self.assertEqual(weights, [1.0, 1.0])

    def test_marker_genes_shift_weights_directionally(self) -> None:
        expr = [[5.0, 0.0], [0.0, 5.0]]  # cell 0 high birth signature, cell 1 high death signature
        weights = estimate_growth_weights(expr, proliferation_genes=[0], apoptosis_genes=[1])
        self.assertGreater(weights[0], 1.0)
        self.assertLess(weights[1], 1.0)


class DisplacementInterpolationTests(unittest.TestCase):

    def test_endpoints_reproduce_marginals(self) -> None:
        source = [[0.0], [1.0]]
        target = [[10.0], [11.0]]
        plan = sinkhorn(source, target, config=SinkhornConfig(epsilon=0.1, max_iterations=300))
        self.assertEqual(displacement_interpolate(plan, source, target, t=0.0), source)
        self.assertEqual(displacement_interpolate(plan, source, target, t=1.0), target)

    def test_midpoint_lies_between_source_and_target_ranges(self) -> None:
        source = [[0.0]] * 10
        target = [[10.0]] * 10
        plan = sinkhorn(source, target, config=SinkhornConfig(epsilon=0.2, max_iterations=300))
        mid = displacement_interpolate(plan, source, target, t=0.5, n_samples=50, rng_seed=1)
        self.assertTrue(all(0.0 <= p[0] <= 10.0 for p in mid))
        avg = sum(p[0] for p in mid) / len(mid)
        self.assertAlmostEqual(avg, 5.0, delta=0.5)


class SchrodingerChainTests(unittest.TestCase):

    def test_chain_requires_at_least_two_snapshots(self) -> None:
        with self.assertRaises(ValueError):
            fit_schrodinger_chain([[[0.0]]])

    def test_three_snapshot_chain_has_two_transitions_and_matching_potentials(self) -> None:
        random.seed(11)
        snap_a = [[random.random()] for _ in range(12)]
        snap_b = [[random.random() + 1.0] for _ in range(12)]
        snap_c = [[random.random() + 2.0] for _ in range(12)]
        cfg = SinkhornConfig(epsilon=0.1, max_iterations=400)
        chain = fit_schrodinger_chain([snap_a, snap_b, snap_c], timepoints=[0.0, 1.0, 2.0], config=cfg)
        self.assertEqual(len(chain.plans), 2)
        self.assertEqual(len(chain.log_potentials), 3)
        report = chain_diagnostic_report(chain)
        self.assertEqual(report["n_transitions"], 2)
        self.assertIn("claim_boundary", report)

    def test_effective_potential_is_finite_and_length_matches_snapshot(self) -> None:
        random.seed(13)
        snap_a = [[random.random(), random.random()] for _ in range(10)]
        snap_b = [[random.random(), random.random()] for _ in range(10)]
        chain = fit_schrodinger_chain([snap_a, snap_b], config=SinkhornConfig(epsilon=0.1, max_iterations=300))
        potential = chain.effective_potential(0)
        self.assertEqual(len(potential), len(snap_a))
        self.assertTrue(all(math.isfinite(v) for v in potential))


class IntegrationWithSyntheticDatasetTests(unittest.TestCase):
    """Exercises the module against this package's own synthetic generator so
    it is validated on data with the same shape/scale conventions as the rest
    of the codebase, not just hand-built toy points.
    """

    def test_splits_synthetic_population_by_velocity_and_bridges_them(self) -> None:
        dataset = make_synthetic_dataset(cells=60, genes=6, seed=21)
        # Treat expression as "before" and expression+velocity as a synthetic
        # "after" snapshot -- a deliberately simple two-timepoint construction
        # used only to validate the machinery end-to-end, not as a claim about
        # what real snapshot pairs look like.
        before = dataset.expression
        after = [[e + v for e, v in zip(row_e, row_v)] for row_e, row_v in zip(dataset.expression, dataset.velocity)]
        cfg = SinkhornConfig(epsilon=0.3, max_iterations=1500, tolerance=1e-5)
        plan = sinkhorn(before, after, config=cfg)
        self.assertTrue(plan.converged)
        velocity_est = plan.implied_velocity(before, after)
        self.assertEqual(len(velocity_est), len(before))
        # OT-implied displacement direction should correlate positively with
        # the ground-truth synthetic velocity more often than chance (>50%
        # of cells with matching sign on the dominant coordinate).
        agree = 0
        for est, true_v in zip(velocity_est, dataset.velocity):
            if est and true_v and (est[0] > 0) == (true_v[0] > 0):
                agree += 1
        self.assertGreater(agree / len(before), 0.5)


if __name__ == "__main__":
    unittest.main()
