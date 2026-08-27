from __future__ import annotations

import math
import unittest

from thermodynamic_waddington.preprocessing import (
    normalize_total, log1p, normalize_expression, scale_velocity_to_expression,
    velocity_in_transformed_space,
)


class NormalizeTotalTests(unittest.TestCase):

    def test_equal_totals_are_unchanged_by_matching_target(self) -> None:
        expression = [[1.0, 2.0, 3.0], [2.0, 4.0, 0.0]]  # both total 6
        normalized, factors = normalize_total(expression, target_sum=6.0)
        for row in normalized:
            self.assertAlmostEqual(sum(row), 6.0, places=6)

    def test_every_cell_totals_to_target_sum_after_normalization(self) -> None:
        expression = [[10.0, 0.0], [1.0, 1.0], [100.0, 300.0]]
        normalized, factors = normalize_total(expression, target_sum=50.0)
        for row in normalized:
            self.assertAlmostEqual(sum(row), 50.0, places=4)

    def test_relative_proportions_within_a_cell_are_preserved(self) -> None:
        expression = [[10.0, 30.0, 60.0]]
        normalized, _ = normalize_total(expression, target_sum=1000.0)
        row = normalized[0]
        self.assertAlmostEqual(row[1] / row[0], 3.0, places=6)
        self.assertAlmostEqual(row[2] / row[0], 6.0, places=6)

    def test_default_target_is_median_of_positive_totals(self) -> None:
        expression = [[10.0], [20.0], [30.0]]  # totals 10, 20, 30 -> median 20
        _, factors = normalize_total(expression)
        # cell with total 20 should have size_factor 1.0 (no change)
        self.assertAlmostEqual(factors[1], 1.0, places=6)

    def test_zero_total_cell_does_not_crash_or_divide_by_zero(self) -> None:
        expression = [[0.0, 0.0], [5.0, 5.0]]
        normalized, factors = normalize_total(expression, target_sum=10.0)
        self.assertTrue(all(math.isfinite(v) for row in normalized for v in row))

    def test_size_factors_scale_with_cell_total(self) -> None:
        expression = [[10.0], [20.0], [40.0]]
        _, factors = normalize_total(expression, target_sum=10.0)
        # factor should be proportional to total / target
        self.assertAlmostEqual(factors[1] / factors[0], 2.0, places=6)
        self.assertAlmostEqual(factors[2] / factors[0], 4.0, places=6)


class Log1pTests(unittest.TestCase):

    def test_log1p_of_zero_is_zero(self) -> None:
        self.assertEqual(log1p([[0.0, 0.0]]), [[0.0, 0.0]])

    def test_log1p_matches_known_values(self) -> None:
        # log1p(e - 1) = 1
        result = log1p([[math.e - 1.0]])
        self.assertAlmostEqual(result[0][0], 1.0, places=6)

    def test_log1p_is_monotonic(self) -> None:
        result = log1p([[1.0, 5.0, 100.0, 1000.0]])[0]
        self.assertTrue(all(result[i] < result[i + 1] for i in range(len(result) - 1)))

    def test_log1p_compresses_a_fixed_absolute_gap_more_at_large_scale(self) -> None:
        # log1p's derivative is 1/(1+x), strictly decreasing -- so the same
        # fixed ADDITIVE step in raw count space produces a smaller gap in
        # log1p space at large x than at small x. (A fixed MULTIPLICATIVE gap
        # is the wrong property to test here: log1p(10x) - log1p(x) actually
        # approaches log(10), a constant, as x grows, so it does not shrink.)
        small_scale_gap = log1p([[11.0]])[0][0] - log1p([[10.0]])[0][0]
        large_scale_gap = log1p([[1001.0]])[0][0] - log1p([[1000.0]])[0][0]
        self.assertGreater(small_scale_gap, large_scale_gap)

    def test_log1p_handles_negative_input_by_flooring_at_zero(self) -> None:
        # should not raise (math.log1p of a negative below -1 would); this
        # module's log1p is documented as expression-only (non-negative)
        result = log1p([[-5.0]])
        self.assertEqual(result[0][0], 0.0)


class NormalizeExpressionTests(unittest.TestCase):

    def test_full_pipeline_matches_expected_scale_after_real_like_counts(self) -> None:
        # mimics the actual failure mode this module was built to fix: raw
        # counts up to several hundred should end up on an O(1)-O(10) scale
        # after normalize+log1p, not left at raw magnitude.
        expression = [[0.0, 50.0, 200.0], [10.0, 0.0, 567.0], [5.0, 5.0, 5.0]]
        normalized, report = normalize_expression(expression)
        flat = [v for row in normalized for v in row]
        self.assertLess(max(flat), 20.0)  # log1p should have compressed this well below the raw 567
        self.assertTrue(report.log_transformed)
        self.assertEqual(report.n_cells, 3)
        self.assertEqual(report.n_genes, 3)

    def test_report_records_pre_and_post_totals_differently(self) -> None:
        expression = [[100.0, 200.0], [10.0, 20.0]]
        _, report = normalize_expression(expression, target_sum=300.0)
        self.assertNotEqual(report.pre_normalization_total_range, report.post_normalization_total_range)

    def test_to_dict_includes_warning_about_scale_assumption(self) -> None:
        _, report = normalize_expression([[1.0, 2.0]])
        payload = report.to_dict()
        self.assertIn("warning", payload)
        self.assertIn("unit-scale", payload["warning"])

    def test_apply_log1p_false_skips_log_transform(self) -> None:
        expression = [[10.0, 20.0], [5.0, 5.0]]
        normalized, report = normalize_expression(expression, apply_log1p=False)
        self.assertFalse(report.log_transformed)
        # without log1p, a cell normalized to target_sum should still sum
        # to that target (only the total-count step ran)
        target = report.target_sum
        for row in normalized:
            self.assertAlmostEqual(sum(row), target, places=3)


class ScaleVelocityTests(unittest.TestCase):

    def test_velocity_scaled_by_same_factors_as_expression(self) -> None:
        velocity = [[10.0, -10.0], [4.0, -4.0]]
        size_factors = [2.0, 0.5]
        scaled = scale_velocity_to_expression(velocity, size_factors)
        self.assertAlmostEqual(scaled[0][0], 5.0)   # 10 / 2.0
        self.assertAlmostEqual(scaled[1][0], 8.0)   # 4 / 0.5

    def test_negative_velocity_stays_negative_after_scaling(self) -> None:
        velocity = [[-6.0]]
        scaled = scale_velocity_to_expression(velocity, [3.0])
        self.assertLess(scaled[0][0], 0.0)

    def test_zero_size_factor_does_not_crash(self) -> None:
        velocity = [[5.0]]
        scaled = scale_velocity_to_expression(velocity, [0.0])
        self.assertTrue(math.isfinite(scaled[0][0]))


class VelocityInTransformedSpaceTests(unittest.TestCase):
    """This is the function that matters most: it was written specifically
    because `scale_velocity_to_expression` alone left real (large-magnitude)
    count-scale velocity incompatible with log1p-transformed expression,
    which silently produced physically meaningless (overflowed) downstream
    results. These tests check the actual scale-consistency property, not
    just that the function runs.
    """

    def test_zero_velocity_gives_zero_transformed_velocity(self) -> None:
        expression = [[10.0, 20.0, 30.0], [5.0, 15.0, 25.0]]
        velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        _, v_transformed, _ = velocity_in_transformed_space(expression, velocity)
        for row in v_transformed:
            for value in row:
                self.assertAlmostEqual(value, 0.0, places=9)

    def test_transformed_velocity_has_comparable_magnitude_to_log_expression(self) -> None:
        # This is the central regression check: even with large raw-count
        # velocity (hundreds, as in real UMI data), the TRANSFORMED velocity
        # must stay on the same order of magnitude as the log1p-transformed
        # expression itself -- not hundreds or thousands of times larger, the
        # exact failure mode that produced 1e19-scale entropy production on
        # real pancreas data before this function existed.
        expression = [[50.0, 200.0, 10.0] for _ in range(5)]
        velocity = [[400.0, -150.0, 300.0] for _ in range(5)]  # large, real-scale magnitudes
        expr_t, vel_t, _ = velocity_in_transformed_space(expression, velocity)
        max_expr = max(abs(v) for row in expr_t for v in row)
        max_vel = max(abs(v) for row in vel_t for v in row)
        # transformed velocity should be within roughly an order of magnitude
        # of the transformed expression scale, not 100-1000x larger
        self.assertLess(max_vel, max_expr * 10 + 1.0)

    def test_positive_velocity_increases_transformed_value_monotonically(self) -> None:
        # log1p is monotonic increasing, so a strictly positive raw velocity
        # (pure upregulation) must give strictly positive transformed velocity
        expression = [[10.0]]
        velocity = [[50.0]]
        _, v_transformed, _ = velocity_in_transformed_space(expression, velocity)
        self.assertGreater(v_transformed[0][0], 0.0)

    def test_negative_velocity_decreases_transformed_value(self) -> None:
        expression = [[100.0]]
        velocity = [[-50.0]]
        _, v_transformed, _ = velocity_in_transformed_space(expression, velocity)
        self.assertLess(v_transformed[0][0], 0.0)

    def test_extrapolated_expression_below_zero_is_floored_not_negative_logged(self) -> None:
        # velocity more negative than current expression would give a
        # negative extrapolated count, which log1p cannot take directly;
        # must not raise and must still return a finite, decreasing value
        expression = [[5.0]]
        velocity = [[-50.0]]  # extrapolated raw count would be -45
        _, v_transformed, _ = velocity_in_transformed_space(expression, velocity)
        self.assertTrue(math.isfinite(v_transformed[0][0]))
        self.assertLessEqual(v_transformed[0][0], 0.0)

    def test_returned_expression_matches_normalize_expression_output(self) -> None:
        expression = [[10.0, 20.0], [30.0, 40.0]]
        velocity = [[1.0, 1.0], [1.0, 1.0]]
        expr_t, _, report_a = velocity_in_transformed_space(expression, velocity, target_sum=100.0)
        expr_direct, report_b = normalize_expression(expression, target_sum=100.0)
        for row_a, row_b in zip(expr_t, expr_direct):
            for a, b in zip(row_a, row_b):
                self.assertAlmostEqual(a, b, places=9)


if __name__ == "__main__":
    unittest.main()
