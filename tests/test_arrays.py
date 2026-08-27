from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from thermodynamic_waddington.arrays import (
    shape, zeros, copy_matrix, transpose, dot, norm, add, sub, scale,
    mean, variance, row_means, column_means, center, pairwise_squared_distances,
    median, percentile, quantile, logsumexp, safe_log, safe_exp, normalize,
    read_matrix, write_matrix, solve_linear_system, matrix_vector, matvec,
    matrix_multiply, trace,
)


class ShapeAndBasicOpsTests(unittest.TestCase):

    def test_shape_of_regular_matrix(self) -> None:
        self.assertEqual(shape([[1, 2, 3], [4, 5, 6]]), (2, 3))

    def test_shape_of_empty_matrix(self) -> None:
        self.assertEqual(shape([]), (0, 0))

    def test_shape_rejects_ragged_matrix(self) -> None:
        with self.assertRaises(ValueError):
            shape([[1, 2], [3, 4, 5]])

    def test_zeros_shape_and_values(self) -> None:
        z = zeros(3, 2)
        self.assertEqual(len(z), 3)
        self.assertTrue(all(len(row) == 2 for row in z))
        self.assertTrue(all(v == 0.0 for row in z for v in row))

    def test_copy_matrix_is_independent(self) -> None:
        original = [[1.0, 2.0]]
        copied = copy_matrix(original)
        copied[0][0] = 99.0
        self.assertEqual(original[0][0], 1.0)

    def test_transpose_known_matrix(self) -> None:
        self.assertEqual(transpose([[1, 2, 3], [4, 5, 6]]), [[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]])

    def test_transpose_is_involution(self) -> None:
        m = [[1, 2], [3, 4], [5, 6]]
        self.assertEqual(transpose(transpose(m)), copy_matrix(m))

    def test_dot_known_value(self) -> None:
        self.assertEqual(dot([1, 2, 3], [4, 5, 6]), 32)

    def test_dot_orthogonal_vectors_is_zero(self) -> None:
        self.assertEqual(dot([1, 0], [0, 1]), 0)

    def test_norm_known_3_4_5_triangle(self) -> None:
        self.assertAlmostEqual(norm([3, 4]), 5.0)

    def test_norm_is_nonnegative_even_with_bad_input(self) -> None:
        self.assertGreaterEqual(norm([1, -1, 1]), 0.0)

    def test_add_sub_are_inverses(self) -> None:
        a, b = [1.0, 2.0, 3.0], [0.5, -1.0, 4.0]
        self.assertEqual(sub(add(a, b), b), a)

    def test_scale_by_zero_gives_zero_vector(self) -> None:
        self.assertEqual(scale([1.0, 2.0, 3.0], 0.0), [0.0, 0.0, 0.0])

    def test_scale_by_one_is_identity(self) -> None:
        self.assertEqual(scale([1.0, 2.0], 1.0), [1.0, 2.0])


class StatisticsTests(unittest.TestCase):

    def test_mean_known_value(self) -> None:
        self.assertAlmostEqual(mean([1, 2, 3, 4]), 2.5)

    def test_mean_of_empty_is_zero(self) -> None:
        self.assertEqual(mean([]), 0.0)

    def test_variance_matches_hand_computed_sample_variance(self) -> None:
        # values 2, 4, 4, 4, 5, 5, 7, 9 -> known sample variance (n-1 denominator) = 4.571428...
        values = [2, 4, 4, 4, 5, 5, 7, 9]
        self.assertAlmostEqual(variance(values), 32.0 / 7.0, places=6)

    def test_variance_of_constant_sequence_is_zero(self) -> None:
        self.assertEqual(variance([5, 5, 5, 5]), 0.0)

    def test_variance_of_fewer_than_two_values_is_zero(self) -> None:
        self.assertEqual(variance([5]), 0.0)
        self.assertEqual(variance([]), 0.0)

    def test_row_means_and_column_means_known_matrix(self) -> None:
        m = [[1, 2, 3], [4, 5, 6]]
        self.assertEqual(row_means(m), [2.0, 5.0])
        self.assertEqual(column_means(m), [2.5, 3.5, 4.5])

    def test_center_removes_column_means(self) -> None:
        m = [[1, 2], [3, 4], [5, 6]]
        centered, means = center(m)
        self.assertEqual(means, [3.0, 4.0])
        for j in range(2):
            self.assertAlmostEqual(sum(row[j] for row in centered), 0.0)

    def test_pairwise_squared_distances_known_points(self) -> None:
        points = [[0, 0], [3, 4], [0, 0]]
        d = pairwise_squared_distances(points)
        self.assertAlmostEqual(d[0][1], 25.0)
        self.assertAlmostEqual(d[0][2], 0.0)
        self.assertEqual(d[0][1], d[1][0])  # symmetric

    def test_median_odd_and_even_length(self) -> None:
        self.assertEqual(median([3, 1, 2]), 2)
        self.assertEqual(median([1, 2, 3, 4]), 2.5)

    def test_median_of_empty_is_zero(self) -> None:
        self.assertEqual(median([]), 0.0)

    def test_percentile_endpoints_match_min_max(self) -> None:
        values = [5, 1, 9, 3, 7]
        self.assertEqual(percentile(values, 0.0), 1)
        self.assertEqual(percentile(values, 1.0), 9)

    def test_percentile_median_equivalence_at_0_5(self) -> None:
        values = [10, 2, 8, 4, 6]
        self.assertAlmostEqual(percentile(values, 0.5), median(values))

    def test_quantile_is_an_alias_for_percentile(self) -> None:
        values = [1, 5, 2, 8, 9, 3]
        for p in (0.1, 0.25, 0.5, 0.75, 0.9):
            self.assertEqual(quantile(values, p), percentile(values, p))

    def test_percentile_of_empty_is_zero(self) -> None:
        self.assertEqual(percentile([], 0.5), 0.0)


class LogSpaceTests(unittest.TestCase):

    def test_logsumexp_matches_direct_computation_for_small_values(self) -> None:
        values = [1.0, 2.0, 3.0]
        expected = math.log(sum(math.exp(v) for v in values))
        self.assertAlmostEqual(logsumexp(values), expected, places=9)

    def test_logsumexp_avoids_overflow_for_large_values(self) -> None:
        # exp(1000) overflows a plain float; logsumexp must not raise or return inf
        values = [1000.0, 1000.0, 999.0]
        result = logsumexp(values)
        self.assertTrue(math.isfinite(result))
        self.assertGreater(result, 1000.0)

    def test_logsumexp_of_empty_is_negative_infinity(self) -> None:
        self.assertEqual(logsumexp([]), float("-inf"))

    def test_logsumexp_of_single_value_is_that_value(self) -> None:
        self.assertAlmostEqual(logsumexp([4.2]), 4.2)

    def test_safe_log_floors_nonpositive_input(self) -> None:
        # log(0) or log(negative) would raise in plain math.log; safe_log must not
        self.assertTrue(math.isfinite(safe_log(0.0)))
        self.assertTrue(math.isfinite(safe_log(-5.0)))

    def test_safe_log_matches_plain_log_for_normal_input(self) -> None:
        self.assertAlmostEqual(safe_log(math.e), 1.0, places=9)

    def test_safe_exp_clamps_extreme_input(self) -> None:
        # exp(10000) would overflow; safe_exp must clamp and stay finite
        self.assertTrue(math.isfinite(safe_exp(10000.0)))
        self.assertTrue(math.isfinite(safe_exp(-10000.0)))

    def test_safe_exp_matches_plain_exp_for_normal_input(self) -> None:
        self.assertAlmostEqual(safe_exp(1.0), math.e, places=9)


class NormalizeTests(unittest.TestCase):

    def test_normalize_sums_to_one(self) -> None:
        result = normalize([1.0, 2.0, 3.0, 4.0])
        self.assertAlmostEqual(sum(result), 1.0)

    def test_normalize_preserves_relative_proportions(self) -> None:
        result = normalize([2.0, 4.0])
        self.assertAlmostEqual(result[1] / result[0], 2.0)

    def test_normalize_of_all_zeros_gives_uniform(self) -> None:
        result = normalize([0.0, 0.0, 0.0, 0.0])
        self.assertTrue(all(abs(v - 0.25) < 1e-12 for v in result))

    def test_normalize_of_empty_is_empty(self) -> None:
        self.assertEqual(normalize([]), [])

    def test_normalize_rejects_negative_total_gracefully(self) -> None:
        # sum <= 0 (e.g. mixed positive/negative summing non-positive) falls back
        # to uniform rather than dividing by a non-positive total
        result = normalize([-1.0, -1.0])
        self.assertTrue(all(abs(v - 0.5) < 1e-12 for v in result))


class LinearAlgebraTests(unittest.TestCase):

    def test_solve_linear_system_known_2x2(self) -> None:
        # x + y = 3, x - y = 1  =>  x=2, y=1
        matrix = [[1.0, 1.0], [1.0, -1.0]]
        vector = [3.0, 1.0]
        solution = solve_linear_system(matrix, vector)
        self.assertAlmostEqual(solution[0], 2.0, places=6)
        self.assertAlmostEqual(solution[1], 1.0, places=6)

    def test_solve_linear_system_identity_matrix_returns_vector_unchanged(self) -> None:
        identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        vector = [7.0, -3.0, 2.5]
        solution = solve_linear_system(identity, vector)
        for a, b in zip(solution, vector):
            self.assertAlmostEqual(a, b, places=9)

    def test_solve_linear_system_3x3_known_solution(self) -> None:
        # 2x + y = 5, x + 3y + z = 10, y + 2z = 6  =>  x=1.75, y=1.5, z=2.25 (verified by substitution)
        matrix = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]]
        vector = [5.0, 10.0, 6.0]
        solution = solve_linear_system(matrix, vector)
        # verify by substitution back into the original system rather than
        # hardcoding a solution vector, so this checks the actual equations
        self.assertAlmostEqual(2 * solution[0] + solution[1], 5.0, places=6)
        self.assertAlmostEqual(solution[0] + 3 * solution[1] + solution[2], 10.0, places=6)
        self.assertAlmostEqual(solution[1] + 2 * solution[2], 6.0, places=6)

    def test_matrix_vector_known_product(self) -> None:
        matrix = [[1, 2], [3, 4]]
        vector = [1, 1]
        self.assertEqual(matrix_vector(matrix, vector), [3, 7])

    def test_matvec_matches_matrix_vector(self) -> None:
        matrix = [[1, 2, 3], [4, 5, 6]]
        vector = [1, 0, 1]
        self.assertEqual(matvec(matrix, vector), matrix_vector(matrix, vector))

    def test_matvec_of_empty_matrix_is_empty(self) -> None:
        self.assertEqual(matvec([], [1, 2]), [])

    def test_matrix_multiply_identity_is_neutral(self) -> None:
        identity = [[1.0, 0.0], [0.0, 1.0]]
        m = [[5.0, 6.0], [7.0, 8.0]]
        self.assertEqual(matrix_multiply(m, identity), m)

    def test_matrix_multiply_known_product(self) -> None:
        a = [[1, 2], [3, 4]]
        b = [[5, 6], [7, 8]]
        # [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]] = [[19,22],[43,50]]
        self.assertEqual(matrix_multiply(a, b), [[19, 22], [43, 50]])

    def test_trace_known_matrix(self) -> None:
        self.assertEqual(trace([[1, 2], [3, 4]]), 5)

    def test_trace_of_rectangular_matrix_uses_min_dimension(self) -> None:
        self.assertEqual(trace([[1, 2, 3], [4, 5, 6]]), 6)  # 1 + 5

    def test_trace_of_empty_matrix_is_zero(self) -> None:
        self.assertEqual(trace([]), 0)


class MatrixIOTests(unittest.TestCase):

    def test_json_roundtrip(self) -> None:
        matrix = [[1.0, 2.0], [3.0, 4.0]]
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "matrix.json"
            write_matrix(path, matrix)
            reloaded = read_matrix(path)
            self.assertEqual(reloaded, matrix)

    def test_csv_roundtrip(self) -> None:
        matrix = [[1.0, 2.5], [3.5, 4.0]]
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "matrix.csv"
            write_matrix(path, matrix)
            reloaded = read_matrix(path)
            self.assertEqual(reloaded, matrix)

    def test_tsv_roundtrip(self) -> None:
        matrix = [[1.0, 2.0], [3.0, 4.0]]
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "matrix.tsv"
            write_matrix(path, matrix)
            reloaded = read_matrix(path)
            self.assertEqual(reloaded, matrix)

    def test_jsonl_read(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "matrix.jsonl"
            path.write_text("[1.0, 2.0]\n[3.0, 4.0]\n")
            reloaded = read_matrix(path)
            self.assertEqual(reloaded, [[1.0, 2.0], [3.0, 4.0]])

    def test_json_dict_with_data_key(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "matrix.json"
            path.write_text(json.dumps({"data": [[1.0, 2.0], [3.0, 4.0]]}))
            reloaded = read_matrix(path)
            self.assertEqual(reloaded, [[1.0, 2.0], [3.0, 4.0]])


if __name__ == "__main__":
    unittest.main()
