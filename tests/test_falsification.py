from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from thermodynamic_waddington.falsification import (
    _rank_correlation,
    run_falsification_suite,
    write_falsification_report,
    FalsificationCase,
)


class RankCorrelationTests(unittest.TestCase):
    """`_rank_correlation` is a from-scratch Spearman rank correlation with
    tie handling (average-rank method) -- every downstream falsification case
    depends on it correctly distinguishing "same ordering", "reversed
    ordering", and "unrelated ordering". If this primitive is subtly wrong,
    every falsification case built on top of it is unreliable in a way that
    would not show up as a crash, only as a silently wrong pass/fail -- worth
    validating directly against known Spearman values.
    """

    def test_identical_sequences_give_correlation_one(self) -> None:
        values = [3.0, 1.0, 4.0, 1.5, 5.0, 9.0]
        self.assertAlmostEqual(_rank_correlation(values, values), 1.0, places=9)

    def test_perfectly_reversed_sequences_give_correlation_negative_one(self) -> None:
        a = [1.0, 2.0, 3.0, 4.0, 5.0]
        b = [5.0, 4.0, 3.0, 2.0, 1.0]
        self.assertAlmostEqual(_rank_correlation(a, b), -1.0, places=9)

    def test_known_spearman_example(self) -> None:
        # Classic textbook example (no ties): x ranks 1..5, y ranks give rho=0.6
        x = [1, 2, 3, 4, 5]
        y = [2, 1, 4, 3, 5]
        # d = rank(x)-rank(y) = [-1,1,-1,1,0], sum(d^2)=4, n=5
        # rho = 1 - 6*sum(d^2)/(n*(n^2-1)) = 1 - 24/120 = 0.8
        self.assertAlmostEqual(_rank_correlation(x, y), 0.8, places=6)

    def test_monotonic_nonlinear_transform_still_gives_correlation_one(self) -> None:
        # Spearman is rank-based, so a monotonic (even nonlinear) transform of
        # the same underlying order must still give exactly 1.0
        a = [1.0, 2.0, 3.0, 4.0, 5.0]
        b = [math.exp(v) for v in a]
        self.assertAlmostEqual(_rank_correlation(a, b), 1.0, places=9)

    def test_ties_are_handled_via_average_rank(self) -> None:
        # a has a tie at positions 1,2 (both value 2.0); should not crash and
        # should produce a finite, bounded correlation
        a = [1.0, 2.0, 2.0, 3.0]
        b = [1.0, 2.0, 3.0, 4.0]
        result = _rank_correlation(a, b)
        self.assertTrue(math.isfinite(result))
        self.assertGreaterEqual(result, -1.0)
        self.assertLessEqual(result, 1.0)

    def test_constant_sequence_gives_zero_not_a_crash(self) -> None:
        # all ranks tied -> zero variance in ranks -> da or db is 0 -> function
        # must return 0.0 rather than dividing by zero
        a = [5.0, 5.0, 5.0, 5.0]
        b = [1.0, 2.0, 3.0, 4.0]
        self.assertEqual(_rank_correlation(a, b), 0.0)

    def test_fewer_than_two_points_gives_nan(self) -> None:
        self.assertTrue(math.isnan(_rank_correlation([1.0], [2.0])))
        self.assertTrue(math.isnan(_rank_correlation([], [])))

    def test_uses_shorter_length_when_sequences_differ(self) -> None:
        a = [1.0, 2.0, 3.0]
        b = [1.0, 2.0, 3.0, 999.0, -999.0]  # extra values should be ignored
        self.assertAlmostEqual(_rank_correlation(a, b), 1.0, places=9)


class FalsificationSuiteTests(unittest.TestCase):
    """These run the real suite (through actual fit_landscape calls, same as
    a real user invoking it) rather than mocking it out -- slower, but this
    module's entire purpose is to catch real regressions in the fitted model,
    so testing it against a mock would defeat the point.
    """

    def test_suite_passes_on_the_current_real_implementation(self) -> None:
        result = run_falsification_suite(seed=101)
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["passed"], result["total"])
        for case in result["cases"]:
            self.assertTrue(case["passed"], f"{case['name']} failed: {case}")

    def test_same_seed_gives_identical_fingerprint(self) -> None:
        first = run_falsification_suite(seed=202)
        second = run_falsification_suite(seed=202)
        self.assertEqual(first["fingerprint"], second["fingerprint"])

    def test_different_seed_can_give_different_fingerprint(self) -> None:
        first = run_falsification_suite(seed=1)
        second = run_falsification_suite(seed=2)
        # not strictly guaranteed to differ for every possible seed pair, but
        # for these two the underlying synthetic datasets differ, so the
        # energy vectors (and therefore statistics) should too
        self.assertNotEqual(first["fingerprint"], second["fingerprint"])

    def test_all_six_cases_present(self) -> None:
        result = run_falsification_suite(seed=303)
        names = {case["name"] for case in result["cases"]}
        self.assertEqual(names, {
            "seed_determinism", "velocity_sign_sensitivity", "velocity_permutation_sensitivity",
            "affine_coordinate_stability", "finite_output", "null_velocity_disclosure",
        })

    def test_every_case_carries_the_software_only_claim_boundary(self) -> None:
        result = run_falsification_suite(seed=404)
        for case in result["cases"]:
            self.assertIn("not biological", case["claim_boundary"])
        self.assertIn("not establish biological causality", result["claim_boundary"])

    def test_write_falsification_report_produces_valid_json_file(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "sub" / "report.json"
            payload = write_falsification_report(output=path, seed=505)
            self.assertTrue(path.exists())
            reloaded = json.loads(path.read_text())
            self.assertEqual(reloaded["fingerprint"], payload["fingerprint"])
            self.assertEqual(reloaded["status"], "passed")


class FalsificationCaseTests(unittest.TestCase):

    def test_to_dict_round_trips_all_fields(self) -> None:
        case = FalsificationCase(
            name="example", hypothesis="h", passed=True, statistic=0.5, threshold=0.9, interpretation="i",
        )
        d = case.to_dict()
        self.assertEqual(d["name"], "example")
        self.assertEqual(d["passed"], True)
        self.assertIn("Software stress test only", d["claim_boundary"])


if __name__ == "__main__":
    unittest.main()
