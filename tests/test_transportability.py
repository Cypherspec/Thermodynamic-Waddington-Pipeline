from __future__ import annotations

import unittest

from thermodynamic_waddington.transportability import (
    StudyObservation,
    TransportabilityConfig,
    evaluate_transportability,
)


class TransportabilityTests(unittest.TestCase):
    def _study(self, study_id: str, shift: float = 0.0, negative_control: bool = False) -> StudyObservation:
        predictions = tuple(float(index) for index in range(30))
        outcomes = tuple(float(index) + shift for index in range(30))
        return StudyObservation(study_id, predictions, outcomes, domain="hematopoietic", provenance={"source": study_id}, negative_control=negative_control)

    def test_independent_studies_are_equal_study_aggregated(self) -> None:
        report = evaluate_transportability([self._study("a"), self._study("b", 2.0)])
        self.assertEqual(report.status, "transportable_predictive_signal_under_audited_benchmark")
        self.assertEqual(report.heterogeneity["studies"], 2.0)
        self.assertAlmostEqual(report.pooled_equal_study["spearman"], 1.0)
        self.assertTrue(report.fingerprint)

    def test_duplicate_studies_are_blocked(self) -> None:
        report = evaluate_transportability([self._study("a"), self._study("a")])
        self.assertIn("duplicate_study_id", report.blockers)
        self.assertEqual(report.status, "transportability_not_established")

    def test_negative_control_is_a_gate(self) -> None:
        predictions = tuple(float(index) for index in range(30))
        outcomes = tuple(float(29 - index) for index in range(30))
        control = StudyObservation("control", predictions, outcomes, provenance={"source": "null"}, negative_control=True)
        report = evaluate_transportability([self._study("a"), control], TransportabilityConfig(max_negative_control_abs_spearman=0.2))
        self.assertIn("negative_control_signal_detected", report.blockers)


if __name__ == "__main__":
    unittest.main()
