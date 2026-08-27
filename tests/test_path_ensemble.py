from __future__ import annotations

import unittest

from thermodynamic_waddington.path_ensemble import (
    PathEnsembleConfig,
    PathWorkSample,
    evaluate_path_ensemble,
)


class PathEnsembleTests(unittest.TestCase):
    def _samples(self) -> list[PathWorkSample]:
        return [
            PathWorkSample(f"f-{index}", "protocol-a", "forward", 1.0 + 0.04 * index, provenance={"source": "measured"})
            for index in range(6)
        ] + [
            PathWorkSample(f"r-{index}", "protocol-a", "reverse", -1.0 - 0.04 * index, provenance={"source": "measured"})
            for index in range(6)
        ]

    def test_complete_protocol_bound_ensemble(self) -> None:
        report = evaluate_path_ensemble(self._samples(), PathEnsembleConfig(bootstrap_rounds=50))
        self.assertEqual(report.status, "protocol_bound_ensemble_diagnostic")
        self.assertEqual(report.protocol_count, 1)
        self.assertEqual(report.replicate_count, 12)
        self.assertEqual(report.protocols[0].forward_replicates, 6)
        self.assertIsNotNone(report.protocols[0].delta_f_forward)
        self.assertIsNotNone(report.protocols[0].crooks_crossing)
        self.assertEqual(report.blockers, ())
        self.assertEqual(len(report.fingerprint), 64)

    def test_missing_reverse_paths_are_blocked(self) -> None:
        report = evaluate_path_ensemble(self._samples()[:6])
        self.assertEqual(report.status, "path_ensemble_not_established")
        self.assertIn("protocol-a:insufficient_reverse_replicates", report.blockers)

    def test_duplicate_trajectory_ids_are_blocked(self) -> None:
        samples = self._samples()
        samples[-1] = PathWorkSample("f-0", "protocol-a", "reverse", -1.0)
        report = evaluate_path_ensemble(samples)
        self.assertIn("duplicate_trajectory_id", report.blockers)

    def test_protocols_are_not_pooled(self) -> None:
        samples = self._samples() + [PathWorkSample(f"other-{index}", "protocol-b", "forward", 0.5) for index in range(6)]
        report = evaluate_path_ensemble(samples)
        self.assertEqual(report.protocol_count, 2)
        self.assertIn("multiple_protocols_stratified_never_pooled", report.warnings)


if __name__ == "__main__":
    unittest.main()
