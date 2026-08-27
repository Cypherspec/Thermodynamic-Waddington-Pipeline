from __future__ import annotations

import unittest

from thermodynamic_waddington.intervention_trial import TrialConfig, TrialObservation, analyze_trial, preregistration


class InterventionTrialTests(unittest.TestCase):
    def _rows(self):
        rows = []
        for donor in ("d1", "d2", "d3"):
            for arm, value in (("control", 0.2), ("treatment", 0.7)):
                for index in range(8):
                    rows.append(TrialObservation(
                        cell_id=f"{donor}-{arm}-{index}",
                        donor_id=donor,
                        intervention=arm,
                        fate_probability=value,
                        viability=0.95,
                        lineage_id=f"{donor}-lineage-{index}",
                        barcode_id=f"{donor}-barcode-{index}",
                        orthogonal_fate="target" if value > 0.5 else "other",
                    ))
        return rows

    def test_missing_measurements_are_blocked(self):
        report = analyze_trial(self._rows())
        self.assertEqual(report.status, "blocked_pending_measured_experiment")
        self.assertIn("preregistration hash is missing", report.blockers)

    def test_complete_fixture_is_donor_level(self):
        config = TrialConfig(bootstrap_rounds=32)
        registration = preregistration(config)
        report = analyze_trial(self._rows(), config, registration["preregistration_hash"])
        self.assertEqual(report.status, "measured_intervention_result_ready_for_replication")
        self.assertFalse(report.aggregate["cell_pooling_used"])
        self.assertEqual(report.aggregate["donor_count"], 3)
        self.assertAlmostEqual(report.aggregate["effect"]["estimate"], 0.5)

    def test_preregistration_is_stable(self):
        config = TrialConfig(bootstrap_rounds=32)
        self.assertEqual(preregistration(config)["preregistration_hash"], preregistration(config)["preregistration_hash"])


if __name__ == "__main__":
    unittest.main()
