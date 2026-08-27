from __future__ import annotations

import unittest

from thermodynamic_waddington.decision_engine import DecisionConfig, decide


class DecisionEngineTests(unittest.TestCase):
    def test_missing_evidence_does_not_upgrade_claim(self) -> None:
        report = decide()
        self.assertEqual(report.decision, "do_not_upgrade_claim")
        self.assertTrue(report.blockers)
        self.assertEqual(len(report.fingerprint), 64)

    def test_complete_evidence_can_advance(self) -> None:
        trial = {
            "status": "measured_trial_complete",
            "preregistration_hash": "abc",
            "aggregate": {
                "donor_count": 4,
                "effect": {"estimate": 0.22},
                "viability": {"control": 0.90, "treatment": 0.88},
            },
            "blockers": [],
            "warnings": [],
        }
        transportability = {
            "metrics": [{"spearman": 0.50}, {"spearman": 0.42}],
            "pooled_equal_study": {"spearman": 0.46},
            "negative_controls": {"max_abs_spearman": 0.05},
            "blockers": [],
        }
        robustness = {
            "summary": {"finite_fraction": 1.0, "minimum_energy_spearman": 0.90},
            "blockers": [],
        }
        path_ensemble = {
            "status": "protocol_bound_ensemble_diagnostic",
            "protocols": [{"protocol_id": "p", "forward_replicates": 5, "reverse_replicates": 5}],
            "blockers": [],
            "warnings": [],
        }
        report = decide(trial=trial, transportability=transportability, robustness=robustness, path_ensemble=path_ensemble, config=DecisionConfig())
        self.assertEqual(report.decision, "advance_to_replicated_intervention")
        self.assertEqual(report.blockers, ())


if __name__ == "__main__":
    unittest.main()
