from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from thermodynamic_waddington.causal_validation import StudyCell, StudyConfig, analyze_study, create_preregistration


class CausalValidationTests(unittest.TestCase):
    def test_missing_control_and_outcome_are_blocked(self):
        report = analyze_study([StudyCell("c1", "d1", "r1", "program", target_probability=0.8)], StudyConfig("Monocyte", "vehicle"))
        self.assertEqual(report.status, "invalid_study_schema")
        self.assertIn("missing_control_intervention", report.blockers)

    def test_donor_paired_effect_needs_all_gates(self):
        cells = []
        for donor, treatment, control in [("d1", 0.80, 0.50), ("d2", 0.82, 0.51), ("d3", 0.79, 0.49)]:
            cells.extend([
                StudyCell(f"{donor}-t", donor, "r1", "program_A", target_probability=treatment, lineage_id=f"{donor}-lt"),
                StudyCell(f"{donor}-c", donor, "r1", "vehicle", target_probability=control, lineage_id=f"{donor}-lc"),
            ])
        report = analyze_study(cells, StudyConfig("Monocyte", "vehicle", minimum_cells_per_group=1, bootstrap_rounds=100))
        self.assertEqual(report.status, "experimentally_supported_candidate")
        self.assertEqual(report.interventions[0].status, "experimentally_supported_candidate")

    def test_preregistration_is_hashed_and_not_claimed_as_done(self):
        payload = create_preregistration(StudyConfig("Monocyte", "vehicle"), ["program_A"])
        self.assertEqual(payload["status"], "ready_for_registration_not_executed")
        self.assertEqual(len(payload["preregistration_hash"]), 64)


if __name__ == "__main__":
    unittest.main()
