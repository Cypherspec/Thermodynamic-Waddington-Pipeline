from __future__ import annotations
import unittest
from thermodynamic_waddington.robustness import evaluate_robustness, run_synthetic_robustness
from thermodynamic_waddington.synthetic import make_synthetic_dataset

class RobustnessTests(unittest.TestCase):
    def test_synthetic_report_is_deterministic(self):
        first = run_synthetic_robustness(seed=12, cells=32, genes=6)
        second = run_synthetic_robustness(seed=12, cells=32, genes=6)
        self.assertEqual(first["fingerprint"], second["fingerprint"])
        self.assertEqual(first["summary"]["scenario_count"], 9)
    def test_failed_scenario_is_a_blocker(self):
        data = make_synthetic_dataset(cells=24, genes=5, seed=3).to_dict()
        report = evaluate_robustness(data["expression"], data["velocity"], scenarios=[{"name": "invalid", "neighbors": 999}])
        self.assertIn("nonfinite_or_failed_sensitivity_scenario", report["blockers"])

if __name__ == "__main__":
    unittest.main()
