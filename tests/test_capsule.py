from __future__ import annotations
import unittest
from thermodynamic_waddington.award_capsule import build_capsule

class CapsuleTests(unittest.TestCase):
    def test_capsule_retains_negative_result_and_boundary(self):
        payload = build_capsule().to_dict()
        self.assertEqual(payload["status"], "review_ready_not_discovery")
        statuses = {claim["status"] for claim in payload["claims"]}
        self.assertIn("falsified_or_not_supported", statuses)
        self.assertIn("not_established", statuses)
        self.assertEqual(len(payload["fingerprint"]), 64)

if __name__ == "__main__":
    unittest.main()
