from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from thermodynamic_waddington.flagship_benchmark import Observation
from thermodynamic_waddington.replication import ReplicationStudy, evaluate_replication, load_studies


class ReplicationTests(unittest.TestCase):
    def _study(self, study_id: str, donor_prefix: str, offset: float = 0.0) -> ReplicationStudy:
        rows = tuple(
            Observation(f"{study_id}-{i}", f"{donor_prefix}-{i % 2}", f"clone-{i}", "day-2", "control", (i + offset) / 10.0, (i + offset) / 10.0, True)
            for i in range(10)
        )
        return ReplicationStudy(study_id, f"source://{study_id}", rows, "assay-v1", f"pre-{study_id}")

    def test_replication_uses_studies_not_cells(self) -> None:
        report = evaluate_replication([self._study("a", "da"), self._study("b", "db")])
        self.assertEqual(report.status, "replicated_predictive_result")
        self.assertFalse(report.meta_analysis["cell_pooling_used"])
        self.assertEqual(report.independent_studies, 2)

    def test_duplicate_study_content_does_not_count_twice(self) -> None:
        first = self._study("a", "d")
        second = ReplicationStudy("b", first.provenance, first.observations, first.assay_version, first.preregistration_hash)
        report = evaluate_replication([first, second])
        self.assertEqual(report.status, "independent_replication_not_established")
        self.assertEqual(report.independent_studies, 1)
        self.assertTrue(report.duplicate_fingerprints)

    def test_json_loader(self) -> None:
        payload = {"studies": [{"study_id": "a", "provenance": "x", "assay_version": "v1", "observations": [{"cell_id": "c", "donor": "d", "clone": "l", "prediction": 0.2, "outcome": 0.3}]}]}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "studies.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            self.assertEqual(load_studies(path)[0].study_id, "a")


if __name__ == "__main__":
    unittest.main()
