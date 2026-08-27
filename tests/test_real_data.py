from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from thermodynamic_waddington.real_data import catalog_payload


class RealDataWorkflowTests(unittest.TestCase):
    def test_catalog_is_explicit_and_has_three_benchmarks(self) -> None:
        payload = catalog_payload()
        self.assertEqual(len(payload["datasets"]), 3)
        self.assertEqual(payload["scientific_policy"], "Real-data downloads are opt-in; no dataset is silently fetched during fitting.")
        self.assertTrue(all(item["source_url"].startswith(("https://", "http://")) for item in payload["datasets"]))

    def test_catalog_roundtrip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "catalog.json"
            path.write_text(json.dumps(catalog_payload()))
            self.assertEqual(len(json.loads(path.read_text())["datasets"]), 3)
