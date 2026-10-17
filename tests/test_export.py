import csv
import json
import tempfile
import unittest

from thermodynamic_waddington import make_synthetic_dataset
from thermodynamic_waddington.export import COLUMNS, export_cells, phase_of


class ExportTests(unittest.TestCase):
    def test_writes_table_and_summary(self):
        ds = make_synthetic_dataset(cells=150, genes=16, seed=3)
        labs = sorted(set(ds.labels))
        with tempfile.TemporaryDirectory() as tmp:
            res = export_cells(ds.expression, ds.velocity, ds.labels, [labs[0]], [labs[-1]], tmp,
                               prefix="demo", dataset="synthetic")
            rows = list(csv.DictReader(res.cells_path.open()))
            summary = json.loads(res.summary_path.read_text())
        self.assertEqual(len(rows), 150)
        self.assertEqual(list(rows[0]), list(COLUMNS))
        ids = {r["tw_id"] for r in rows}
        for r in rows:
            q = float(r["tw_committor"])
            self.assertTrue(0.0 <= q <= 1.0)
            self.assertEqual(r["tw_phase"], phase_of(q))
            self.assertTrue(0.0 <= float(r["tw_cyclic_fraction"]) <= 1.0)
            for d in filter(None, r["tw_downstream"].split("|")):
                self.assertIn(d, ids)
                self.assertGreater(float(next(x for x in rows if x["tw_id"] == d)["tw_committor"]), q)
            self.assertIn(r["tw_stage"], r["tw_text"])
        self.assertEqual(summary["format"], "thermodynamic-waddington-export")
        self.assertEqual(sum(p["n_cells"] for p in summary["phases"]), 150)
        self.assertEqual(len(summary["commitment"]["pmf_q"]), len(summary["commitment"]["pmf_kt"]))
        self.assertIn("noise_matched_excess", summary["irreversibility"])

    def test_rejects_mismatched_inputs(self):
        with self.assertRaises(ValueError):
            export_cells([[0.0, 1.0]] * 5, [[0.0, 1.0]] * 4, ["a"] * 5, ["a"], ["a"], tempfile.gettempdir())

    def test_phase_bands(self):
        self.assertEqual(phase_of(0.0), "progenitor")
        self.assertEqual(phase_of(0.5), "transition state")
        self.assertEqual(phase_of(1.0), "committed")


if __name__ == "__main__":
    unittest.main()
