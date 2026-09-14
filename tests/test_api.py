import json
import unittest

from thermodynamic_waddington import AnalysisReport, FitConfig, analyze, make_synthetic_dataset


def _cfg():
    return FitConfig(
        neighbors=12, dimensions=6, seed=5, enable_cycle_decomposition=False,
        bootstrap_replicates=2, entropy_production_bootstrap_replicates=2,
        entropy_production_permutation_replicates=2,
    )


class TestAnalyze(unittest.TestCase):
    def setUp(self):
        self.ds = make_synthetic_dataset(cells=120, genes=16, seed=5)

    def test_basic_report(self):
        rep = analyze(self.ds.expression, self.ds.velocity, config=_cfg())
        self.assertIsInstance(rep, AnalysisReport)
        self.assertEqual(rep.n_cells, 120)
        self.assertGreater(rep.n_edges, 0)

    def test_committor_fields_present(self):
        labs = sorted(set(self.ds.labels))
        rep = analyze(self.ds.expression, self.ds.velocity, labels=self.ds.labels,
                      source_labels=[labs[0]], target_labels=[labs[-1]], config=_cfg())
        self.assertIsNotNone(rep.committor_order)
        self.assertIsNotNone(rep.commitment_barrier_kt)
        self.assertGreaterEqual(rep.commitment_barrier_kt, 0.0)

    def test_serializable(self):
        rep = analyze(self.ds.expression, self.ds.velocity, config=_cfg())
        json.dumps(rep.to_dict())

    def test_partial_labels_warn(self):
        labs = sorted(set(self.ds.labels))
        rep = analyze(self.ds.expression, self.ds.velocity, labels=self.ds.labels,
                      source_labels=[labs[0]], config=_cfg())
        self.assertTrue(rep.warnings)


if __name__ == "__main__":
    unittest.main()
