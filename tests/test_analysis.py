import unittest

from thermodynamic_waddington.analysis import report
from thermodynamic_waddington.graph import Edge


class AnalysisReportTests(unittest.TestCase):
    def test_report_runs_on_a_small_chain(self):
        energies = [3.0, 2.0, 1.0, 0.0]
        edges = [Edge(i, i + 1, 1.0, 0.9, work=energies[i + 1] - energies[i]) for i in range(3)]
        velocities = [[1.0, 0.0]] * 4
        result = report(edges, energies, velocities, attractors=[3], draws=32)
        self.assertEqual(result.basin_sizes, {3: 4})
        self.assertEqual(len(result.transition_flux), 1)
        self.assertIn("edge_work", result.uncertainty)
        self.assertEqual(result.uncertainty["edge_work"]["draws"], 32)
        self.assertGreater(result.score, 0.0)


if __name__ == "__main__":
    unittest.main()
