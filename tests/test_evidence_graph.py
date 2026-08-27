from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from thermodynamic_waddington.evidence_graph import EvidenceEdge, EvidenceGraph, EvidenceNode, build_evidence_graph


class EvidenceGraphTests(unittest.TestCase):
    def test_causal_claim_requires_measured_intervention(self):
        graph = EvidenceGraph()
        graph.add_node(EvidenceNode("claim", "claim", "causal", "causal"))
        graph.add_node(EvidenceNode("protocol", "protocol", "protocol", "present"))
        graph.add_edge(EvidenceEdge("protocol", "claim", "specifies-required-evidence", "blocked"))
        self.assertTrue(any("measured intervention" in error for error in graph.validate()))

    def test_measured_intervention_can_support_causal_claim(self):
        graph = EvidenceGraph()
        graph.add_node(EvidenceNode("claim", "claim", "causal", "causal"))
        graph.add_node(EvidenceNode("data", "measured_intervention", "data", "present"))
        graph.add_edge(EvidenceEdge("data", "claim", "supports"))
        self.assertEqual(graph.validate(), [])

    def test_workspace_graph_is_explicitly_bounded(self):
        graph = build_evidence_graph(Path(__file__).parents[1])
        payload = graph.to_dict()
        self.assertEqual(payload["status"], "blocked_evidence_graph")
        self.assertTrue(any("measured intervention" in item for item in payload["blockers"]))
        self.assertGreaterEqual(payload["counts"]["nodes"], 8)


if __name__ == "__main__":
    unittest.main()
