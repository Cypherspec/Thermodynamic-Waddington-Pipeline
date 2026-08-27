from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.model import fit_landscape
from thermodynamic_waddington.streaming_queue import IncrementalWindow, stream_values
from thermodynamic_waddington.synthetic import make_synthetic_dataset
from thermodynamic_waddington.visual_diagnostics import diagnose_fit, diagnostics_payload
from thermodynamic_waddington.visual_queue import queue_from_fit_stages


class VisualWorkflowTests(unittest.TestCase):
    def test_incremental_window_commits_complete_and_tail_windows(self) -> None:
        window = IncrementalWindow("cells", 3)
        self.assertIsNone(window.append([1.0, 2.0]))
        checkpoint = window.append([3.0, 4.0])
        self.assertIsNotNone(checkpoint)
        tail = window.flush()
        self.assertIsNotNone(tail)
        self.assertEqual([checkpoint.row_start, checkpoint.row_end], [0, 3])
        self.assertEqual([tail.row_start, tail.row_end], [3, 4])

    def test_stream_values_is_deterministic_in_partitioning(self) -> None:
        checkpoints = list(stream_values(range(7), 3, "demo"))
        self.assertEqual([(item.row_start, item.row_end) for item in checkpoints], [(0, 3), (3, 6), (6, 7)])
        self.assertTrue(all(item.digest.startswith("fnv1a:") for item in checkpoints))

    def test_visual_diagnostics_are_serializable(self) -> None:
        dataset = make_synthetic_dataset(cells=32, genes=6, seed=17)
        fit = fit_landscape(dataset.expression, dataset.velocity, FitConfig(neighbors=8, dimensions=3, bootstrap_replicates=3))
        payload = diagnostics_payload(fit)
        json.dumps(payload)
        self.assertEqual(len(diagnose_fit(fit)), 5)
        self.assertIn(payload["overall"], {"good", "warning", "critical"})

    def test_visual_queue_completes(self) -> None:
        queue = queue_from_fit_stages([("a", "A", {"n": 1}), ("b", "B", {"n": 2})])
        self.assertEqual(queue.jobs[0].state, "complete")
        self.assertEqual(queue.jobs[0].progress, 1.0)
        self.assertEqual(len(queue.events), 5)


if __name__ == "__main__":
    unittest.main()
