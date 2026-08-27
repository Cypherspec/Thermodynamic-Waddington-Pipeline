from __future__ import annotations

import math
import tempfile
from pathlib import Path

import unittest

from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.benchmarking import rank_correlation
from thermodynamic_waddington.causal import intervention_scan, rank_targets
from thermodynamic_waddington.geometry import summarize
from thermodynamic_waddington.observability import fingerprint
from thermodynamic_waddington.topology.persistence import lower_star_pairs
from thermodynamic_waddington.graph import build_knn
from thermodynamic_waddington.model import fit_landscape
from thermodynamic_waddington.synthetic import make_synthetic_dataset
from thermodynamic_waddington.protocols.jarzynski_protocol import WorkSample, estimate_work, protocol_diagnostic
from thermodynamic_waddington.protocols.jarzynski_protocol import WorkSample, estimate_work
from thermodynamic_waddington.protocols.jarzynski_protocol import WorkSample, estimate_work
from thermodynamic_waddington.protocols.jarzynski_protocol import WorkSample, estimate_work, protocol_diagnostic
from thermodynamic_waddington.protocols.jarzynski_protocol import WorkSample, estimate_work, protocol_diagnostic


class CoreTests(unittest.TestCase):

    def test_knn_excludes_self_and_has_expected_degree(self) -> None:
        graph = build_knn([[0.0], [1.0], [2.0], [3.0]], 2)
        self.assertTrue(all(len(indices) == 2 for indices in graph.neighbors))
        self.assertTrue(all(i not in indices for i, indices in enumerate(graph.neighbors)))

    def test_synthetic_generation_is_seeded(self) -> None:
        first = make_synthetic_dataset(cells=30, genes=6, seed=11)
        second = make_synthetic_dataset(cells=30, genes=6, seed=11)
        self.assertEqual(first.expression, second.expression)
        self.assertEqual(first.velocity, second.velocity)

    def test_fit_is_finite_and_has_diagnostics(self) -> None:
        dataset = make_synthetic_dataset(cells=80, genes=8, seed=4)
        config = FitConfig(neighbors=8, dimensions=5, bootstrap_replicates=4, seed=4)
        fit = fit_landscape(dataset.expression, dataset.velocity, config, dataset.embedding, dataset.labels, dataset.cell_ids, dataset.lineage_outcomes)
        self.assertEqual(len(fit.energies), 80)
        self.assertTrue(all(math.isfinite(value) for value in fit.energies))
        self.assertIn("lineage_fate_calibration", fit.diagnostics)
        self.assertGreaterEqual(float(fit.diagnostics["coverage"]), 0)
        self.assertLessEqual(float(fit.diagnostics["coverage"]), 1)

    def test_scientific_audit_is_explicit_about_protocol_identification(self) -> None:
        dataset = make_synthetic_dataset(48, 6, 23)
        config = FitConfig(neighbors=8, dimensions=4, bootstrap_replicates=3, trajectory_steps=3, trajectory_replicates=3)
        fit = fit_landscape(dataset.expression, dataset.velocity, config, dataset.embedding, dataset.labels, dataset.cell_ids, dataset.lineage_outcomes, dataset.metadata)
        audit = fit.diagnostics["scientific_audit"]
        self.assertGreater(audit["path_work"]["count"], 0)
        self.assertFalse(audit["path_protocol"]["protocol_defined"])
        self.assertEqual(audit["diffusion_tensor"]["cells"], 48)
        self.assertEqual(audit["current"]["cells"], 48)

    def test_roundtrip_json(self) -> None:
        dataset = make_synthetic_dataset(cells=30, genes=6, seed=3)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "dataset.json"
            dataset.save(path)
            self.assertTrue(path.exists())
            self.assertGreater(len(path.read_text()), 100)
    def test_scientific_audit_is_explicit_and_jarzynski_module_is_conservative(self) -> None:
        dataset = make_synthetic_dataset(cells=32, genes=6, seed=8)
        fit = fit_landscape(dataset.expression, dataset.velocity, FitConfig(neighbors=6, dimensions=3, bootstrap=2, trajectory_steps=2, trajectory_replicates=2), dataset.embedding)
        audit = fit.diagnostics["scientific_audit"]
        self.assertFalse(audit["claims"]["physical_free_energy_in_kT"])
        self.assertGreater(len(audit["blockers"]), 2)
        samples = [WorkSample(float(index) / 4.0, "forward", f"trajectory-{index}") for index in range(8)]
        estimate = estimate_work(samples)
        self.assertEqual(estimate.status, "diagnostic_only_insufficient_protocol_ensemble")
        self.assertEqual(protocol_diagnostic(samples)["status"], "requires_explicit_protocol_and_more_trajectories")

    def test_advanced_diagnostics_are_deterministic(self) -> None:
        dataset = make_synthetic_dataset(cells=32, genes=6, seed=23)
        fit = fit_landscape(dataset.expression, dataset.velocity, config=FitConfig(neighbors=6, dimensions=3, bootstrap_replicates=2), embedding=dataset.embedding, metadata={"velocity_status": "observed"})
        self.assertEqual(fit.diagnostics["scientific_audit"]["velocity_observed"], True)
        self.assertEqual(len(fit.diagnostics["jarzynski_diagnostics"]), 32)
        self.assertTrue(all(math.isfinite(float(item["free_energy"])) for item in fit.diagnostics["jarzynski_diagnostics"]))

    def test_scientific_audit_and_jarzynski_diagnostics_are_serializable(self) -> None:
        dataset = make_synthetic_dataset(cells=32, genes=6, seed=23)
        fit = fit_landscape(dataset.expression, dataset.velocity, config=FitConfig(neighbors=6, dimensions=3, bootstrap_replicates=2, max_paths=12), embedding=dataset.embedding, metadata={"velocity_status": "observed"})
        audit = fit.diagnostics["scientific_audit"]
        self.assertEqual(audit["velocity_observed"], True)
        self.assertIn("warnings", audit)
        diagnostics = fit.diagnostics["jarzynski_diagnostics"]
        self.assertEqual(len(diagnostics), 32)
        self.assertTrue(all(math.isfinite(float(item["free_energy"])) for item in diagnostics))

        dataset = make_synthetic_dataset(cells=24, genes=6, seed=19)
        fit = fit_landscape(dataset.expression, dataset.velocity, config=FitConfig(neighbors=6, dimensions=3, bootstrap_replicates=2), embedding=dataset.embedding)
        edges = []
        for item in fit.edges:
            from thermodynamic_waddington.graph import Edge
            edges.append(Edge(item["source"], item["target"], item["distance"], item["alignment"]))
        metrics = summarize(fit.embedding, dataset.velocity, edges)
        self.assertTrue(math.isfinite(metrics.mean_speed))
        self.assertEqual(rank_correlation([1.0, 2.0], [1.0, 2.0]), 1.0)
        effects = intervention_scan(edges, fit.energies, damping=0.2)
        self.assertEqual(rank_targets(effects, limit=3), rank_targets(effects, limit=3))
        pairs = lower_star_pairs(fit.energies, [(edge.source, edge.target) for edge in edges])
        self.assertTrue(all(pair.persistence >= 0 for pair in pairs))
        self.assertEqual(fingerprint({"seed": 19}), fingerprint({"seed": 19}))


    def test_scientific_audit_and_protocol_estimator_are_explicit(self) -> None:
        dataset = make_synthetic_dataset(cells=30, genes=6, seed=7)
        fit = fit_landscape(dataset.expression, dataset.velocity, FitConfig(neighbors=6, dimensions=4, bootstrap=2, seed=7), metadata={"velocity_status": "observed"})
        audit = fit.diagnostics["scientific_audit"]
        self.assertEqual(audit["status"], "audited_effective_landscape")
        self.assertIn("controlled repeated trajectories or an explicitly justified path measure", audit["required_for_physical_claim"])
        estimate = estimate_work([WorkSample(1.0, "pulse", "a"), WorkSample(1.5, "pulse", "b")])
        self.assertTrue(math.isfinite(estimate.free_energy))
        self.assertGreaterEqual(estimate.effective_sample_size, 1.0)
        diagnostic = protocol_diagnostic([WorkSample(1.0, "pulse", "a")])
        self.assertFalse(diagnostic["reverse_samples"])
        with self.assertRaises(ValueError):
            estimate_work([WorkSample(1.0, "a"), WorkSample(1.0, "b")])
    def test_protocol_aware_jarzynski_estimator_is_stable_and_audited(self) -> None:
        samples = [WorkSample(1.0, "pull-a", "t0"), WorkSample(1.5, "pull-a", "t1"), WorkSample(2.0, "pull-a", "t2"), WorkSample(0.8, "pull-a", "r0", "reverse")]
        estimate = estimate_work(samples, temperature=1.0)
        self.assertTrue(math.isfinite(estimate.free_energy))
        self.assertGreater(estimate.effective_sample_size, 0)
        self.assertEqual(protocol_diagnostic(samples)["status"], "ready")

    def test_scientific_audit_flags_unobserved_velocity(self) -> None:
        dataset = make_synthetic_dataset(cells=24, genes=6, seed=8)
        fit = fit_landscape(dataset.expression, dataset.velocity, FitConfig(neighbors=6, dimensions=3, bootstrap=2), metadata={"velocity_status": "not_observed"})
        audit = fit.diagnostics["scientific_audit"]
        self.assertFalse(audit["velocity_observed"])
        self.assertTrue(audit["warnings"])

    def test_scientific_audit_exposes_unvalidated_claims(self) -> None:
        dataset = make_synthetic_dataset(cells=24, genes=6, seed=7)
        fit = fit_landscape(dataset.expression, dataset.velocity, FitConfig(neighbors=6, dimensions=3, bootstrap=2), metadata={"velocity_status": "observed"})
        audit = fit.diagnostics["scientific_audit"]
        self.assertEqual(audit["protocol_declared"], False)
        self.assertIn("protocol-defined work ensemble", audit["blockers"])
        self.assertIn("physical_free_energy_in_kT", audit["claims"])

    def test_protocol_aware_jarzynski_estimator_is_stable(self) -> None:
        samples = [WorkSample(value, "protocol-a", f"path-{index}") for index, value in enumerate([1.0, 1.2, 0.8, 1.1, 0.9])]
        estimate = estimate_work(samples, temperature=1.0)
        self.assertTrue(math.isfinite(estimate.free_energy))
        self.assertGreater(estimate.effective_sample_size, 1.0)
        self.assertEqual(estimate.protocols, ("protocol-a",))

    def test_scientific_scope_audit_is_explicit(self) -> None:
        dataset = make_synthetic_dataset(cells=24, genes=6, seed=22)
        fit = fit_landscape(dataset.expression, dataset.velocity, FitConfig(neighbors=5, dimensions=3, bootstrap=2), metadata={"velocity_status": "observed"})
        audit = fit.diagnostics["scientific_audit"]
        self.assertEqual(audit["status"], "provisional")
        self.assertFalse(audit["claims"]["physical_free_energy_in_kT"])
        self.assertGreater(len(audit["blockers"]), 0)

    def test_protocol_jarzynski_estimator_reports_effective_sample_size(self) -> None:
        estimate = estimate_work([WorkSample(1.0, "pull", "a"), WorkSample(1.2, "pull", "b"), WorkSample(0.9, "pull", "c")])
        self.assertEqual(estimate.samples, 3)
        self.assertGreater(estimate.effective_sample_size, 1.0)
        self.assertEqual(estimate.protocols, ("pull",))
