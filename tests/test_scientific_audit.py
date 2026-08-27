from __future__ import annotations

import unittest

from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.model import fit_landscape
from thermodynamic_waddington.protocols.jarzynski_protocol import WorkSample, estimate_work, protocol_diagnostic
from thermodynamic_waddington.synthetic import make_synthetic_dataset


class ScientificAuditTests(unittest.TestCase):
    def test_audit_marks_effective_landscape_and_reports_path_degeneracy(self) -> None:
        dataset = make_synthetic_dataset(cells=36, genes=6, seed=12)
        fit = fit_landscape(
            dataset.expression,
            dataset.velocity,
            FitConfig(neighbors=6, dimensions=4, bootstrap=2, max_paths=8),
            dataset.embedding,
            dataset.labels,
            dataset.cell_ids,
            dataset.lineage_outcomes,
            dataset.metadata,
        )
        audit = fit.diagnostics["scientific_audit"]
        self.assertEqual(audit["status"], "audited_effective_landscape")
        self.assertIn("required_for_physical_claim", audit)
        self.assertEqual(len(fit.diagnostics["jarzynski_diagnostics"]), len(fit.energies))

    def test_missing_velocity_is_explicit(self) -> None:
        dataset = make_synthetic_dataset(cells=24, genes=6, seed=3)
        metadata = dict(dataset.metadata)
        metadata["velocity_status"] = "not_observed"
        fit = fit_landscape(dataset.expression, dataset.velocity, FitConfig(neighbors=5, dimensions=3, bootstrap=2), metadata=metadata)
        self.assertFalse(fit.diagnostics["scientific_audit"]["velocity_observed"])
        self.assertTrue(fit.diagnostics["scientific_audit"]["warnings"])

    def test_protocol_estimator_requires_one_protocol(self) -> None:
        result = estimate_work([WorkSample(1.0, protocol="p"), WorkSample(2.0, protocol="p")])
        self.assertGreater(result.effective_sample_size, 1.0)
        self.assertEqual(result.status, "finite_sample")
        self.assertTrue(protocol_diagnostic([WorkSample(1.0, protocol="p")])["single_protocol"])
        with self.assertRaises(ValueError):
            estimate_work([WorkSample(1.0, protocol="a"), WorkSample(1.0, protocol="b")])


if __name__ == "__main__":
    unittest.main()
