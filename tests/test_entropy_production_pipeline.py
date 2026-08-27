from __future__ import annotations

import json
import unittest

from thermodynamic_waddington.config import FitConfig
from thermodynamic_waddington.model import fit_landscape
from thermodynamic_waddington.synthetic import make_synthetic_dataset


# neighbors defaults to 24 and dimensions to 6 in FitConfig; every fixture
# here uses enough cells/genes to satisfy both without needing to override
# config just to get past validation, matching how the rest of tests/ sizes
# its synthetic datasets.
CELLS = 40
GENES = 8


class EntropyProductionPipelineIntegrationTests(unittest.TestCase):
    """These exercise entropy_production.py exactly as a real caller reaches
    it -- through fit_landscape -- rather than by calling the module
    directly, so a wiring regression (wrong keys passed, dict/dataclass
    mismatch between model.py and entropy_production.py, config fields not
    threaded through) fails here even if the unit tests for the module in
    isolation still pass.
    """

    def test_entropy_production_present_and_well_formed_when_velocity_observed(self) -> None:
        dataset = make_synthetic_dataset(cells=CELLS, genes=GENES, seed=41)
        config = FitConfig(entropy_production_bootstrap_replicates=5, entropy_production_permutation_replicates=5)
        fit = fit_landscape(dataset.expression, dataset.velocity, config=config)
        ep = fit.diagnostics["entropy_production"]
        self.assertIn("entropy_production_rate", ep)
        self.assertGreaterEqual(ep["entropy_production_rate"], 0.0)
        self.assertIn("claim_boundary", ep)
        self.assertEqual(ep["n_cells"], CELLS)

    def test_cross_validation_present_and_numeric_when_velocity_observed(self) -> None:
        dataset = make_synthetic_dataset(cells=CELLS, genes=GENES, seed=42)
        config = FitConfig(entropy_production_bootstrap_replicates=3, entropy_production_permutation_replicates=3)
        fit = fit_landscape(dataset.expression, dataset.velocity, config=config)
        cross = fit.diagnostics["entropy_production_cross_validation"]
        self.assertIsNotNone(cross)
        self.assertIn("graph_based_entropy_production", cross)
        self.assertIn("jarzynski_mean_dissipation", cross)
        self.assertIsInstance(cross["jarzynski_mean_dissipation"], float)

    def test_marked_unavailable_when_velocity_not_observed(self) -> None:
        dataset = make_synthetic_dataset(cells=CELLS, genes=GENES, seed=43)
        zero_velocity = [[0.0] * len(row) for row in dataset.velocity]
        fit = fit_landscape(dataset.expression, zero_velocity, metadata={"velocity_status": "not_observed"})
        ep = fit.diagnostics["entropy_production"]
        self.assertEqual(ep["status"], "unavailable")
        self.assertIsNone(fit.diagnostics["entropy_production_cross_validation"])

    def test_disabled_via_config_flag_even_with_velocity(self) -> None:
        dataset = make_synthetic_dataset(cells=CELLS, genes=GENES, seed=44)
        fit = fit_landscape(dataset.expression, dataset.velocity, config=FitConfig(enable_entropy_production=False))
        ep = fit.diagnostics["entropy_production"]
        self.assertEqual(ep["status"], "unavailable")

    def test_fit_result_is_still_json_serializable_with_entropy_production_included(self) -> None:
        # LandscapeFit.save() round-trips through json.dumps -- confirm the new
        # diagnostics payload (which includes an EntropyProductionConfig-derived
        # dict and possibly None values) doesn't break that contract.
        dataset = make_synthetic_dataset(cells=CELLS, genes=GENES, seed=45)
        config = FitConfig(entropy_production_bootstrap_replicates=2, entropy_production_permutation_replicates=2)
        fit = fit_landscape(dataset.expression, dataset.velocity, config=config)
        payload = json.dumps(fit.to_dict())
        reloaded = json.loads(payload)
        self.assertIn("entropy_production", reloaded["diagnostics"])


if __name__ == "__main__":
    unittest.main()
