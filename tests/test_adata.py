import importlib.util
import unittest

import numpy as np

from thermodynamic_waddington import AnalysisReport, FitConfig, analyze_adata

_HAS_ANNDATA = importlib.util.find_spec("anndata") is not None


def _cfg():
    return FitConfig(
        neighbors=8, dimensions=4, seed=3, enable_cycle_decomposition=False,
        bootstrap_replicates=2, entropy_production_bootstrap_replicates=2,
        entropy_production_permutation_replicates=3,
    )


@unittest.skipUnless(_HAS_ANNDATA, "anndata not installed")
class TestAnalyzeAdata(unittest.TestCase):
    def _adata(self):
        import anndata as ad
        import pandas as pd
        rng = np.random.default_rng(0)
        n, g = 60, 10
        spliced = rng.gamma(2.0, 1.0, size=(n, g))
        unspliced = spliced + rng.normal(0.0, 0.3, size=(n, g))
        labels = ["A"] * 30 + ["B"] * 30
        obs = pd.DataFrame({"clusters": labels})
        return ad.AnnData(X=spliced, layers={"spliced": spliced, "unspliced": unspliced}, obs=obs)

    def test_runs_and_writes_back(self):
        a = self._adata()
        rep = analyze_adata(a, label_key="clusters", source=["A"], target=["B"],
                            config=_cfg(), n_top_genes=8)
        self.assertIsInstance(rep, AnalysisReport)
        self.assertIn("tw_energy", a.obs.columns)
        self.assertIn("tw_committor", a.obs.columns)
        self.assertIn("tw", a.uns)
        self.assertEqual(len(a.obs["tw_energy"]), a.n_obs)

    def test_proxy_velocity_from_layers(self):
        a = self._adata()  # no explicit velocity layer -> proxy from spliced/unspliced
        rep = analyze_adata(a, config=_cfg(), n_top_genes=8)
        self.assertEqual(rep.n_cells, a.n_obs)

    def test_missing_velocity_raises(self):
        import anndata as ad
        rng = np.random.default_rng(1)
        a = ad.AnnData(X=rng.normal(size=(40, 6)))  # no layers at all
        with self.assertRaises(ValueError):
            analyze_adata(a, config=_cfg())


if __name__ == "__main__":
    unittest.main()
