"""Advanced-analysis subpackage.

Note: this subpackage previously contained a large amount of templated
placeholder code (multiple ~788-line files that were structurally identical
copies of each other with only a domain label changed) which has been
removed. Only genuinely distinct, functioning modules remain: `core`
(shared numeric/audit primitives) and `fate_control` (intervention-scoring
engine used by discovery_engine.py). If you are looking for
`AdvancedRun`/`run_advanced` or the fate-graph/energy-landscape "engines"
that used to be exported here, those were part of the removed placeholder
code and were never functionally distinct from each other.
"""
from .core import ObservationContract, EvidenceLedger, OnlineMoments, ReproducibleEngine
from .fate_control import FateControlConfig, FateControlEngine, Intervention, run_fate_control

__all__ = [
    'ObservationContract', 'EvidenceLedger', 'OnlineMoments', 'ReproducibleEngine',
    'FateControlConfig', 'FateControlEngine', 'Intervention', 'run_fate_control',
]
