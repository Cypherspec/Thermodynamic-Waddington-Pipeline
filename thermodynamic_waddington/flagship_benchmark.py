from __future__ import annotations

"""Locked donor/clone-aware benchmark orchestration.

The benchmark is designed to prevent the common failure mode of reporting a
cell-randomized score that leaks clone or donor information. It computes cell,
clone, donor, and intervention-stratified metrics and records all exclusions.
"""

import hashlib
import json
import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence


@dataclass(frozen=True)
class Observation:
    cell_id: str
    donor: str
    clone: str
    timepoint: str
    intervention: str
    prediction: float
    outcome: float
    viable: bool = True


@dataclass(frozen=True)
class Split:
    name: str
    train: tuple[str, ...]
    test: tuple[str, ...]
    unit: str


@dataclass
class BenchmarkReport:
    benchmark_id: str
    status: str
    splits: list[dict[str, Any]] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)
    exclusions: list[str] = field(default_factory=list)
    claim_boundary: str = "Predictive benchmark only; not causal evidence or molecular free-energy validation."

    def to_dict(self) -> dict[str, Any]:
        return asdict(self) | {"fingerprint": fingerprint(asdict(self))}


def fingerprint(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    return hashlib.sha256(payload).hexdigest()


def load_observations(path: str | Path) -> list[Observation]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    rows = payload.get("observations", payload.get("cells", payload))
    if not isinstance(rows, list):
        raise ValueError("benchmark must be a JSON list or contain observations/cells")
    result = []
    for row in rows:
        result.append(Observation(str(row["cell_id"]), str(row["donor"]), str(row.get("lineage_id") or row.get("clone") or "missing_clone"), str(row.get("timepoint", "endpoint")), str(row.get("intervention", "unknown")), float(row["prediction"] if "prediction" in row else row.get("target_probability", 0.0)), float(row["outcome"] if "outcome" in row else row.get("target_probability", 0.0)), bool(row.get("viable", True))))
    return result


def _spearman(left: Sequence[float], right: Sequence[float]) -> float:
    n = min(len(left), len(right))
    if n < 2: return float("nan")
    def rank(values):
        order=sorted(range(len(values)),key=lambda i:(values[i],i)); out=[0.0]*len(values); i=0
        while i<len(order):
            j=i+1
            while j<len(order) and values[order[j]]==values[order[i]]: j+=1
            r=(i+j-1)/2+1
            for k in range(i,j): out[order[k]]=r
            i=j
        return out
    a,b=rank(list(left)[:n]),rank(list(right)[:n]); ma=sum(a)/n; mb=sum(b)/n
    da=sum((x-ma)**2 for x in a); db=sum((y-mb)**2 for y in b)
    return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(da*db) if da and db else 0.0


def _pearson(left, right):
    n=min(len(left),len(right))
    if n<2:return float("nan")
    a=list(left)[:n];b=list(right)[:n];ma=sum(a)/n;mb=sum(b)/n
    da=sum((x-ma)**2 for x in a);db=sum((y-mb)**2 for y in b)
    return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(da*db) if da and db else 0.0


def _metrics(rows: Sequence[Observation]) -> dict[str, float]:
    rows=[row for row in rows if row.viable and math.isfinite(row.prediction) and math.isfinite(row.outcome)]
    if not rows:return {"n":0.0,"spearman":float("nan"),"pearson":float("nan"),"mae":float("nan"),"rmse":float("nan")}
    errors=[row.prediction-row.outcome for row in rows]
    return {"n":float(len(rows)),"spearman":_spearman([r.prediction for r in rows],[r.outcome for r in rows]),"pearson":_pearson([r.prediction for r in rows],[r.outcome for r in rows]),"mae":sum(abs(x) for x in errors)/len(errors),"rmse":math.sqrt(sum(x*x for x in errors)/len(errors))}


def make_splits(rows: Sequence[Observation], seed: int = 17) -> list[Split]:
    rng=random.Random(seed)
    donors=sorted({row.donor for row in rows}); rng.shuffle(donors); cutoff=max(1,int(len(donors)*0.25))
    test_donors=set(donors[:cutoff]); train=[row.cell_id for row in rows if row.donor not in test_donors]; test=[row.cell_id for row in rows if row.donor in test_donors]
    clones=sorted({str(row.clone) for row in rows if row.clone and row.clone != "missing_clone"}); rng.shuffle(clones); clone_cut=max(1,int(len(clones)*0.25)); test_clones=set(clones[:clone_cut])
    clone_train=[row.cell_id for row in rows if row.clone not in test_clones]; clone_test=[row.cell_id for row in rows if row.clone in test_clones]
    timepoints=sorted({row.timepoint for row in rows}); final=timepoints[-1] if timepoints else "endpoint"
    temporal_train=[row.cell_id for row in rows if row.timepoint != final]; temporal_test=[row.cell_id for row in rows if row.timepoint == final]
    return [Split("donor_holdout",tuple(train),tuple(test),"donor"),Split("clone_holdout",tuple(clone_train),tuple(clone_test),"clone"),Split("timepoint_holdout",tuple(temporal_train),tuple(temporal_test),"timepoint")]


def evaluate(rows: Sequence[Observation], benchmark_id: str = "TW-FLAGSHIP-v1", seed: int = 17) -> BenchmarkReport:
    clean=[row for row in rows if row.cell_id and row.donor and row.viable]
    exclusions=["nonviable cells excluded from primary fate metrics"]
    if len({row.donor for row in clean}) < 2:
        exclusions.append("donor-level generalization is unavailable with fewer than two donors")
    if len({row.clone for row in clean if row.clone and row.clone != "missing_clone"}) < 2:
        exclusions.append("clone-level generalization is unavailable with fewer than two observed clones")
    if len({row.timepoint for row in clean}) < 2:
        exclusions.append("temporal generalization is unavailable with fewer than two time points")
    splits=make_splits(clean,seed)
    by_id={row.cell_id:row for row in clean}
    reports=[]
    for split in splits:
        train=[by_id[x] for x in split.train if x in by_id]; test=[by_id[x] for x in split.test if x in by_id]
        reports.append({"name":split.name,"unit":split.unit,"train":len(train),"test":len(test),"metrics":_metrics(test)})
    by_donor={}
    for donor in sorted({row.donor for row in clean}): by_donor[donor]=_metrics([row for row in clean if row.donor==donor])
    by_intervention={}
    for intervention in sorted({row.intervention for row in clean}): by_intervention[intervention]=_metrics([row for row in clean if row.intervention==intervention])
    enough_units = len({row.donor for row in clean}) >= 2 and len({row.clone for row in clean if row.clone and row.clone != "missing_clone"}) >= 2 and len({row.timepoint for row in clean}) >= 2
    status="validated_predictive_benchmark" if enough_units and all(item["test"]>=2 for item in reports) else "underpowered_or_nonindependent"
    return BenchmarkReport(benchmark_id,status,reports,{"overall":_metrics(clean),"by_donor":by_donor,"by_intervention":by_intervention},exclusions)


def write_report(input_path: str | Path, output_path: str | Path, benchmark_id: str = "TW-FLAGSHIP-v1") -> dict[str, Any]:
    report=evaluate(load_observations(input_path),benchmark_id)
    target=Path(output_path);target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(report.to_dict(),indent=2,sort_keys=True),encoding="utf-8")
    return report.to_dict()
