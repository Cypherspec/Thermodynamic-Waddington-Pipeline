from __future__ import annotations

"""Deterministic evidence graph for auditable computational biology releases.

The graph makes a scientific claim inspectable as a chain from raw provenance
through transformations, benchmarks, negative results, and preregistered
experiments. It is an accountability layer, not a claim amplifier: missing
measured evidence remains an explicit blocker.
"""

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping


@dataclass(frozen=True)
class EvidenceNode:
    node_id: str
    kind: str
    label: str
    status: str
    path: str | None = None
    provenance: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class EvidenceEdge:
    source: str
    target: str
    relation: str
    status: str = "supported"
    detail: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class EvidenceGraph:
    schema: str = "thermodynamic-waddington/evidence-graph-v1"
    nodes: list[EvidenceNode] = field(default_factory=list)
    edges: list[EvidenceEdge] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    claim_boundary: str = "Computational evidence is not causal biological evidence."

    def add_node(self, node: EvidenceNode) -> None:
        if any(existing.node_id == node.node_id for existing in self.nodes):
            return
        self.nodes.append(node)

    def add_edge(self, edge: EvidenceEdge) -> None:
        if edge.source not in {node.node_id for node in self.nodes}:
            self.blockers.append(f"edge source missing: {edge.source}")
            return
        if edge.target not in {node.node_id for node in self.nodes}:
            self.blockers.append(f"edge target missing: {edge.target}")
            return
        self.edges.append(edge)

    def _reachable(self, source: str, relation: str | None = None) -> set[str]:
        seen = {source}
        frontier = [source]
        while frontier:
            current = frontier.pop()
            for edge in self.edges:
                if edge.source != current or (relation is not None and edge.relation != relation):
                    continue
                if edge.target not in seen:
                    seen.add(edge.target)
                    frontier.append(edge.target)
        return seen

    def validate(self) -> list[str]:
        errors = list(dict.fromkeys(self.blockers))
        ids = [node.node_id for node in self.nodes]
        if len(ids) != len(set(ids)):
            errors.append("duplicate node identifiers")
        for edge in self.edges:
            if edge.source not in ids or edge.target not in ids:
                errors.append(f"dangling edge: {edge.source}->{edge.target}")
        claim_nodes = [node for node in self.nodes if node.kind == "claim"]
        for claim in claim_nodes:
            incoming = [edge for edge in self.edges if edge.target == claim.node_id]
            if not incoming:
                errors.append(f"claim has no evidence edge: {claim.node_id}")
            if claim.status in {"causal", "discovery", "mechanistic"}:
                measured = any(self.nodes_by_id().get(edge.source, EvidenceNode("", "", "", "")).kind == "measured_intervention" for edge in incoming)
                if not measured:
                    errors.append(f"strong claim lacks measured intervention evidence: {claim.node_id}")
        return list(dict.fromkeys(errors))

    def nodes_by_id(self) -> dict[str, EvidenceNode]:
        return {node.node_id: node for node in self.nodes}

    def to_dict(self) -> dict[str, Any]:
        errors = self.validate()
        payload = {
            "schema": self.schema,
            "status": "valid_evidence_graph" if not errors else "blocked_evidence_graph",
            "nodes": [node.to_dict() for node in self.nodes],
            "edges": [edge.to_dict() for edge in self.edges],
            "blockers": errors,
            "warnings": self.warnings,
            "claim_boundary": self.claim_boundary,
            "counts": {"nodes": len(self.nodes), "edges": len(self.edges), "claims": sum(node.kind == "claim" for node in self.nodes)},
        }
        payload["fingerprint"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
        return payload


def _read_json(root: Path, relative: str) -> dict[str, Any] | None:
    path = root / relative
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def _artifact_node(root: Path, relative: str, kind: str = "artifact") -> EvidenceNode:
    path = root / relative
    present = path.is_file() and path.stat().st_size > 0
    digest = None
    if present:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return EvidenceNode(
        node_id=f"artifact:{relative}",
        kind=kind,
        label=relative,
        status="present" if present else "missing",
        path=relative,
        provenance="workspace artifact",
        metadata={"bytes": path.stat().st_size if present else 0, "sha256": digest},
    )


def build_evidence_graph(root: str | Path = ".") -> EvidenceGraph:
    base = Path(root).resolve()
    graph = EvidenceGraph()
    artifacts = {
        "experiments/reproduction_run.json": "reproduction",
        "experiments/evidence_audit.json": "audit",
        "experiments/falsification_report.json": "falsification",
        "experiments/robustness_report.json": "robustness",
        "experiments/research_capsule.json": "capsule",
        "experiments/wetlab_validation_protocol.json": "protocol",
        "experiments/measured_study_manifest.template.json": "study_template",
        "README.md": "documentation",
    }
    for relative, kind in artifacts.items():
        node = _artifact_node(base, relative, kind)
        graph.add_node(node)
    present = {node.node_id for node in graph.nodes if node.status == "present"}

    claim = EvidenceNode(
        "claim:public-heldout-predictive-result",
        "claim",
        "public held-out predictive benchmark",
        "computational_pattern",
        provenance="public benchmark artifact and locked evaluation",
    )
    causal = EvidenceNode(
        "claim:causal-biological-discovery",
        "claim",
        "causal biological discovery",
        "causal",
        provenance="requires measured randomized intervention data",
    )
    graph.add_node(claim)
    graph.add_node(causal)

    def edge(source: str, target: str, relation: str, status: str = "supported", detail: str = "") -> None:
        if source in present or source in {node.node_id for node in graph.nodes}:
            graph.add_edge(EvidenceEdge(source, target, relation, status, detail))

    capsule_id = "artifact:experiments/research_capsule.json"
    benchmark_id = "artifact:experiments/reproduction_run.json"
    audit_id = "artifact:experiments/evidence_audit.json"
    falsification_id = "artifact:experiments/falsification_report.json"
    protocol_id = "artifact:experiments/wetlab_validation_protocol.json"
    template_id = "artifact:experiments/measured_study_manifest.template.json"
    edge(benchmark_id, claim.node_id, "supports", detail="reproducible computational benchmark")
    edge(audit_id, claim.node_id, "audits", detail="evidence audit bounds the claim")
    edge(falsification_id, claim.node_id, "stress-tests", detail="adversarial implementation checks")
    edge(capsule_id, claim.node_id, "packages", detail="review capsule records the benchmark")
    edge(protocol_id, causal.node_id, "specifies-required-evidence", "blocked", "protocol exists but is not measured data")
    edge(template_id, causal.node_id, "defines-input-schema", "blocked", "template is not an executed study")

    if protocol_id not in present:
        graph.blockers.append("wet-lab validation protocol artifact is missing")
    if template_id not in present:
        graph.blockers.append("measured study manifest template is missing")
    graph.blockers.append("causal biological discovery is blocked until measured intervention data are ingested")
    graph.warnings.append("Public LARRY and synthetic results are observational or computational evidence, not wet-lab intervention evidence.")
    return graph


def write_evidence_graph(path: str | Path = "experiments/evidence_graph.json", root: str | Path = ".") -> dict[str, Any]:
    payload = build_evidence_graph(root).to_dict()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
