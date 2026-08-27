from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from .arrays import dot, norm, sub


@dataclass(frozen=True)
class LineageScore:
    source: int
    target: int
    velocity_alignment: float
    displacement: float
    fate_probability: float


def score_lineage_links(points: Sequence[Sequence[float]], velocities: Sequence[Sequence[float]], links: Sequence[tuple[int, int]]) -> list[LineageScore]:
    scores: list[LineageScore] = []
    for source, target in links:
        displacement = sub(points[target], points[source])
        velocity = velocities[source]
        denominator = max(1e-12, norm(displacement) * norm(velocity))
        alignment = dot(displacement, velocity) / denominator
        probability = 1.0 / (1.0 + pow(2.718281828, -4.0 * alignment))
        scores.append(LineageScore(source, target, alignment, norm(displacement), probability))
    return scores


def transition_matrix(scores: Sequence[LineageScore], n_cells: int) -> list[list[float]]:
    matrix = [[0.0 for _ in range(n_cells)] for _ in range(n_cells)]
    for score in scores:
        matrix[score.source][score.target] = score.fate_probability
    for row in matrix:
        total = sum(row)
        if total:
            for index in range(len(row)):
                row[index] /= total
    return matrix


def lineage_summary(labels: Sequence[str], outcomes: Sequence[str] | None) -> dict[str, object]:
    if outcomes is None:
        return {"available": False, "n": 0}
    counts: dict[str, int] = {}
    for outcome in outcomes:
        counts[str(outcome)] = counts.get(str(outcome), 0) + 1
    total = sum(counts.values()) or 1
    agreement = sum(1 for label, outcome in zip(labels, outcomes) if label == outcome) / max(1, min(len(labels), len(outcomes)))
    return {"available": True, "n": len(outcomes), "outcome_counts": counts, "label_outcome_agreement": agreement, "outcome_probabilities": {key: value / total for key, value in counts.items()}}
