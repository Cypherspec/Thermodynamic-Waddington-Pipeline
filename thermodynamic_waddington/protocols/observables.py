from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..arrays import dot, norm, sub


@dataclass(frozen=True)
class TransitionObservable:
    source: int
    target: int
    displacement: float
    velocity_projection: float
    transverse_velocity: float
    signed_work: float
    activity: float

    def to_dict(self) -> dict[str, float | int]:
        return {
            "source": self.source,
            "target": self.target,
            "displacement": self.displacement,
            "velocity_projection": self.velocity_projection,
            "transverse_velocity": self.transverse_velocity,
            "signed_work": self.signed_work,
            "activity": self.activity,
        }


def transition_observables(points: Sequence[Sequence[float]], velocities: Sequence[Sequence[float]], edges: Sequence[tuple[int, int]], diffusion: Sequence[float]) -> list[TransitionObservable]:
    result = []
    for source, target in edges:
        displacement = sub(points[target], points[source])
        distance = norm(displacement)
        speed = norm(velocities[source])
        if distance == 0.0:
            projection = 0.0
            transverse = speed
        else:
            projection = dot(velocities[source], displacement) / distance
            transverse_sq = max(0.0, speed * speed - projection * projection)
            transverse = transverse_sq ** 0.5
        result.append(TransitionObservable(source, target, distance, projection, transverse, -projection * distance, speed * speed / max(2.0 * diffusion[source], 1e-12)))
    return result


def irreversibility_score(observables: Sequence[TransitionObservable]) -> float:
    if not observables:
        return 0.0
    forward = sum(max(0.0, observable.signed_work) for observable in observables)
    reverse = sum(max(0.0, -observable.signed_work) for observable in observables)
    return (forward - reverse) / max(1.0, forward + reverse)
