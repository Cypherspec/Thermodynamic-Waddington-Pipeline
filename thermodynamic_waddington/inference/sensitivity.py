from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..arrays import mean


@dataclass(frozen=True)
class SensitivityProfile:
    parameter: str
    values: tuple[float, ...]
    output_means: tuple[float, ...]
    output_spreads: tuple[float, ...]
    normalized_slope: float

    def to_dict(self) -> dict[str, object]:
        return {"parameter": self.parameter, "values": list(self.values), "output_means": list(self.output_means), "output_spreads": list(self.output_spreads), "normalized_slope": self.normalized_slope}


def profile(parameter: str, values: Sequence[float], outputs: Sequence[Sequence[float]]) -> SensitivityProfile:
    means = tuple(mean(output) if output else 0.0 for output in outputs)
    spreads = tuple((max(output) - min(output)) if output else 0.0 for output in outputs)
    if len(means) > 1 and values[-1] != values[0]:
        slope = (means[-1] - means[0]) / (values[-1] - values[0])
    else:
        slope = 0.0
    return SensitivityProfile(parameter, tuple(values), means, spreads, slope)
