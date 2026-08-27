from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

from .arrays import mean, percentile
from .visual_queue import QueueEvent, VisualQueue, VisualJob


@dataclass(frozen=True)
class StreamCheckpoint:
    checkpoint_id: str
    source: str
    row_start: int
    row_end: int
    emitted_at: float
    digest: str
    statistics: dict[str, float]

    def to_dict(self) -> dict[str, Any]:
        return {"checkpoint_id": self.checkpoint_id, "source": self.source, "row_start": self.row_start, "row_end": self.row_end, "emitted_at": self.emitted_at, "digest": self.digest, "statistics": self.statistics}


@dataclass
class IncrementalWindow:
    source: str
    window_size: int
    row_offset: int = 0
    values: list[float] = field(default_factory=list)
    checkpoints: list[StreamCheckpoint] = field(default_factory=list)

    def append(self, values: Sequence[float]) -> StreamCheckpoint | None:
        self.values.extend(float(value) for value in values)
        if len(self.values) < self.window_size:
            return None
        window = self.values[:self.window_size]
        self.values = self.values[self.window_size:]
        start = self.row_offset
        self.row_offset += len(window)
        checkpoint = StreamCheckpoint(f"{self.source}:{len(self.checkpoints):06d}", self.source, start, self.row_offset, time.time(), _digest(window), {"mean": mean(window), "p50": percentile(window, 0.5), "p95": percentile(window, 0.95), "n": float(len(window))})
        self.checkpoints.append(checkpoint)
        return checkpoint

    def flush(self) -> StreamCheckpoint | None:
        if not self.values:
            return None
        window = self.values
        self.values = []
        start = self.row_offset
        self.row_offset += len(window)
        checkpoint = StreamCheckpoint(f"{self.source}:{len(self.checkpoints):06d}", self.source, start, self.row_offset, time.time(), _digest(window), {"mean": mean(window), "p50": percentile(window, 0.5), "p95": percentile(window, 0.95), "n": float(len(window))})
        self.checkpoints.append(checkpoint)
        return checkpoint

    def to_dict(self) -> dict[str, Any]:
        return {"source": self.source, "window_size": self.window_size, "row_offset": self.row_offset, "buffered": self.values, "checkpoints": [checkpoint.to_dict() for checkpoint in self.checkpoints]}


def _digest(values: Sequence[float]) -> str:
    state = 2166136261
    for value in values:
        for byte in repr(round(float(value), 12)).encode():
            state ^= byte
            state = (state * 16777619) & 0xFFFFFFFF
    return f"fnv1a:{state:08x}"


def stream_values(values: Iterable[float], window_size: int, source: str = "stream") -> Iterator[StreamCheckpoint]:
    window = IncrementalWindow(source, max(1, window_size))
    for value in values:
        checkpoint = window.append([value])
        if checkpoint is not None:
            yield checkpoint
    checkpoint = window.flush()
    if checkpoint is not None:
        yield checkpoint


def attach_stream_to_queue(queue: VisualQueue, job: VisualJob, checkpoints: Sequence[StreamCheckpoint], stage: str = "stream") -> None:
    for checkpoint in checkpoints:
        queue.record(job, stage, f"Committed checkpoint {checkpoint.checkpoint_id}", checkpoint.to_dict())


def save_checkpoints(checkpoints: Sequence[StreamCheckpoint], path: str | Path) -> None:
    Path(path).write_text(json.dumps([checkpoint.to_dict() for checkpoint in checkpoints], indent=2, sort_keys=True))
