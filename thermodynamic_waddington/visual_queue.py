from __future__ import annotations

import html
import json
import math
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

from .arrays import mean, percentile


@dataclass(frozen=True)
class QueueColor:
    name: str
    background: str
    foreground: str
    accent: str

    def css(self) -> str:
        return f"--queue-bg:{self.background};--queue-fg:{self.foreground};--queue-accent:{self.accent};"


PHOSPHOR = QueueColor("phosphor", "#07111d", "#e4f6fb", "#39e6bd")
AMBER = QueueColor("amber", "#171006", "#fff4d0", "#ffbf5b")
ICE = QueueColor("ice", "#06121a", "#e8fbff", "#6ed7ff")


@dataclass(frozen=True)
class QueueEvent:
    event_id: str
    stage: str
    state: str
    timestamp: float
    progress: float
    message: str
    payload: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class VisualJob:
    job_id: str
    title: str
    description: str
    total_steps: int
    current_step: int = 0
    state: str = "queued"
    events: list[QueueEvent] = field(default_factory=list)
    started_at: float | None = None
    completed_at: float | None = None
    error: str | None = None

    @property
    def progress(self) -> float:
        if self.total_steps <= 0:
            return 1.0
        return max(0.0, min(1.0, self.current_step / self.total_steps))

    def emit(self, stage: str, state: str, message: str, payload: dict[str, Any] | None = None) -> QueueEvent:
        event = QueueEvent(f"{self.job_id}:{len(self.events):05d}", stage, state, time.time(), self.progress, message, payload or {})
        self.events.append(event)
        self.state = state
        return event

    def step(self, stage: str, message: str, payload: dict[str, Any] | None = None) -> QueueEvent:
        self.current_step = min(self.total_steps, self.current_step + 1)
        return self.emit(stage, "running" if self.current_step < self.total_steps else "complete", message, payload)

    def fail(self, stage: str, error: Exception | str) -> QueueEvent:
        self.error = str(error)
        return self.emit(stage, "failed", self.error, {"error_type": type(error).__name__ if isinstance(error, Exception) else "Error"})

    def to_dict(self) -> dict[str, Any]:
        return {"job_id": self.job_id, "title": self.title, "description": self.description, "total_steps": self.total_steps, "current_step": self.current_step, "progress": self.progress, "state": self.state, "events": [event.to_dict() for event in self.events], "started_at": self.started_at, "completed_at": self.completed_at, "error": self.error}


@dataclass
class VisualQueue:
    """Deterministic, replayable execution queue for visual analytics stages."""
    jobs: list[VisualJob] = field(default_factory=list)
    events: list[QueueEvent] = field(default_factory=list)
    color: QueueColor = PHOSPHOR

    def enqueue(self, job_id: str, title: str, description: str, total_steps: int) -> VisualJob:
        if any(job.job_id == job_id for job in self.jobs):
            raise ValueError(f"job already exists: {job_id}")
        job = VisualJob(job_id, title, description, max(1, total_steps))
        event = job.emit("queue", "queued", "Job entered the visual queue")
        self.jobs.append(job)
        self.events.append(event)
        return job

    def start(self, job: VisualJob) -> QueueEvent:
        if job not in self.jobs:
            raise ValueError("job is not registered with this queue")
        job.started_at = time.time()
        event = job.emit("queue", "running", "Job started")
        self.events.append(event)
        return event

    def record(self, job: VisualJob, stage: str, message: str, payload: dict[str, Any] | None = None) -> QueueEvent:
        event = job.step(stage, message, payload)
        self.events.append(event)
        return event

    def finish(self, job: VisualJob, message: str = "Job completed") -> QueueEvent:
        job.current_step = job.total_steps
        job.completed_at = time.time()
        event = job.emit("queue", "complete", message)
        self.events.append(event)
        return event

    def fail(self, job: VisualJob, stage: str, error: Exception | str) -> QueueEvent:
        event = job.fail(stage, error)
        self.events.append(event)
        return event

    def run(self, job: VisualJob, stages: Sequence[tuple[str, str, Callable[[], dict[str, Any] | None]]]) -> VisualJob:
        self.start(job)
        try:
            for stage, message, operation in stages:
                payload = operation() or {}
                self.record(job, stage, message, payload)
            self.finish(job)
        except Exception as error:
            self.fail(job, "exception", error)
            raise
        return job

    def state_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for job in self.jobs:
            counts[job.state] = counts.get(job.state, 0) + 1
        return counts

    def to_dict(self) -> dict[str, Any]:
        return {"theme": asdict(self.color), "jobs": [job.to_dict() for job in self.jobs], "events": [event.to_dict() for event in self.events], "state_counts": self.state_counts()}

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True))


def _bar(progress: float, width: int = 32) -> str:
    filled = int(max(0.0, min(1.0, progress)) * width)
    return "[" + "█" * filled + "·" * (width - filled) + "]"


def queue_terminal_summary(queue: VisualQueue) -> str:
    lines = ["THERMODYNAMIC WADDINGTON / VISUAL QUEUE", "=" * 54]
    for job in queue.jobs:
        lines.append(f"{job.state:9s} {_bar(job.progress)} {job.progress:6.1%} {job.title}")
        if job.error:
            lines.append(f"  ERROR: {job.error}")
        if job.events:
            lines.append(f"  last: {job.events[-1].message}")
    lines.append("-" * 54)
    lines.append("states: " + ", ".join(f"{key}={value}" for key, value in sorted(queue.state_counts().items())))
    return "\n".join(lines)


def _event_rows(queue: VisualQueue) -> str:
    rows = []
    for event in queue.events[-200:]:
        rows.append(f'<tr data-state="{html.escape(event.state)}"><td>{html.escape(event.stage)}</td><td><span class="state {html.escape(event.state)}">{html.escape(event.state)}</span></td><td>{event.progress:.1%}</td><td>{html.escape(event.message)}</td></tr>')
    return "".join(rows)


def render_queue_html(queue: VisualQueue, output: str | Path, title: str = "Thermodynamic Waddington / Visual Queue") -> None:
    payload = json.dumps(queue.to_dict(), separators=(",", ":"))
    cards = []
    for job in queue.jobs:
        cards.append(f'<article class="job" data-state="{html.escape(job.state)}"><div class="job-head"><div><span class="eyebrow">{html.escape(job.job_id)}</span><h2>{html.escape(job.title)}</h2></div><span class="state {html.escape(job.state)}">{html.escape(job.state)}</span></div><p>{html.escape(job.description)}</p><div class="progress"><i style="width:{job.progress * 100:.2f}%"></i></div><div class="job-foot"><span>{job.current_step}/{job.total_steps} stages</span><span>{job.progress:.1%}</span></div></article>')
    document = f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title><style>
    :root{{{queue.color.css()}}}*{{box-sizing:border-box}}body{{margin:0;background:radial-gradient(circle at 85% -10%,#15354a 0,#07111d 45rem);color:var(--queue-fg);font-family:Inter,system-ui,sans-serif}}main{{max-width:1250px;margin:auto;padding:32px 20px 70px}}.eyebrow{{font:11px ui-monospace,monospace;letter-spacing:.15em;color:var(--queue-accent);text-transform:uppercase}}h1{{font-size:clamp(34px,6vw,72px);letter-spacing:-.07em;line-height:.95;margin:15px 0}}h2{{font-size:18px;margin:7px 0}}p{{color:#8faebb;line-height:1.55}}.hero{{display:flex;justify-content:space-between;gap:20px;align-items:end;margin-bottom:28px}}.hero-copy{{max-width:690px}}.badge{{border:1px solid rgba(150,220,240,.22);border-radius:999px;padding:10px 13px;color:var(--queue-accent);font:11px ui-monospace,monospace;white-space:nowrap}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:12px}}.job{{border:1px solid rgba(150,220,240,.16);border-radius:16px;padding:17px;background:rgba(8,23,37,.72);box-shadow:0 16px 46px #0003}}.job-head,.job-foot{{display:flex;justify-content:space-between;gap:15px;align-items:start}}.job-foot{{font:11px ui-monospace,monospace;color:#7193a0;margin-top:9px}}.state{{font:10px ui-monospace,monospace;text-transform:uppercase;letter-spacing:.12em;border-radius:99px;padding:6px 8px;border:1px solid currentColor}}.state.complete{{color:#39e6bd}}.state.running{{color:#6ed7ff}}.state.queued{{color:#9aabb4}}.state.failed{{color:#ff7777}}.progress{{height:7px;background:#102b3a;border-radius:99px;overflow:hidden}}.progress i{{display:block;height:100%;border-radius:inherit;background:linear-gradient(90deg,#2475ca,var(--queue-accent))}}section{{margin-top:30px}}.panel{{border:1px solid rgba(150,220,240,.16);border-radius:16px;background:rgba(8,23,37,.65);padding:15px;overflow:auto}}table{{border-collapse:collapse;width:100%;font:12px ui-monospace,monospace}}td,th{{padding:9px;text-align:left;border-bottom:1px solid rgba(150,220,240,.10)}}th{{color:var(--queue-accent);font-size:10px;text-transform:uppercase}}button{{background:#102b3a;border:1px solid rgba(150,220,240,.28);color:var(--queue-fg);padding:8px 10px;border-radius:8px;cursor:pointer}}button:hover{{border-color:var(--queue-accent)}}pre{{font:11px ui-monospace,monospace;white-space:pre-wrap;color:#a9d1db;line-height:1.5}}footer{{color:#567887;font:10px ui-monospace,monospace;margin-top:25px}}@media(max-width:700px){{.hero{{display:block}}.badge{{display:inline-block;margin-top:18px}}}}</style></head><body><main><div class="hero"><div class="hero-copy"><span class="eyebrow">live instrument / visual queue</span><h1>{html.escape(title)}</h1><p>Stage-aware orchestration for reproducible landscape inference. Every stage is visible, serializable, and audit-friendly.</p></div><span class="badge">{len(queue.events)} EVENTS · {len(queue.jobs)} JOBS</span></div><div class="grid">{''.join(cards)}</div><section><div style="display:flex;justify-content:space-between;align-items:center"><h2>Event stream</h2><button onclick="filterEvents('all')">all events</button></div><div class="panel"><table><thead><tr><th>stage</th><th>state</th><th>progress</th><th>message</th></tr></thead><tbody id="events">{_event_rows(queue)}</tbody></table></div></section><section><h2>Serialized queue state</h2><div class="panel"><pre id="payload"></pre></div></section><footer>Thermodynamic Waddington · visual queue state is an execution artifact, not a scientific conclusion.</footer></main><script>const DATA={payload};document.getElementById('payload').textContent=JSON.stringify(DATA,null,2);function filterEvents(state){{document.querySelectorAll('#events tr').forEach(row=>{{row.style.display=state==='all'||row.dataset.state===state?'':'none'}})}};</script></body></html>'''
    Path(output).write_text(document)


def queue_from_fit_stages(stages: Sequence[tuple[str, str, dict[str, Any]]], job_id: str = "landscape-fit") -> VisualQueue:
    queue = VisualQueue()
    job = queue.enqueue(job_id, "Effective landscape fit", "A staged, inspectable fit pipeline", len(stages))
    queue.start(job)
    for name, message, payload in stages:
        queue.record(job, name, message, payload)
    queue.finish(job)
    return queue
