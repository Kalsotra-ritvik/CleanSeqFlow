"""Terminal output and run logs.

* Human-readable, step-by-step messages go to the terminal (stdout; warnings
  and errors to stderr). ``--quiet`` keeps only warnings, errors and the
  final summary; ``--verbose`` adds command lines and repeated explanations.
* Everything shown, plus debug detail, is appended to ``logs/pipeline.log``.
* Structured events (one JSON object per line) go to ``logs/events.jsonl``
  for machine processing.

All methods are thread-safe so concurrent samples do not interleave blocks.
"""

from __future__ import annotations

import json
import os
import sys
import threading
from pathlib import Path
from typing import Any, TextIO

from .explain import STAGE_LABELS, STEPS, stage_header
from .utils import utc_now_iso

_RESET = "\033[0m"
_COLORS = {
    "bold": "\033[1m",
    "green": "\033[32m",
    "yellow": "\033[33m",
    "red": "\033[31m",
    "cyan": "\033[36m",
    "dim": "\033[2m",
}
_STATUS_COLOR = {"COMPLETED": "green", "SKIPPED": "cyan", "FAILED": "red", "RUNNING": "yellow"}


class Console:
    def __init__(
        self,
        quiet: bool = False,
        verbose: bool = False,
        stream: TextIO | None = None,
        err_stream: TextIO | None = None,
        color: bool | None = None,
    ) -> None:
        self.quiet = quiet
        self.verbose = verbose and not quiet
        self.out = stream or sys.stdout
        self.err = err_stream or sys.stderr
        if color is None:
            color = (
                self.out.isatty()
                and os.environ.get("NO_COLOR") is None
                and os.environ.get("TERM") != "dumb"
            )
        self.color = color
        self._lock = threading.Lock()
        self._log: TextIO | None = None
        self._events: TextIO | None = None
        self._explained: set[int] = set()
        self.warnings: list[dict[str, Any]] = []

    # Log files ---------------------------------------------------------------
    def attach_logs(self, text_log: Path, events_log: Path) -> None:
        text_log.parent.mkdir(parents=True, exist_ok=True)
        self._log = open(text_log, "a", encoding="utf-8")  # noqa: SIM115 - closed in close()
        self._events = open(events_log, "a", encoding="utf-8")  # noqa: SIM115

    def close(self) -> None:
        with self._lock:
            for fh in (self._log, self._events):
                if fh is not None:
                    fh.close()
            self._log = self._events = None

    def event(self, kind: str, **fields: Any) -> None:
        with self._lock:
            if self._events is not None:
                record = {"time": utc_now_iso(), "event": kind, **fields}
                self._events.write(json.dumps(record, default=str) + "\n")
                self._events.flush()

    def _write_log(self, text: str) -> None:
        if self._log is not None:
            stamp = utc_now_iso()
            for line in text.splitlines() or [""]:
                self._log.write(f"{stamp} {line}\n")
            self._log.flush()

    # Formatting ----------------------------------------------------------------
    def style(self, text: str, *styles: str) -> str:
        if not self.color or not styles:
            return text
        return "".join(_COLORS[s] for s in styles) + text + _RESET

    def _emit(self, text: str, *, level: str = "info", force: bool = False) -> None:
        """Print *text* (unless quiet) and always write it to the text log."""
        with self._lock:
            self._write_log(text if level == "info" else f"[{level.upper()}] {text}")
            if level in ("warning", "error"):
                print(text, file=self.err, flush=True)
            elif force or not self.quiet:
                print(text, file=self.out, flush=True)

    # Public API ----------------------------------------------------------------
    def info(self, message: str, sample: str | None = None) -> None:
        self._emit(f"[{sample}] {message}" if sample else message)

    def detail(self, message: str) -> None:
        """Shown only with --verbose; always logged."""
        with self._lock:
            self._write_log(message)
            if self.verbose:
                print(self.style(message, "dim"), file=self.out, flush=True)

    def summary(self, message: str) -> None:
        """Shown even with --quiet."""
        self._emit(message, force=True)

    def warning(self, message: str, sample: str | None = None) -> None:
        self.warnings.append({"sample": sample, "message": message})
        prefix = self.style("WARNING", "yellow", "bold")
        self._emit(f"{prefix}{f' [{sample}]' if sample else ''}: {message}", level="warning")
        self.event("warning", sample=sample, message=message)

    def error(self, message: str, hint: str | None = None, sample: str | None = None) -> None:
        prefix = self.style("ERROR", "red", "bold")
        text = f"{prefix}{f' [{sample}]' if sample else ''}: {message}"
        if hint:
            text += f"\n  What to do: {hint}"
        self._emit(text, level="error")
        self.event("error", sample=sample, message=message, hint=hint)

    def step_start(
        self,
        stage: str,
        sample: str | None,
        position: tuple[int, int] | None,
        inputs: list[str],
        tool: str | None,
        version: str | None,
        params: list[str],
        outputs: list[str],
    ) -> None:
        info = stage_header(stage)
        lines = [
            "",
            self.style(
                f"[STEP {info.number}/{len(STEPS)}] {info.title} — {STAGE_LABELS[stage]}", "bold"
            ),
        ]
        first = info.number not in self._explained
        if first or self.verbose:
            lines.append(f"Purpose: {info.purpose}")
            lines.append(f"  {info.details}")
            self._explained.add(info.number)
        if sample is not None and position is not None:
            lines.append(f"Sample: {sample} ({position[0]} of {position[1]})")
        if tool:
            lines.append(f"Tool: {tool} {version or '(version unknown)'}")
        if inputs:
            lines.append("Input:")
            lines.extend(f"  {i}" for i in inputs)
        if params:
            lines.append("Parameters:")
            lines.extend(f"  {p}" for p in params)
        if outputs:
            lines.append("Output:")
            lines.extend(f"  {o}" for o in outputs)
        lines.append("Status: " + self.style("RUNNING", "yellow"))
        self._emit("\n".join(lines))
        self.event("stage_start", stage=stage, sample=sample)

    def step_end(
        self,
        stage: str,
        sample: str | None,
        status: str,
        duration_s: float | None = None,
        notes: list[str] | None = None,
        reason: str | None = None,
    ) -> None:
        label = STAGE_LABELS[stage]
        head = f"{label}{f' [{sample}]' if sample else ''}: " + self.style(
            status, _STATUS_COLOR.get(status, "bold")
        )
        if duration_s is not None:
            head += f" ({duration_s:.1f} s)"
        if reason:
            head += f" — {reason}"
        lines = [head]
        lines.extend(f"  • {n}" for n in (notes or []))
        self._emit("\n".join(lines))
        self.event(
            "stage_end",
            stage=stage,
            sample=sample,
            status=status,
            duration_s=duration_s,
            reason=reason,
        )
