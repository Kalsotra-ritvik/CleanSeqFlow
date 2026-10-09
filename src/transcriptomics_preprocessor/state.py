"""Run journal: the status of every stage for every sample, persisted atomically.

A stage is only considered reusable (skipped on ``--resume``) when *all* of
the following hold:

1. its recorded status is ``COMPLETED``;
2. its recorded fingerprint equals the fingerprint computed now from the
   stage's parameters, the tool version and the identity (size and
   modification time) of its input files;
3. every recorded output file still exists with the recorded size and
   modification time.

Because fingerprints include the identity of input files, regenerating an
upstream output automatically invalidates everything downstream.
"""

from __future__ import annotations

import contextlib
import fcntl
import json
import os
import threading
from collections.abc import Iterator
from enum import Enum
from pathlib import Path
from typing import Any

from . import __version__
from .errors import OutputDirectoryError, StateError
from .utils import atomic_write_json, file_signature, utc_now_iso

SCHEMA_VERSION = 1


class Status(str, Enum):
    NOT_STARTED = "NOT_STARTED"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    SKIPPED = "SKIPPED"


#: Per-sample stages, in execution order.
SAMPLE_STAGES: tuple[str, ...] = (
    "download",
    "convert",
    "validate_raw",
    "fastqc_raw",
    "fastp",
    "validate_clean",
    "fastqc_clean",
    "publish",
)
#: Run-level stages that combine all samples.
AGGREGATE_STAGES: tuple[str, ...] = ("multiqc_raw", "multiqc_clean", "report")


def empty_state() -> dict[str, Any]:
    now = utc_now_iso()
    return {
        "schema_version": SCHEMA_VERSION,
        "tool": "transcriptomics-preprocessor",
        "tool_version": __version__,
        "created": now,
        "updated": now,
        "samples": {},
        "aggregate": {},
        "runs": [],
    }


def _stage_record() -> dict[str, Any]:
    return {"status": Status.NOT_STARTED.value}


class RunState:
    """Thread-safe wrapper around the JSON journal file."""

    def __init__(self, path: Path, data: dict[str, Any]) -> None:
        self.path = path
        self.data = data
        self._lock = threading.RLock()

    # Loading / saving -------------------------------------------------------
    @classmethod
    def load(cls, path: Path) -> RunState:
        try:
            text = path.read_text(encoding="utf-8")
            data = json.loads(text)
        except FileNotFoundError:
            raise
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise StateError(
                f"the run journal {path} is unreadable or corrupt ({exc})",
                hint="A backup may exist next to it (state.json.bak). Otherwise restart with "
                "--overwrite; downloaded SRA files are reused by prefetch.",
            ) from exc
        if not isinstance(data, dict) or data.get("schema_version") != SCHEMA_VERSION:
            raise StateError(
                f"the run journal {path} has an unsupported format "
                f"(schema_version={data.get('schema_version') if isinstance(data, dict) else '?'})",
                hint="Restart with --overwrite, or use the program version that created it.",
            )
        for key in ("samples", "aggregate", "runs"):
            if not isinstance(data.get(key), (dict if key != "runs" else list)):
                raise StateError(f"the run journal {path} is missing section '{key}'")
        return cls(path, data)

    @classmethod
    def new(cls, path: Path) -> RunState:
        return cls(path, empty_state())

    def save(self) -> None:
        with self._lock:
            self.data["updated"] = utc_now_iso()
            if self.path.exists():
                with contextlib.suppress(OSError):
                    atomic_write_json(
                        self.path.with_suffix(".json.bak"),
                        json.loads(self.path.read_text(encoding="utf-8")),
                    )
            atomic_write_json(self.path, self.data)

    # Access helpers ---------------------------------------------------------
    def sample_entry(self, name: str) -> dict[str, Any]:
        with self._lock:
            samples = self.data["samples"]
            entry = samples.setdefault(name, {"stages": {}})
            entry.setdefault("stages", {})
            return entry  # type: ignore[no-any-return]

    def stage(self, sample: str | None, stage: str) -> dict[str, Any]:
        with self._lock:
            if sample is None:
                return self.data["aggregate"].setdefault(stage, _stage_record())  # type: ignore[no-any-return]
            return self.sample_entry(sample)["stages"].setdefault(stage, _stage_record())  # type: ignore[no-any-return]

    def status(self, sample: str | None, stage: str) -> Status:
        return Status(self.stage(sample, stage).get("status", Status.NOT_STARTED.value))

    def update(self, sample: str | None, stage: str, **fields: Any) -> None:
        with self._lock:
            rec = self.stage(sample, stage)
            for k, v in fields.items():
                if isinstance(v, Status):
                    v = v.value
                rec[k] = v
            self.save()

    def reset(self, sample: str | None, stage: str) -> None:
        with self._lock:
            target = (
                self.data["aggregate"] if sample is None else self.sample_entry(sample)["stages"]
            )
            target[stage] = _stage_record()

    def mark_interrupted(self) -> None:
        """Turn any RUNNING records into FAILED (used after an interruption)."""
        with self._lock:
            records = [
                rec for entry in self.data["samples"].values() for rec in entry["stages"].values()
            ] + list(self.data["aggregate"].values())
            for rec in records:
                if rec.get("status") == Status.RUNNING.value:
                    rec["status"] = Status.FAILED.value
                    rec["message"] = "interrupted before completion"
                    rec["finished"] = utc_now_iso()
            self.save()


def outputs_record(paths: dict[str, Path], root: Path) -> dict[str, dict[str, Any]]:
    """Describe output files (relative path + signature) for the journal."""
    out: dict[str, dict[str, Any]] = {}
    for label, p in paths.items():
        rel = os.path.relpath(p, root)
        out[label] = {"path": rel, **file_signature(p)}
    return out


def outputs_intact(record: dict[str, Any], root: Path) -> tuple[bool, str]:
    """Check that the outputs recorded for a stage still exist unchanged."""
    outputs = record.get("outputs")
    if not isinstance(outputs, dict):
        return False, "no outputs recorded"
    for label, info in outputs.items():
        try:
            p = root / info["path"]
            sig = file_signature(p)
        except FileNotFoundError:
            return False, f"output '{label}' is missing ({info.get('path')})"
        except (KeyError, TypeError):
            return False, f"output record '{label}' is malformed"
        if sig["size"] != info.get("size") or sig["mtime_ns"] != info.get("mtime_ns"):
            return False, f"output '{label}' changed since it was recorded ({info.get('path')})"
    return True, "outputs intact"


def can_reuse(record: dict[str, Any], fingerprint: str, root: Path) -> tuple[bool, str]:
    """Decide whether a previously completed stage can be skipped."""
    if record.get("status") != Status.COMPLETED.value:
        return False, f"previous status {record.get('status', Status.NOT_STARTED.value)}"
    if record.get("fingerprint") != fingerprint:
        return False, "parameters, tool version or inputs changed since the last run"
    return outputs_intact(record, root)


@contextlib.contextmanager
def exclusive_lock(lock_path: Path) -> Iterator[None]:
    """Hold an exclusive advisory lock so two runs cannot share one output directory."""
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    fh = open(lock_path, "a+")  # noqa: SIM115 - closed in finally
    try:
        try:
            fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise OutputDirectoryError(
                f"another run is already using this output directory (lock: {lock_path})",
                hint="Wait for it to finish, or choose a different --outdir.",
            ) from exc
        fh.seek(0)
        fh.truncate()
        fh.write(f"pid {os.getpid()} since {utc_now_iso()}\n")
        fh.flush()
        yield
    finally:
        with contextlib.suppress(OSError):
            fcntl.flock(fh.fileno(), fcntl.LOCK_UN)
        fh.close()
