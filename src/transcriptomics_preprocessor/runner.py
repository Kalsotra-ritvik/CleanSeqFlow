"""Safe execution of external programs.

* Commands are always argument lists; ``shell=True`` is never used.
* Every command, its working directory, exit code and duration are appended to
  a per-step log file together with the program's own output.
* Children run in their own process group so that an interruption terminates
  the program *and* anything it spawned (e.g. the JVM started by FastQC).
* A process-wide stop event lets signal handlers cancel running commands.
"""

from __future__ import annotations

import contextlib
import os
import shlex
import signal
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path

from .errors import InterruptedRun, ToolExecutionError
from .utils import utc_now_iso

#: Set by the signal handlers; checked by running commands and the scheduler.
STOP_EVENT = threading.Event()

_ACTIVE: set[subprocess.Popen[bytes]] = set()
_ACTIVE_LOCK = threading.Lock()

_POLL_SECONDS = 0.2
_TERMINATE_GRACE_SECONDS = 10.0


@dataclass(frozen=True)
class CommandResult:
    """Outcome of an external command."""

    args: tuple[str, ...]
    returncode: int
    duration_s: float
    log_path: Path | None
    output_tail: str
    timed_out: bool = False

    @property
    def command_line(self) -> str:
        return shlex.join(self.args)


def _tail(path: Path | None, max_bytes: int = 4000) -> str:
    if path is None or not path.exists():
        return ""
    with open(path, "rb") as fh:
        fh.seek(0, os.SEEK_END)
        size = fh.tell()
        fh.seek(max(0, size - max_bytes))
        return fh.read().decode("utf-8", errors="replace")


def _kill_group(proc: subprocess.Popen[bytes], sig: int) -> None:
    with contextlib.suppress(ProcessLookupError, PermissionError):
        os.killpg(proc.pid, sig)


def _terminate(proc: subprocess.Popen[bytes]) -> None:
    """Terminate a child process group: SIGTERM, then SIGKILL after a grace period."""
    if proc.poll() is not None:
        return
    _kill_group(proc, signal.SIGTERM)
    try:
        proc.wait(timeout=_TERMINATE_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        _kill_group(proc, signal.SIGKILL)
        proc.wait()


def terminate_all() -> None:
    """Terminate every command currently started by :func:`run_command`."""
    with _ACTIVE_LOCK:
        procs = list(_ACTIVE)
    for proc in procs:
        _terminate(proc)


def run_command(
    args: list[str],
    *,
    log_path: Path | None,
    timeout: float | None = None,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    stdout_path: Path | None = None,
    check: bool = True,
    description: str = "",
) -> CommandResult:
    """Run *args* and wait for it to finish.

    Parameters
    ----------
    args:
        Program and arguments. Each element is passed verbatim (no shell).
    log_path:
        File to which a header plus the program's stdout/stderr are appended.
        ``None`` discards output (used only for trivial version probes).
    timeout:
        Seconds before the command is terminated. ``None`` means no limit.
    stdout_path:
        If given, the program's stdout (binary data) is written to this file
        and only stderr goes to the log.
    check:
        Raise :class:`ToolExecutionError` on a non-zero exit code or timeout.

    Raises
    ------
    InterruptedRun
        If :data:`STOP_EVENT` is set while the command runs.
    ToolExecutionError
        If the program cannot be started, or (with ``check``) fails.
    """
    if not args or not all(isinstance(a, str) for a in args):
        raise TypeError("args must be a non-empty list of strings")
    if STOP_EVENT.is_set():
        raise InterruptedRun("run interrupted before starting: " + (description or args[0]))

    log_fh = None
    out_fh = None
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_fh = open(log_path, "ab")  # noqa: SIM115 - closed in finally
        header = (
            f"\n### {utc_now_iso()} {description}\n### cwd: {cwd or Path.cwd()}\n"
            f"### command: {shlex.join(args)}\n"
        )
        log_fh.write(header.encode())
        log_fh.flush()
    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        out_fh = open(stdout_path, "wb")  # noqa: SIM115 - closed in finally

    stdout_target = out_fh if out_fh is not None else (log_fh or subprocess.DEVNULL)
    stderr_target = log_fh if log_fh is not None else subprocess.DEVNULL
    start = time.monotonic()
    timed_out = False
    returncode: int | None = None
    try:
        try:
            proc = subprocess.Popen(
                args,
                cwd=cwd,
                env=env,
                stdin=subprocess.DEVNULL,
                stdout=stdout_target,
                stderr=stderr_target,
                start_new_session=True,
            )
        except FileNotFoundError as exc:
            raise ToolExecutionError(
                f"program not found: {args[0]}",
                hint="Run 'transcriptomics-preprocess doctor' to check the installation.",
                log_path=str(log_path) if log_path else None,
            ) from exc
        except PermissionError as exc:
            raise ToolExecutionError(
                f"program is not executable: {args[0]} ({exc.strerror})",
                hint="Check file permissions or re-run ./install.sh --repair.",
            ) from exc

        with _ACTIVE_LOCK:
            _ACTIVE.add(proc)
        try:
            while True:
                try:
                    proc.wait(timeout=_POLL_SECONDS)
                    break
                except subprocess.TimeoutExpired:
                    pass
                if STOP_EVENT.is_set():
                    _terminate(proc)
                    raise InterruptedRun(
                        "run interrupted while running: " + (description or args[0])
                    )
                if timeout is not None and time.monotonic() - start > timeout:
                    timed_out = True
                    _terminate(proc)
                    break
        finally:
            with _ACTIVE_LOCK:
                _ACTIVE.discard(proc)
        returncode = proc.returncode
    finally:
        duration = time.monotonic() - start
        if log_fh is not None:
            log_fh.write(
                f"### exit code: {returncode if returncode is not None else 'n/a'} "
                f"(duration {duration:.1f} s{', TIMED OUT' if timed_out else ''})\n".encode()
            )
            log_fh.close()
        if out_fh is not None:
            out_fh.close()

    assert returncode is not None
    result = CommandResult(
        args=tuple(args),
        returncode=returncode,
        duration_s=duration,
        log_path=log_path,
        output_tail=_tail(log_path),
        timed_out=timed_out,
    )
    if check and (timed_out or returncode != 0):
        reason = (
            f"timed out after {timeout:.0f} s"
            if timed_out
            else f"exited with code {returncode}"
            + (f" (killed by signal {-returncode})" if returncode < 0 else "")
        )
        raise ToolExecutionError(
            f"{description or args[0]} {reason}",
            hint=f"See the log for details: {log_path}" if log_path else None,
            returncode=returncode,
            log_path=str(log_path) if log_path else None,
            output_tail=result.output_tail,
        )
    return result


def probe_output(args: list[str], timeout: float = 60.0) -> tuple[int, str]:
    """Run a quick command (e.g. ``--version``) and return (exit code, combined output)."""
    try:
        cp = subprocess.run(
            args,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    except FileNotFoundError:
        return 127, f"{args[0]}: not found"
    except PermissionError as exc:
        return 126, f"{args[0]}: {exc.strerror}"
    except subprocess.TimeoutExpired:
        return 124, f"{args[0]}: timed out after {timeout:.0f} s"
    return cp.returncode, cp.stdout.decode("utf-8", errors="replace")
