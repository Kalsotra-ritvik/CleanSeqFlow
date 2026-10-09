from __future__ import annotations

import json
import os
import sys
import threading
import time
from pathlib import Path

import pytest

from transcriptomics_preprocessor.errors import (
    InterruptedRun,
    OutputDirectoryError,
    StateError,
    ToolExecutionError,
)
from transcriptomics_preprocessor.runner import STOP_EVENT, probe_output, run_command
from transcriptomics_preprocessor.state import (
    RunState,
    Status,
    can_reuse,
    exclusive_lock,
    outputs_intact,
    outputs_record,
)
from transcriptomics_preprocessor.utils import (
    atomic_write_text,
    human_bytes,
    parse_size,
    safe_rmtree,
    stable_hash,
)

PY = sys.executable


# ----------------------------------------------------------------- state
def test_state_roundtrip_and_backup(tmp_path: Path) -> None:
    p = tmp_path / "state.json"
    st = RunState.new(p)
    st.update("s1", "download", status=Status.COMPLETED, fingerprint="abc")
    st.update("s1", "convert", status=Status.RUNNING)
    again = RunState.load(p)
    assert again.status("s1", "download") == Status.COMPLETED
    assert again.status("s1", "fastp") == Status.NOT_STARTED
    assert (tmp_path / "state.json.bak").exists()
    again.mark_interrupted()
    assert RunState.load(p).status("s1", "convert") == Status.FAILED


@pytest.mark.parametrize(
    "content", ["{not json", '{"schema_version": 99}', "[]", '{"schema_version": 1}']
)
def test_corrupt_state_rejected(tmp_path: Path, content: str) -> None:
    p = tmp_path / "state.json"
    p.write_text(content)
    with pytest.raises(StateError):
        RunState.load(p)


def test_reuse_requires_fingerprint_and_intact_outputs(tmp_path: Path) -> None:
    f = tmp_path / "out.txt"
    f.write_text("hello")
    rec = {
        "status": "COMPLETED",
        "fingerprint": "fp1",
        "outputs": outputs_record({"x": f}, tmp_path),
    }
    assert can_reuse(rec, "fp1", tmp_path)[0]
    assert not can_reuse(rec, "fp2", tmp_path)[0]  # parameters/inputs changed
    assert not can_reuse({**rec, "status": "FAILED"}, "fp1", tmp_path)[0]
    f.write_text("hello!")  # modified (size changes)
    ok, why = can_reuse(rec, "fp1", tmp_path)
    assert not ok and "changed" in why
    f.unlink()
    ok, why = outputs_intact(rec, tmp_path)
    assert not ok and "missing" in why
    assert not outputs_intact({"status": "COMPLETED"}, tmp_path)[0]


def test_same_size_modification_detected_by_mtime(tmp_path: Path) -> None:
    f = tmp_path / "out.txt"
    f.write_text("aaaa")
    rec = {"status": "COMPLETED", "fingerprint": "x", "outputs": outputs_record({"x": f}, tmp_path)}
    time.sleep(0.01)
    f.write_text("bbbb")
    os.utime(f, ns=(time.time_ns(), time.time_ns()))
    assert not can_reuse(rec, "x", tmp_path)[0]


def test_lock_prevents_concurrent_runs(tmp_path: Path) -> None:
    lock = tmp_path / "run.lock"
    with (
        exclusive_lock(lock),
        pytest.raises(OutputDirectoryError, match="another run"),
        exclusive_lock(lock),
    ):
        pass
    with exclusive_lock(lock):  # released again
        pass


def test_stable_hash_order_independent() -> None:
    assert stable_hash({"a": 1, "b": [1, 2]}) == stable_hash({"b": [1, 2], "a": 1})
    assert stable_hash({"a": 1}) != stable_hash({"a": 2})


# ----------------------------------------------------------------- utils
def test_sizes() -> None:
    assert (
        parse_size("1K") == 1024
        and parse_size("20G") == 20 * 1024**3
        and parse_size("1.5M") == int(1.5 * 1024**2)
    )
    assert parse_size("100GiB") == 100 * 1024**3
    for bad in ("", "abc", "-1G", "nanG"):
        with pytest.raises(ValueError):
            parse_size(bad)
    assert (
        human_bytes(1536) == "1.5 KiB"
        and human_bytes(None) == "unknown"
        and human_bytes(10) == "10 B"
    )


def test_safe_rmtree_guards(tmp_path: Path) -> None:
    root = tmp_path / "root"
    (root / "sub").mkdir(parents=True)
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "keep.txt").write_text("x")
    link = root / "link"
    link.symlink_to(outside)
    safe_rmtree(link, root)  # removes the link only
    assert (outside / "keep.txt").exists()
    with pytest.raises(PermissionError):
        safe_rmtree(outside, root)
    with pytest.raises(PermissionError):
        safe_rmtree(root, root)
    safe_rmtree(root / "sub", root)
    assert not (root / "sub").exists()


def test_atomic_write(tmp_path: Path) -> None:
    p = tmp_path / "a" / "b.txt"
    atomic_write_text(p, "x")
    assert p.read_text() == "x"
    assert not [q for q in p.parent.iterdir() if q.name.endswith(".tmp")]


# ----------------------------------------------------------------- runner
def test_run_command_success_logs(tmp_path: Path) -> None:
    log = tmp_path / "l.log"
    r = run_command([PY, "-c", "print('hi')"], log_path=log, description="say hi")
    assert r.returncode == 0
    text = log.read_text()
    assert "say hi" in text and "hi" in text and "exit code: 0" in text


def test_run_command_failure_raises_with_tail(tmp_path: Path) -> None:
    log = tmp_path / "l.log"
    with pytest.raises(ToolExecutionError) as e:
        run_command([PY, "-c", "import sys; print('boom'); sys.exit(3)"], log_path=log)
    assert (
        e.value.returncode == 3
        and "boom" in e.value.output_tail
        and "exited with code 3" in e.value.message
    )
    r = run_command([PY, "-c", "import sys; sys.exit(2)"], log_path=log, check=False)
    assert r.returncode == 2


def test_run_command_timeout(tmp_path: Path) -> None:
    t0 = time.monotonic()
    with pytest.raises(ToolExecutionError, match="timed out"):
        run_command([PY, "-c", "import time; time.sleep(30)"], log_path=tmp_path / "l", timeout=0.5)
    assert time.monotonic() - t0 < 15


def test_run_command_missing_program(tmp_path: Path) -> None:
    with pytest.raises(ToolExecutionError, match="program not found"):
        run_command(["definitely-not-a-program-xyz"], log_path=tmp_path / "l")
    assert probe_output(["definitely-not-a-program-xyz"])[0] == 127


def test_run_command_interrupt_kills_process_group(tmp_path: Path) -> None:
    pidfile = tmp_path / "child.pid"
    script = (
        "import subprocess,sys,time;"
        f"c=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)']);"
        f"open({str(pidfile)!r},'w').write(str(c.pid));time.sleep(60)"
    )
    threading.Timer(1.0, STOP_EVENT.set).start()
    with pytest.raises(InterruptedRun):
        run_command([PY, "-c", script], log_path=tmp_path / "l")
    child = int(pidfile.read_text())
    time.sleep(0.5)
    # The grandchild must not be running any more. It may linger as a zombie ("Z") until
    # the system's init process reaps it, which is outside our control.
    stat = Path(f"/proc/{child}/stat")
    state = stat.read_text().split()[2] if stat.exists() else "gone"
    assert state in ("Z", "X", "gone"), f"grandchild still running (state {state})"


def test_run_command_no_shell_injection(tmp_path: Path) -> None:
    marker = tmp_path / "pwned"
    run_command(
        [PY, "-c", "import sys; print(sys.argv[1])", f"; touch {marker}"], log_path=tmp_path / "l"
    )
    assert not marker.exists()


def test_stdout_capture(tmp_path: Path) -> None:
    out = tmp_path / "o.bin"
    run_command(
        [PY, "-c", "import sys; sys.stdout.write('data')"], log_path=tmp_path / "l", stdout_path=out
    )
    assert out.read_text() == "data"


def test_state_json_is_valid_json(tmp_path: Path) -> None:
    st = RunState.new(tmp_path / "s.json")
    st.update(None, "multiqc_raw", status=Status.SKIPPED)
    assert (
        json.loads((tmp_path / "s.json").read_text())["aggregate"]["multiqc_raw"]["status"]
        == "SKIPPED"
    )
