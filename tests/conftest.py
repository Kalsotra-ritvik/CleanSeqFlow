"""Shared fixtures.

Integration tests use the real FastQC, fastp, MultiQC and pigz when they are
on PATH (tests marked ``tools`` are skipped otherwise) and always use the
test doubles in ``tests/fixtures/fake_bin`` for ``prefetch`` and
``fasterq-dump``, so no network access to NCBI is needed.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
FAKE_BIN = ROOT / "tests" / "fixtures" / "fake_bin"
REAL_TOOLS = ("fastqc", "fastp", "multiqc", "pigz")


def tools_available() -> bool:
    return all(shutil.which(t) for t in REAL_TOOLS)


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if tools_available():
        return
    skip = pytest.mark.skip(
        reason=f"needs {', '.join(REAL_TOOLS)} on PATH (run inside the installed environment)"
    )
    for item in items:
        if "tools" in item.keywords:
            item.add_marker(skip)


@pytest.fixture(autouse=True)
def _reset_stop_event() -> Iterator[None]:
    from transcriptomics_preprocessor.runner import STOP_EVENT

    STOP_EVENT.clear()
    yield
    STOP_EVENT.clear()


@pytest.fixture
def fake_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> dict[str, str]:
    """Environment with the fake SRA tools first on PATH; returns a copy for subprocesses."""
    for key in list(os.environ):
        if key.startswith("FAKE_"):
            monkeypatch.delenv(key, raising=False)
    monkeypatch.setenv("PATH", f"{FAKE_BIN}{os.pathsep}{os.environ['PATH']}")
    # Make the package importable by the fake tools even when not installed.
    src = str(ROOT / "src")
    monkeypatch.setenv(
        "PYTHONPATH",
        src + (os.pathsep + os.environ["PYTHONPATH"] if os.environ.get("PYTHONPATH") else ""),
    )
    monkeypatch.setenv("JAVA_TOOL_OPTIONS", "")
    return dict(os.environ)


def run_cli(
    args: list[str],
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
    timeout: float = 600,
) -> subprocess.CompletedProcess[str]:
    """Run the CLI in a subprocess (as a user would) and capture output."""
    return subprocess.run(
        [sys.executable, "-m", "transcriptomics_preprocessor", *args],
        env=env,
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )


def write_metadata(
    path: Path, rows: list[str], header: str = "sample,accession,layout,condition"
) -> Path:
    path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path
