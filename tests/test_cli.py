"""CLI behaviour that does not require external tools."""

from __future__ import annotations

from pathlib import Path

import pytest

from transcriptomics_preprocessor import __version__
from transcriptomics_preprocessor.cli import main

from .conftest import run_cli, write_metadata


def test_help_lists_commands_and_examples(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["--help"]) == 0
    out = capsys.readouterr().out
    for word in ("run", "check-input", "doctor", "explain", "--dry-run", "--resume", "exit codes"):
        assert word in out


def test_run_help_documents_options(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["run", "--help"]) == 0
    out = capsys.readouterr().out
    for opt in (
        "--metadata",
        "--accessions",
        "--layout",
        "--outdir",
        "--threads",
        "--jobs",
        "--tmpdir",
        "--overwrite",
        "--remove-sra",
        "--min-length",
        "--quiet",
        "--verbose",
        "--on-error",
    ):
        assert opt in out


def test_version() -> None:
    r = run_cli(["--version"])
    assert r.returncode == 0 and __version__ in r.stdout


def test_no_arguments_prints_help_and_fails() -> None:
    r = run_cli([])
    assert r.returncode == 2 and "usage" in r.stdout


@pytest.mark.parametrize(
    "args",
    [
        ["frobnicate"],
        ["run", "--outdir", "x"],  # no input
        [
            "run",
            "--accessions",
            "SRR1234567",
            "--metadata",
            "m.csv",
            "--outdir",
            "x",
        ],  # both inputs
        ["run", "--accessions", "SRR1234567", "--layout", "triple", "--outdir", "x"],
        ["run", "--accessions", "SRR1234567", "--layout", "paired"],  # no outdir
        ["run", "--accessions", "SRR1234567", "--layout", "paired", "-o", "x", "--threads", "0"],
        ["run", "--accessions", "SRR1234567", "--layout", "paired", "-o", "x", "--threads", "abc"],
        [
            "run",
            "--accessions",
            "SRR1234567",
            "--layout",
            "paired",
            "-o",
            "x",
            "--resume",
            "--overwrite",
        ],
        [
            "run",
            "--accessions",
            "SRR1234567",
            "--layout",
            "paired",
            "-o",
            "x",
            "--quiet",
            "--verbose",
        ],
        ["run", "--accessions", "SRR1234567", "--layout", "paired", "-o", "x", "--retries", "-1"],
    ],
)
def test_usage_errors_exit_2(args: list[str]) -> None:
    r = run_cli(args)
    assert r.returncode == 2, r.stderr


def test_semantic_validation_errors_exit_2(tmp_path: Path) -> None:
    r = run_cli(
        ["run", "--accessions", "GSE243105", "--layout", "paired", "-o", str(tmp_path / "o")]
    )
    assert r.returncode == 2 and "GEO Series" in r.stderr and "Run Selector" in r.stderr
    assert not (tmp_path / "o").exists()  # nothing created for invalid input
    r = run_cli(
        [
            "run",
            "--accessions",
            "SRR1234567",
            "--layout",
            "paired",
            "-o",
            str(tmp_path / "o"),
            "--threads",
            "2",
            "--jobs",
            "3",
        ]
    )
    assert r.returncode == 2 and "jobs (3) cannot exceed threads (2)" in r.stderr
    r = run_cli(
        [
            "run",
            "--accessions",
            "SRR1234567",
            "--layout",
            "paired",
            "-o",
            str(tmp_path / "o"),
            "--max-download-size",
            "huge",
        ]
    )
    assert r.returncode == 2 and "download.max_size" in r.stderr
    r = run_cli(
        [
            "run",
            "--accessions",
            "SRR1234567",
            "--layout",
            "paired",
            "-o",
            str(tmp_path / "o"),
            "--adapter-r1",
            "NOTDNA!",
        ]
    )
    assert r.returncode == 2 and "adapter" in r.stderr


def test_check_input_valid(tmp_path: Path) -> None:
    m = write_metadata(
        tmp_path / "m.csv", ["a,SRR1234567,paired,ctrl", "b,SRR7654321,single,treat"]
    )
    r = run_cli(["check-input", "--metadata", str(m)])
    assert r.returncode == 0
    assert "2 valid sample(s)" in r.stdout and "paired-end: 1, single-end: 1" in r.stdout


def test_check_input_invalid_lists_every_problem(tmp_path: Path) -> None:
    m = write_metadata(tmp_path / "m.csv", ["a,SRR1234567,paired,x", "a,SRP000001,paired,y"])
    r = run_cli(["check-input", "--metadata", str(m)])
    assert r.returncode == 2
    assert "appears 2 times" in r.stderr and "SRA Study" in r.stderr


def test_check_input_with_config(tmp_path: Path) -> None:
    c = tmp_path / "c.yaml"
    c.write_text("fastp:\n  length_required: 0\n")
    r = run_cli(
        ["check-input", "--accessions", "SRR1234567", "--layout", "single", "--config", str(c)]
    )
    assert r.returncode == 2 and "fastp.length_required" in r.stderr


@pytest.mark.parametrize("topic", ["steps", "qc", "limitations", "all"])
def test_explain(topic: str) -> None:
    r = run_cli(["explain", topic])
    assert r.returncode == 0 and r.stdout.strip()
    if topic in ("qc", "all"):
        assert "duplication" in r.stdout.lower() and "PASS/WARN/FAIL" in r.stdout


def test_doctor_reports_missing_tools(tmp_path: Path) -> None:
    import os

    env = dict(os.environ, PATH=str(tmp_path))  # no tools reachable
    r = run_cli(["doctor"], env=env)
    assert r.returncode == 3
    assert "FAIL" in r.stdout and "not found on PATH" in r.stdout


def test_doctor_json(tmp_path: Path) -> None:
    import json
    import os

    r = run_cli(
        ["doctor", "--json", "--outdir", str(tmp_path / "o")],
        env=dict(os.environ, PATH=str(tmp_path)),
    )
    data = json.loads(r.stdout)
    assert data["ok"] is False and "checks" in data and "fastp" in data["tools"]


def test_run_missing_dependency_exit_3(tmp_path: Path) -> None:
    import os

    r = run_cli(
        ["run", "--accessions", "SRR1234567", "--layout", "paired", "-o", str(tmp_path / "o")],
        env=dict(os.environ, PATH=str(tmp_path)),
    )
    assert r.returncode == 3
    assert "missing or broken" in r.stderr and "install.sh" in r.stderr
