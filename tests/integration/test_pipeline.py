"""End-to-end tests: real FastQC, fastp, MultiQC and pigz; fake prefetch/fasterq-dump.

Every test runs the installed CLI in a subprocess, exactly as a user would.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from typing import Any

import pytest

from ..conftest import ROOT, run_cli, write_metadata

pytestmark = [pytest.mark.tools, pytest.mark.slow]

META_ROWS = [
    "ctrl_1,SRR26053372,paired,siControl",
    "ogt_1,SRR26053373,paired,siOGT",
    "se_1,SRR26053380,single,siOGT",
]


@pytest.fixture
def env(fake_env: dict[str, str]) -> dict[str, str]:
    e = dict(fake_env)
    e["FAKE_LAYOUT_SRR26053380"] = "single"
    return e


def _run(
    tmp_path: Path,
    env: dict[str, str],
    extra: list[str] | None = None,
    rows: list[str] | None = None,
    outdir: Path | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    meta = write_metadata(tmp_path / "meta.csv", rows or META_ROWS)
    out = outdir or tmp_path / "results dir (test) & more"
    r = run_cli(
        [
            "run",
            "--metadata",
            str(meta),
            "--outdir",
            str(out),
            "--threads",
            "2",
            "--no-color",
            *(extra or []),
        ],
        env=env,
    )
    return r, out


def _summary(out: Path) -> dict[str, Any]:
    return json.loads((out / "07_reports" / "summary.json").read_text())


def _sample(summary: dict[str, Any], name: str) -> dict[str, Any]:
    return next(s for s in summary["samples"] if s["sample"] == name)


def _count(path: Path) -> int:
    import gzip

    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def _wrapper(tmp_path: Path, tool: str, fail_when: str, env: dict[str, str]) -> dict[str, str]:
    """Put a wrapper for *tool* first on PATH that fails when *fail_when* is in its arguments."""
    real = shutil.which(tool, path=env["PATH"])
    assert real, tool
    d = tmp_path / f"wrap_{tool}"
    d.mkdir()
    script = d / tool
    script.write_text(
        "#!/bin/sh\n"
        f'case "$*" in *{fail_when}*) echo "injected failure of {tool}" >&2; exit 1;; esac\n'
        f'exec "{real}" "$@"\n'
    )
    script.chmod(0o755)
    e = dict(env)
    e["PATH"] = f"{d}{os.pathsep}{env['PATH']}"
    return e


# --------------------------------------------------------------------------- happy path
def test_full_run_paired_and_single(tmp_path: Path, env: dict[str, str]) -> None:
    env["FAKE_ORPHANS_SRR26053373"] = "37"
    tmpdir = tmp_path / "my tmp"
    r, out = _run(tmp_path, env, ["--tmpdir", str(tmpdir), "--project", "GSE-like test"])
    assert r.returncode == 0, r.stdout + r.stderr
    # Progress explanations shown
    for text in (
        "[STEP 1/6] DOWNLOAD SRA DATA",
        "[STEP 4/6] FASTP PREPROCESSING",
        "Purpose:",
        "RUN SUMMARY",
    ):
        assert text in r.stdout

    s = _summary(out)
    assert s["counts"] == {"attempted": 3, "completed": 3, "failed": 0, "skipped": 0}
    for name, paired in (("ctrl_1", True), ("ogt_1", True), ("se_1", False)):
        smp = _sample(s, name)
        fp_json = json.loads((out / "04_fastp" / name / f"{name}.fastp.json").read_text())
        rows = {row["metric"]: row for row in smp["comparison"]["rows"]}
        finals = sorted((out / "06_final_fastq").glob(f"{name}*.fastq.gz"))
        raws = sorted((out / "02_raw_fastq").glob(f"{name}*.fastq.gz"))
        assert len(finals) == (2 if paired else 1)
        # report values == independent counts == fastp's own numbers
        raw_reads = sum(_count(p) for p in raws)
        clean_reads = sum(_count(p) for p in finals)
        assert (
            rows["Reads"]["before"]
            == raw_reads
            == fp_json["summary"]["before_filtering"]["total_reads"]
        )
        assert (
            rows["Reads"]["after"]
            == clean_reads
            == fp_json["summary"]["after_filtering"]["total_reads"]
        )
        if paired:
            assert _count(finals[0]) == _count(finals[1])  # R1/R2 stay synchronised
            assert rows["Read pairs"]["after"] == clean_reads // 2
        assert rows["Q30 bases"]["after"] == round(
            100 * fp_json["summary"]["after_filtering"]["q30_rate"], 2
        )
        assert smp["comparison"]["consistency_issues"] == []
        # adapter values come from the right FastQC report
        mate = "R1" if paired else "reads"
        stem = f"{name}_R1" if paired else name
        with zipfile.ZipFile(out / "03_raw_qc" / "fastqc" / name / f"{stem}_fastqc.zip") as zf:
            assert f"{stem}_fastqc/fastqc_data.txt" in zf.namelist()
        assert rows[f"Max adapter content ({mate})"]["before"] is not None
        # metadata preserved
        assert smp["metadata"] == {"condition": {"ctrl_1": "siControl"}.get(name, "siOGT")}

    # raw and cleaned MultiQC reports are separate files
    raw_mqc = out / "03_raw_qc" / "multiqc" / "raw_multiqc_report.html"
    clean_mqc = out / "05_clean_qc" / "multiqc" / "clean_multiqc_report.html"
    assert raw_mqc.stat().st_size > 10000 and clean_mqc.stat().st_size > 10000
    raw_inputs = (out / "03_raw_qc" / "multiqc" / "raw_multiqc_report_inputs.txt").read_text()
    assert "04_fastp" not in raw_inputs and "05_clean_qc" not in raw_inputs
    clean_inputs = (out / "05_clean_qc" / "multiqc" / "clean_multiqc_report_inputs.txt").read_text()
    assert "fastp.json" in clean_inputs and "03_raw_qc" not in clean_inputs

    # unpaired reads kept separately, never mixed into R1/R2
    assert _count(out / "02_raw_fastq" / "unpaired" / "ogt_1_unpaired.fastq.gz") == 37
    assert _count(out / "02_raw_fastq" / "ogt_1_R1.fastq.gz") == 1500

    # manifest and checksums
    manifest = list(csv.DictReader((out / "06_final_fastq" / "MANIFEST.csv").open()))
    assert len(manifest) == 5
    for row in manifest:
        data = (out / "06_final_fastq" / row["path"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == row["sha256"]
        assert int(row["reads"]) == _count(out / "06_final_fastq" / row["path"])
    sm = list(csv.DictReader((out / "07_reports" / "sample_manifest.csv").open()))
    assert {r["sample"]: r["condition"] for r in sm}["ctrl_1"] == "siControl"

    # raw data untouched, temp dirs cleaned, no partial files
    assert (out / "01_sra" / "SRR26053372" / "SRR26053372.sra").is_file()
    assert not list(out.rglob("*.partial*")) and not list(out.rglob(".work_*"))
    assert not any(tmpdir.iterdir())
    # logs and run info
    events = [json.loads(line) for line in (out / "logs" / "events.jsonl").read_text().splitlines()]
    assert (
        sum(1 for e in events if e["event"] == "stage_end" and e["status"] == "COMPLETED")
        == 3 * 8 + 2
    )
    versions = json.loads((out / "pipeline_info" / "software_versions.json").read_text())
    assert (
        versions["tools"]["fastp"]["version"]
        == _sample(s, "ctrl_1")["stages"]["fastp"]["tool_version"]
    )
    assert versions["tools"]["fastp"]["version"] == fp_json["summary"]["fastp_version"]
    html = (out / "07_reports" / "summary_report.html").read_text()
    assert "GSE-like test" in html and "All samples completed" in html and "PASS/WARN/FAIL" in html


def test_resume_reuses_and_invalidates(tmp_path: Path, env: dict[str, str]) -> None:
    calls = tmp_path / "prefetch_calls.txt"
    env["FAKE_PREFETCH_COUNT_FILE"] = str(calls)
    r, out = _run(tmp_path, env, rows=META_ROWS[:2])
    assert r.returncode == 0, r.stderr
    final = out / "06_final_fastq" / "ctrl_1_R1.fastq.gz"
    before = final.stat().st_mtime_ns
    assert len(calls.read_text().splitlines()) == 2

    # 1) identical resume: everything reused, no download, outputs untouched
    r, _ = _run(tmp_path, env, ["--resume"], rows=META_ROWS[:2], outdir=out)
    assert r.returncode == 0
    assert r.stdout.count(": REUSED") == 2 * 8 + 2
    assert len(calls.read_text().splitlines()) == 2
    assert final.stat().st_mtime_ns == before

    # 2) without --resume: refused, nothing changed
    r, _ = _run(tmp_path, env, rows=META_ROWS[:2], outdir=out)
    assert r.returncode == 2 and "--resume" in r.stderr

    # 3) changed fastp parameter: fastp and downstream rerun, upstream reused
    r, _ = _run(tmp_path, env, ["--resume", "--min-length", "40"], rows=META_ROWS[:2], outdir=out)
    assert r.returncode == 0
    assert "SRA-to-FASTQ conversion [ctrl_1]: REUSED" in r.stdout
    assert "FastQC (raw reads) [ctrl_1]: REUSED" in r.stdout
    assert "fastp preprocessing [ctrl_1]: COMPLETED" in r.stdout
    assert "FastQC (cleaned reads) [ctrl_1]: COMPLETED" in r.stdout
    assert "MultiQC (raw reads): REUSED" in r.stdout
    fp = _sample(_summary(out), "ctrl_1")["stages"]["fastp"]["metrics"]["command_line"]
    assert fp[fp.index("--length_required") + 1] == "40"
    assert _summary(out)["config"]["fastp"]["length_required"] == 40

    # 4) a supposedly completed output went missing: only affected stages rerun
    (out / "02_raw_fastq" / "ogt_1_R2.fastq.gz").unlink()
    r, _ = _run(tmp_path, env, ["--resume", "--min-length", "40"], rows=META_ROWS[:2], outdir=out)
    assert r.returncode == 0
    assert "SRA-to-FASTQ conversion [ogt_1]: COMPLETED" in r.stdout
    assert "fastp preprocessing [ctrl_1]: REUSED" in r.stdout
    assert len(calls.read_text().splitlines()) == 2  # SRA was still valid: not downloaded again

    # 5) a final file modified by hand is detected and regenerated
    final.write_bytes(b"garbage")
    r, _ = _run(tmp_path, env, ["--resume", "--min-length", "40"], rows=META_ROWS[:2], outdir=out)
    assert r.returncode == 0
    assert "final FASTQ publication [ctrl_1]: COMPLETED" in r.stdout
    assert _count(final) > 0

    # 6) --overwrite recomputes everything
    r, _ = _run(tmp_path, env, ["--overwrite"], rows=META_ROWS[:2], outdir=out)
    assert r.returncode == 0 and ": REUSED" not in r.stdout
    assert list((out / "pipeline_info").glob("state.before-overwrite.*.json"))


def test_refuses_foreign_nonempty_directory(tmp_path: Path, env: dict[str, str]) -> None:
    out = tmp_path / "mine"
    out.mkdir()
    (out / "thesis.docx").write_text("precious")
    for extra in ([], ["--overwrite"], ["--resume"]):
        r, _ = _run(tmp_path, env, extra, outdir=out)
        assert r.returncode == 2 and "not created by this tool" in r.stderr
    assert sorted(p.name for p in out.iterdir()) == ["thesis.docx"]


# --------------------------------------------------------------------------- failures
def test_one_sample_fails_others_complete(tmp_path: Path, env: dict[str, str]) -> None:
    env["FAKE_FASTERQ_FAIL"] = "SRR26053373"
    r, out = _run(tmp_path, env)
    assert r.returncode == 4  # partial failure
    s = _summary(out)
    assert s["counts"] == {"attempted": 3, "completed": 2, "failed": 1, "skipped": 0}
    ogt = _sample(s, "ogt_1")
    assert ogt["status"] == "FAILED" and ogt["failure"]["stage"] == "convert"
    assert ogt["stages"]["fastp"]["status"] == "SKIPPED"
    assert "corrupt" in ogt["failure"]["hint"] or "SRA" in ogt["failure"]["hint"]
    assert not list((out / "06_final_fastq").glob("ogt_1*"))
    assert "ogt_1: FAILED" in r.stdout and "ERROR [ogt_1]" in r.stderr
    status = list(csv.DictReader((out / "07_reports" / "processing_status.csv").open()))
    assert {row["sample"]: row["overall_status"] for row in status} == {
        "ctrl_1": "COMPLETED",
        "ogt_1": "FAILED",
        "se_1": "COMPLETED",
    }
    # fix the cause and resume: only the failed sample is processed
    env.pop("FAKE_FASTERQ_FAIL")
    r, _ = _run(tmp_path, env, ["--resume"], outdir=out)
    assert r.returncode == 0
    assert "fastp preprocessing [ctrl_1]: REUSED" in r.stdout
    assert "fastp preprocessing [ogt_1]: COMPLETED" in r.stdout


def test_on_error_stop(tmp_path: Path, env: dict[str, str]) -> None:
    env["FAKE_FASTERQ_FAIL"] = "SRR26053372"
    r, out = _run(tmp_path, env, ["--on-error", "stop"])
    assert r.returncode == 5
    s = _summary(out)
    assert _sample(s, "ctrl_1")["status"] == "FAILED"
    assert _sample(s, "ogt_1")["stages"]["download"]["status"] == "SKIPPED"
    assert "on_error: stop" in _sample(s, "ogt_1")["stages"]["download"]["message"]


def test_network_failure_is_retried_then_reported(tmp_path: Path, env: dict[str, str]) -> None:
    calls = tmp_path / "calls.txt"
    env.update(FAKE_PREFETCH_NETWORK_FAIL="SRR26053372", FAKE_PREFETCH_COUNT_FILE=str(calls))
    cfg = tmp_path / "c.yaml"
    cfg.write_text("download:\n  retries: 2\n  retry_wait_seconds: 0\n")
    r, out = _run(tmp_path, env, ["--config", str(cfg)], rows=META_ROWS[:1])
    assert r.returncode == 5
    assert calls.read_text().splitlines() == ["SRR26053372"] * 3
    msg = _sample(_summary(out), "ctrl_1")["failure"]["message"]
    assert "network" in msg and "3 attempt" in msg


def test_invalid_accession_not_retried(tmp_path: Path, env: dict[str, str]) -> None:
    calls = tmp_path / "calls.txt"
    env.update(FAKE_PREFETCH_PERMANENT_FAIL="SRR26053372", FAKE_PREFETCH_COUNT_FILE=str(calls))
    r, out = _run(tmp_path, env, rows=META_ROWS[:1])
    assert r.returncode == 5
    assert len(calls.read_text().splitlines()) == 1
    assert "invalid" in _sample(_summary(out), "ctrl_1")["failure"]["message"]


def test_download_reported_success_but_no_file(tmp_path: Path, env: dict[str, str]) -> None:
    env["FAKE_PREFETCH_NO_OUTPUT"] = "SRR26053372"
    r, out = _run(tmp_path, env, rows=META_ROWS[:1])
    assert r.returncode == 5
    assert "no SRR26053372.sra" in _sample(_summary(out), "ctrl_1")["failure"]["message"]


@pytest.mark.parametrize(
    ("variable", "row", "stage", "fragment"),
    [
        # the last record loses 30 characters: detected at that record
        ("FAKE_FASTERQ_TRUNCATE", META_ROWS[0], "validate_raw", "record 1500"),
        ("FAKE_FASTERQ_MISMATCH", META_ROWS[0], "validate_raw", "different numbers of reads"),
        ("FAKE_FASTERQ_NOSPACE", META_ROWS[0], "convert", "exited with code 3"),
    ],
)
def test_conversion_problems(
    tmp_path: Path, env: dict[str, str], variable: str, row: str, stage: str, fragment: str
) -> None:
    env[variable] = "SRR26053372"
    r, out = _run(tmp_path, env, rows=[row])
    assert r.returncode == 5
    f = _sample(_summary(out), "ctrl_1")["failure"]
    assert f["stage"] == stage and fragment in f["message"]
    if variable == "FAKE_FASTERQ_NOSPACE":
        assert "disk space" in f["hint"]
    assert not list((out / "06_final_fastq").glob("*.fastq.gz"))


def test_layout_mismatch(tmp_path: Path, env: dict[str, str]) -> None:
    r, out = _run(tmp_path, env, rows=["x,SRR26053380,paired,c"])  # data are single-end
    assert r.returncode == 5
    f = _sample(_summary(out), "x")["failure"]
    assert "declared paired-end" in f["message"] and "single" in f["hint"]
    r, out = _run(tmp_path, env, rows=["y,SRR26053372,single,c"], outdir=tmp_path / "o2")
    assert "declared single-end" in _sample(_summary(out), "y")["failure"]["message"]


@pytest.mark.parametrize(("tool", "stage"), [("fastp", "fastp"), ("fastqc", "fastqc_raw")])
def test_tool_failure_for_one_sample(
    tmp_path: Path, env: dict[str, str], tool: str, stage: str
) -> None:
    e = _wrapper(tmp_path, tool, "ogt_1_R1", env)
    r, out = _run(tmp_path, e, rows=META_ROWS[:2])
    assert r.returncode == 4
    s = _sample(_summary(out), "ogt_1")
    assert s["failure"]["stage"] == stage and "exited with code 1" in s["failure"]["message"]
    assert _sample(_summary(out), "ctrl_1")["status"] == "COMPLETED"
    log = out / "logs" / "commands" / "ogt_1" / f"{stage}.log"
    assert "injected failure" in log.read_text()


def test_multiqc_failure_after_successful_cleaning(tmp_path: Path, env: dict[str, str]) -> None:
    e = _wrapper(tmp_path, "multiqc", "multiqc", env)
    r, out = _run(tmp_path, e, rows=META_ROWS[:1])
    assert r.returncode == 4
    s = _summary(out)
    assert _sample(s, "ctrl_1")["status"] == "COMPLETED"
    assert s["aggregates"]["multiqc_raw"]["status"] == "FAILED"
    assert (out / "06_final_fastq" / "ctrl_1_R1.fastq.gz").is_file()
    assert (out / "07_reports" / "summary_report.html").is_file()


def test_interrupt_and_resume(tmp_path: Path, env: dict[str, str]) -> None:
    env["FAKE_FASTERQ_SLEEP"] = "60"
    meta = write_metadata(tmp_path / "meta.csv", META_ROWS[:1])
    out = tmp_path / "out"
    proc = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "transcriptomics_preprocessor",
            "run",
            "--metadata",
            str(meta),
            "--outdir",
            str(out),
            "--threads",
            "2",
            "--no-color",
        ],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    log = out / "logs" / "commands" / "ctrl_1" / "convert.log"
    deadline = time.time() + 60
    while not (log.exists() and "fasterq-dump" in log.read_text()):
        assert time.time() < deadline, "conversion never started"
        time.sleep(0.2)
    t0 = time.time()
    proc.send_signal(signal.SIGINT)
    stdout, stderr = proc.communicate(timeout=60)
    assert proc.returncode == 130, stdout + stderr
    assert time.time() - t0 < 30  # did not wait for the 60 s tool
    state = json.loads((out / "pipeline_info" / "state.json").read_text())
    stages = state["samples"]["ctrl_1"]["stages"]
    assert stages["download"]["status"] == "COMPLETED"
    assert stages["convert"]["status"] == "FAILED"
    assert "RUNNING" not in json.dumps(state)
    assert "--resume" in stdout
    assert (out / "07_reports" / "summary_report.html").exists()  # partial report written
    assert not list((out / "02_raw_fastq").glob("*.fastq.gz"))  # no partial output published

    env.pop("FAKE_FASTERQ_SLEEP")
    r = run_cli(
        ["run", "--metadata", str(meta), "--outdir", str(out), "--threads", "2", "--resume"],
        env=env,
    )
    assert r.returncode == 0
    assert "SRA download [ctrl_1]: REUSED" in r.stdout


def test_corrupt_state_file(tmp_path: Path, env: dict[str, str]) -> None:
    r, out = _run(tmp_path, env, rows=META_ROWS[:1])
    assert r.returncode == 0
    (out / "pipeline_info" / "state.json").write_text("{broken")
    r, _ = _run(tmp_path, env, ["--resume"], rows=META_ROWS[:1], outdir=out)
    assert r.returncode == 5 and "corrupt" in r.stderr and "--overwrite" in r.stderr
    r, _ = _run(tmp_path, env, ["--overwrite"], rows=META_ROWS[:1], outdir=out)
    assert r.returncode == 0


# --------------------------------------------------------------------------- options
def test_dry_run_creates_nothing(tmp_path: Path, env: dict[str, str]) -> None:
    out = tmp_path / "planned"
    r, _ = _run(tmp_path, env, ["--dry-run", "--min-length", "30"], outdir=out)
    assert r.returncode == 0, r.stdout + r.stderr
    assert not out.exists()
    for text in (
        "DRY RUN",
        "prefetch SRR26053372",
        "fasterq-dump",
        "--length_required 30",
        "multiqc",
    ):
        assert text in r.stdout


def test_dry_run_on_existing_run_shows_reuse_plan(tmp_path: Path, env: dict[str, str]) -> None:
    r, out = _run(tmp_path, env, rows=META_ROWS[:1])
    state_before = (out / "pipeline_info" / "state.json").read_bytes()
    r, _ = _run(
        tmp_path,
        env,
        ["--dry-run", "--resume", "--min-length", "50"],
        rows=META_ROWS[:1],
        outdir=out,
    )
    assert r.returncode == 0
    assert "SRA download: reuse" in r.stdout and "fastp preprocessing: run" in r.stdout
    assert (out / "pipeline_info" / "state.json").read_bytes() == state_before
    r, _ = _run(tmp_path, env, ["--dry-run"], rows=META_ROWS[:1], outdir=out)
    assert r.returncode == 2 and "needs --resume or --overwrite" in r.stdout


def test_remove_sra_option(tmp_path: Path, env: dict[str, str]) -> None:
    calls = tmp_path / "calls.txt"
    env["FAKE_PREFETCH_COUNT_FILE"] = str(calls)
    r, out = _run(tmp_path, env, ["--remove-sra"], rows=META_ROWS[:1])
    assert r.returncode == 0
    assert not (out / "01_sra" / "SRR26053372" / "SRR26053372.sra").exists()
    assert (out / "02_raw_fastq" / "ctrl_1_R1.fastq.gz").exists()  # raw FASTQ kept
    r, _ = _run(tmp_path, env, ["--resume", "--remove-sra"], rows=META_ROWS[:1], outdir=out)
    assert r.returncode == 0 and len(calls.read_text().splitlines()) == 1  # no re-download needed
    (
        out / "02_raw_fastq" / "ctrl_1_R1.fastq.gz"
    ).unlink()  # conversion must be redone -> SRA needed
    r, _ = _run(tmp_path, env, ["--resume", "--remove-sra"], rows=META_ROWS[:1], outdir=out)
    assert r.returncode == 0 and len(calls.read_text().splitlines()) == 2


def test_parallel_jobs_and_quiet(tmp_path: Path, env: dict[str, str]) -> None:
    r, out = _run(tmp_path, env, ["--threads", "4", "--jobs", "3", "--quiet"])
    assert r.returncode == 0, r.stderr
    assert "[STEP" not in r.stdout and "RUN SUMMARY" in r.stdout
    assert _summary(out)["counts"]["completed"] == 3
    assert (out / "logs" / "pipeline.log").read_text().count("[STEP") >= 24  # still logged


def test_verbose_shows_commands(tmp_path: Path, env: dict[str, str]) -> None:
    r, _ = _run(tmp_path, env, ["--verbose"], rows=META_ROWS[2:])
    assert r.returncode == 0
    assert r.stdout.count("Purpose:") >= 6  # explanations repeated in verbose mode


def test_quick_validation_reports_counts_as_unavailable(
    tmp_path: Path, env: dict[str, str]
) -> None:
    cfg = tmp_path / "q.yaml"
    cfg.write_text("validation:\n  mode: quick\n  quick_records: 1000\n")  # file has 1500 reads
    r, out = _run(tmp_path, env, ["--config", str(cfg)], rows=META_ROWS[2:])
    assert r.returncode == 0
    rows = {x["metric"]: x for x in _sample(_summary(out), "se_1")["comparison"]["rows"]}
    assert rows["Reads"]["before"] is None and rows["Reads"]["change"] == "not available"
    assert rows["Bases"]["before"] == 150000  # fastp's own count still reported


def test_explicit_adapter_fixes_single_end_detection(tmp_path: Path, env: dict[str, str]) -> None:
    r, out = _run(
        tmp_path, env, ["--adapter-r1", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"], rows=META_ROWS[2:]
    )
    assert r.returncode == 0
    rows = {x["metric"]: x for x in _sample(_summary(out), "se_1")["comparison"]["rows"]}
    assert (
        rows["Max adapter content (reads)"]["after"] < rows["Max adapter content (reads)"]["before"]
    )
    assert rows["Max adapter content (reads)"]["after"] < 5


def test_accessions_input(tmp_path: Path, env: dict[str, str]) -> None:
    out = tmp_path / "acc"
    r = run_cli(
        [
            "run",
            "--accessions",
            "SRR26053372",
            "--layout",
            "paired",
            "--outdir",
            str(out),
            "--threads",
            "2",
        ],
        env=env,
    )
    assert r.returncode == 0
    assert (out / "06_final_fastq" / "SRR26053372_R1.fastq.gz").is_file()


@pytest.mark.skipif(os.geteuid() == 0, reason="permission checks are not enforced for root")
def test_unwritable_output_location(tmp_path: Path, env: dict[str, str]) -> None:
    ro = tmp_path / "ro"
    ro.mkdir()
    ro.chmod(0o500)
    try:
        r, _ = _run(tmp_path, env, outdir=ro / "out")
        assert r.returncode == 2 and "cannot create" in r.stderr
    finally:
        ro.chmod(0o700)


def test_package_source_has_no_shell_true() -> None:
    """No subprocess call in the package passes shell=True (checked on the syntax tree)."""
    import ast

    src = ROOT / "src" / "transcriptomics_preprocessor"
    for p in src.glob("*.py"):
        for node in ast.walk(ast.parse(p.read_text())):
            if isinstance(node, ast.Call):
                for kw in node.keywords:
                    assert not (
                        kw.arg == "shell" and isinstance(kw.value, ast.Constant) and kw.value.value
                    ), p
    # and the scanner itself works:
    tree = ast.parse("subprocess.run('x', shell=True)")
    call = next(n for n in ast.walk(tree) if isinstance(n, ast.Call))
    assert any(k.arg == "shell" for k in call.keywords)


def test_interrupt_with_parallel_jobs(tmp_path: Path, env: dict[str, str]) -> None:
    """SIGTERM while two samples convert concurrently: both stop, journal consistent, resume works."""
    env["FAKE_FASTERQ_SLEEP"] = "60"
    meta = write_metadata(tmp_path / "meta.csv", META_ROWS[:2])
    out = tmp_path / "out"
    proc = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "transcriptomics_preprocessor",
            "run",
            "--metadata",
            str(meta),
            "--outdir",
            str(out),
            "--threads",
            "2",
            "--jobs",
            "2",
            "--no-color",
        ],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    logs = [out / "logs" / "commands" / s / "convert.log" for s in ("ctrl_1", "ogt_1")]
    deadline = time.time() + 90
    while not all(p.exists() and "fasterq-dump" in p.read_text() for p in logs):
        assert time.time() < deadline, "conversions never started"
        time.sleep(0.2)
    proc.send_signal(signal.SIGTERM)
    stdout, stderr = proc.communicate(timeout=60)
    assert proc.returncode == 130, stdout + stderr
    state = json.loads((out / "pipeline_info" / "state.json").read_text())
    assert "RUNNING" not in json.dumps(state)
    for s in ("ctrl_1", "ogt_1"):
        assert state["samples"][s]["stages"]["convert"]["status"] == "FAILED"
    env.pop("FAKE_FASTERQ_SLEEP")
    r = run_cli(
        [
            "run",
            "--metadata",
            str(meta),
            "--outdir",
            str(out),
            "--threads",
            "2",
            "--jobs",
            "2",
            "--resume",
        ],
        env=env,
    )
    assert r.returncode == 0
    assert r.stdout.count("SRA download [") == 2 and "REUSED" in r.stdout
