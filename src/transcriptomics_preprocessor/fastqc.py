"""FastQC stage (used on raw reads and again on cleaned reads)."""

from __future__ import annotations

from pathlib import Path

from .context import StageContext, StageResult
from .errors import IntegrityError
from .layout import OutputLayout
from .metrics import MetricsError, read_fastqc_zip
from .runner import run_command
from .utils import remove_if_exists, safe_rmtree
from .validation import Sample


def expected_outputs(fastq: Path, outdir: Path) -> tuple[Path, Path]:
    stem = OutputLayout.fastqc_stem(fastq)
    return outdir / f"{stem}_fastqc.html", outdir / f"{stem}_fastqc.zip"


def build_command(
    fastqs: list[Path], outdir: Path, tmp: Path, threads: int, memory_mb: int
) -> list[str]:
    return [
        "fastqc",
        "--outdir",
        str(outdir),
        "--threads",
        str(max(1, min(threads, len(fastqs)))),
        "--memory",
        str(memory_mb),
        "--dir",
        str(tmp),
        "--noextract",
        "--quiet",
        *[str(f) for f in fastqs],
    ]


def run(
    sample: Sample, fastqs: list[Path], outdir: Path, stage: str, ctx: StageContext
) -> StageResult:
    """Run FastQC on a sample's files and parse the resulting reports."""
    result = StageResult()
    outdir.mkdir(parents=True, exist_ok=True)
    labels = ["R1", "R2"] if sample.paired else ["reads"]
    expected = [expected_outputs(f, outdir) for f in fastqs]
    for html, zf in expected:  # never let a stale report masquerade as a new one
        remove_if_exists(html)
        remove_if_exists(zf)
    tmp = ctx.layout.tmp_for(ctx.tmp_root, sample, stage)
    if tmp.exists():
        safe_rmtree(tmp, ctx.tmp_root)
    tmp.mkdir(parents=True)
    run_command(
        build_command(fastqs, outdir, tmp, ctx.threads, int(ctx.config["fastqc.memory_mb"])),
        log_path=ctx.layout.command_log(sample, stage),
        timeout=float(ctx.config["fastqc.timeout_hours"]) * 3600,
        description=f"FastQC {sample.name} ({stage})",
    )
    safe_rmtree(tmp, ctx.tmp_root)
    for label, fq, (html, zf) in zip(labels, fastqs, expected, strict=True):
        for p in (html, zf):
            if not p.is_file() or p.stat().st_size == 0:
                raise IntegrityError(
                    f"FastQC finished but {p.name} is missing or empty",
                    hint=f"See {ctx.layout.command_log(sample, stage)}",
                )
        try:
            m = read_fastqc_zip(zf)
        except MetricsError as exc:
            raise IntegrityError(f"FastQC report for {fq.name} cannot be parsed: {exc}") from exc
        result.outputs[f"{label}_html"] = html
        result.outputs[f"{label}_zip"] = zf
        result.metrics[label] = m.to_dict()
        flags = m.module_status
        n_fail = sum(1 for v in flags.values() if v == "FAIL")
        n_warn = sum(1 for v in flags.values() if v == "WARN")
        result.notes.append(
            f"{fq.name}: {m.total_sequences if m.total_sequences is not None else 'n/a'} sequences; "
            f"FastQC flags: {len(flags) - n_fail - n_warn} PASS, {n_warn} WARN, {n_fail} FAIL "
            "(flags are screening signals, not a verdict on usability)."
        )
    return result
