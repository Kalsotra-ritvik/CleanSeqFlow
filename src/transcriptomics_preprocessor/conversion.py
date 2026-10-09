"""Stage 2: convert the downloaded SRA archive to FASTQ with ``fasterq-dump``.

``fasterq-dump --split-3`` writes, for paired-end runs, ``<acc>_1.fastq`` and
``<acc>_2.fastq`` (mates in the same order) plus, if some spots have only one
biological read, an unsuffixed ``<acc>.fastq`` with those mate-less reads. For
single-end runs it writes ``<acc>.fastq``. Technical reads (barcodes,
adapters recorded as separate reads) are skipped.

The files are compressed with ``pigz`` and renamed to
``<sample>_R1.fastq.gz`` / ``<sample>_R2.fastq.gz`` (paired) or
``<sample>.fastq.gz`` (single). Uncompressed intermediates are written to a
private work directory and removed after successful compression; the
downloaded SRA file is never modified.
"""

from __future__ import annotations

import contextlib
import os
from pathlib import Path

from .context import StageContext, StageResult
from .errors import IntegrityError, ToolExecutionError
from .fastq_validation import validate_fastq
from .runner import run_command
from .utils import free_bytes, human_bytes, remove_if_exists, safe_rmtree
from .validation import Sample

_SPACE_HINTS = ("no space left", "disk-limit", "disk limit", "insufficient", "storage exhausted")


def build_command(sra_file: Path, work_dir: Path, tmp_dir: Path, threads: int) -> list[str]:
    return [
        "fasterq-dump",
        str(sra_file),
        "--outdir",
        str(work_dir),
        "--temp",
        str(tmp_dir),
        "--threads",
        str(max(1, threads)),
        "--split-3",
        "--skip-technical",
    ]


def classify_outputs(work_dir: Path) -> tuple[Path | None, Path | None, Path | None]:
    """Return (mate 1, mate 2, unsuffixed) FASTQ files produced by fasterq-dump."""
    r1 = r2 = single = None
    extra: list[str] = []
    for p in sorted(work_dir.glob("*.fastq")):
        if p.name.endswith("_1.fastq"):
            r1 = p
        elif p.name.endswith("_2.fastq"):
            r2 = p
        elif p.stem.rsplit("_", 1)[-1].isdigit() and "_" in p.stem:
            extra.append(p.name)
        else:
            single = p
    if extra:
        raise IntegrityError(
            f"fasterq-dump produced unexpected additional read files ({', '.join(extra)}); "
            "runs with more than two biological reads per spot are not supported",
        )
    return r1, r2, single


def _check_space(sample: Sample, sra_size: int, ctx: StageContext, result: StageResult) -> None:
    factor = float(ctx.config["conversion.space_factor"])
    need = int(sra_size * factor)
    out_free = free_bytes(ctx.layout.raw_fastq_dir)
    tmp_free = free_bytes(ctx.tmp_root)
    same_fs = False
    with contextlib.suppress(OSError):
        same_fs = os.stat(ctx.layout.raw_fastq_dir).st_dev == os.stat(ctx.tmp_root).st_dev
    required_out = need * 2 if same_fs else need
    for where, free, required in (
        ("output directory", out_free, required_out),
        ("temporary directory", None if same_fs else tmp_free, need),
    ):
        if free is not None and free < required:
            msg = (
                f"{where} has {human_bytes(free)} free; a rule-of-thumb estimate for converting this "
                f"{human_bytes(sra_size)} SRA file is ~{human_bytes(required)} "
                f"(conversion.space_factor={factor:g}). Conversion may fail; consider --tmpdir on "
                "another disk or freeing space."
            )
            result.warnings.append(msg)
            ctx.console.warning(msg, sample=sample.name)


def _compress(src: Path, dest: Path, ctx: StageContext, sample: Sample, log: Path) -> None:
    partial = dest.with_name(dest.name + ".partial")
    remove_if_exists(partial)
    level = int(ctx.config["conversion.compression_level"])
    run_command(
        ["pigz", "-p", str(ctx.threads), f"-{level}", "-c", str(src)],
        log_path=log,
        stdout_path=partial,
        timeout=float(ctx.config["conversion.timeout_hours"]) * 3600,
        description=f"pigz {src.name}",
    )
    if partial.stat().st_size == 0:
        raise IntegrityError(f"compression produced an empty file: {partial}")
    os.replace(partial, dest)


def run(sample: Sample, sra_file: Path, ctx: StageContext) -> StageResult:
    result = StageResult()
    layout = ctx.layout
    layout.raw_fastq_dir.mkdir(parents=True, exist_ok=True)
    work = layout.conversion_work_dir(sample)
    tmp = layout.tmp_for(ctx.tmp_root, sample, "fasterq-dump")
    # Both directories are private to this sample and stage: start clean.
    safe_rmtree(work, layout.root)
    if tmp.exists():
        safe_rmtree(tmp, ctx.tmp_root)
    work.mkdir(parents=True)
    tmp.mkdir(parents=True)
    log = layout.command_log(sample, "convert")

    _check_space(sample, sra_file.stat().st_size, ctx, result)
    try:
        run_command(
            build_command(sra_file, work, tmp, ctx.threads),
            log_path=log,
            timeout=float(ctx.config["conversion.timeout_hours"]) * 3600,
            description=f"fasterq-dump {sample.accession}",
        )
    except ToolExecutionError as exc:
        low = exc.output_tail.lower()
        if any(h in low for h in _SPACE_HINTS):
            exc.hint = (
                "fasterq-dump ran out of disk space. Free space, or point --tmpdir to a larger "
                f"disk, then re-run with --resume. Log: {log}"
            )
        else:
            exc.hint = (
                f"See {log}. A corrupt or incomplete SRA download can cause this; delete "
                f"{layout.sra_accession_dir(sample)} and re-run with --resume to download again."
            )
        raise

    r1, r2, single = classify_outputs(work)
    if r1 is None and r2 is None and single is None:
        raise IntegrityError(
            f"fasterq-dump finished but produced no FASTQ files for {sample.accession}",
            hint=f"The run may contain no biological reads. See {log}.",
        )
    data_paired = r1 is not None and r2 is not None
    if (r1 is None) != (r2 is None):
        raise IntegrityError(
            f"fasterq-dump produced only one mate file for {sample.accession}; cannot form pairs",
            hint=f"See {log}; this run may have an unusual read structure.",
        )
    if sample.paired and not data_paired:
        raise IntegrityError(
            f"{sample.name} was declared paired-end, but {sample.accession} contains single-end "
            "reads only",
            hint="Set layout to 'single' for this sample (check LibraryLayout in the SRA record).",
        )
    if not sample.paired and data_paired:
        raise IntegrityError(
            f"{sample.name} was declared single-end, but {sample.accession} contains paired reads",
            hint="Set layout to 'paired' for this sample (check LibraryLayout in the SRA record).",
        )

    finals = layout.raw_fastq(sample)
    sources = [r1, r2] if data_paired else [single]
    for src, dest in zip(sources, finals, strict=True):
        assert src is not None
        if src.stat().st_size == 0:
            raise IntegrityError(f"fasterq-dump wrote an empty file: {src.name}")
        _compress(src, dest, ctx, sample, log)
        result.outputs["R1" if dest is finals[0] else "R2"] = dest
        remove_if_exists(src)
    if not sample.paired:
        result.outputs = {"reads": finals[0]}

    if data_paired and single is not None:
        unpaired_reads = validate_fastq(single).records
        result.metrics["unpaired_reads"] = unpaired_reads
        if ctx.config["conversion.keep_unpaired"]:
            layout.unpaired_dir.mkdir(parents=True, exist_ok=True)
            dest = layout.unpaired_fastq(sample)
            _compress(single, dest, ctx, sample, log)
            result.outputs["unpaired"] = dest
            where = f"kept separately in {layout.rel(dest)}"
        else:
            where = "discarded (conversion.keep_unpaired is false)"
        result.notes.append(
            f"{unpaired_reads} reads without a mate were {where}; they are not part of R1/R2."
        )
        remove_if_exists(single)
    else:
        result.metrics["unpaired_reads"] = 0

    safe_rmtree(work, layout.root)
    safe_rmtree(tmp, ctx.tmp_root)
    for label, p in result.outputs.items():
        result.metrics[f"{label}_bytes"] = p.stat().st_size
    result.notes.append(
        "Converted to "
        + ", ".join(layout.rel(p) for k, p in result.outputs.items() if k != "unpaired")
    )
    return result
