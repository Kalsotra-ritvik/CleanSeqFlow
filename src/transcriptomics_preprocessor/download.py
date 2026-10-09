"""Stage 1: download SRA runs with ``prefetch`` (SRA Toolkit).

An SRA file is NCBI's compressed archive format for a sequencing run; it is
not FASTQ. ``prefetch`` downloads it (resuming partial downloads and
verifying the result against NCBI's checksum), and ``fasterq-dump`` later
converts it into FASTQ read files.
"""

from __future__ import annotations

from pathlib import Path

from .context import StageContext, StageResult
from .errors import IntegrityError, InterruptedRun, ToolExecutionError
from .runner import STOP_EVENT, run_command
from .utils import free_bytes, human_bytes, parse_size
from .validation import Sample

#: Substrings (lower case) in prefetch output indicating a network/connectivity
#: problem. Checked FIRST: e.g. sra-tools 3.4.1 reports an unreachable NCBI as
#: "connection not found ... cannot resolve remote location of '<acc>'", which
#: must not be mistaken for an invalid accession. (Observed on 2026-10-09 with
#: NCBI blocked by a firewall.)
_NETWORK_ERROR_MARKERS = (
    "connection not found",
    "network system module",
    "timeout",
    "timed out",
    "connection refused",
    "connection reset",
    "could not resolve host",
    "couldn't resolve host",
    "transfer incomplete",
    "ssl",
    "tls",
)
#: Substrings (lower case) suggesting a permanent problem (not retried). Matching
#: is best-effort: an unrecognised failure is treated as transient and retried.
_PERMANENT_ERROR_MARKERS = (
    "invalid accession",
    "is not a valid",
    "access denied",
    "authorization required",
    "requires authorization",
    "dbgap",
    "larger than maximum allowed",
    "exceeds maximum",
    "item not found",
)
_NETWORK_HINT = (
    "Check the internet connection and that NCBI (https://www.ncbi.nlm.nih.gov) is reachable; "
    "then re-run with --resume (prefetch continues partial downloads)."
)


def build_command(sample: Sample, ctx: StageContext) -> list[str]:
    max_kb = parse_size(ctx.config["download.max_size"]) // 1024
    return [
        "prefetch",
        sample.accession,
        "--output-directory",
        str(ctx.layout.sra_dir),
        "--max-size",
        str(max_kb),
        "--verify",
        "yes",
        "--resume",
        "yes",
    ]


def find_sra_file(sample: Sample, ctx: StageContext) -> Path | None:
    """Locate the archive prefetch produced (``.sra`` or SRA-lite ``.sralite``)."""
    d = ctx.layout.sra_accession_dir(sample)
    for name in (f"{sample.accession}.sra", f"{sample.accession}.sralite"):
        p = d / name
        if p.is_file():
            return p
    return None


def classify_failure(output: str) -> str:
    """Return ``'network'``, ``'permanent'`` or ``'unknown'`` for prefetch output."""
    low = output.lower()
    if any(marker in low for marker in _NETWORK_ERROR_MARKERS):
        return "network"
    if any(marker in low for marker in _PERMANENT_ERROR_MARKERS):
        return "permanent"
    return "unknown"


def run(sample: Sample, ctx: StageContext) -> StageResult:
    result = StageResult()
    free = free_bytes(ctx.layout.sra_dir)
    min_free = parse_size(ctx.config["download.min_free_space"])
    if free is not None and free < min_free:
        result.warnings.append(
            f"only {human_bytes(free)} free at {ctx.layout.sra_dir} (warning level "
            f"{human_bytes(min_free)}). SRA files plus decompressed FASTQ and temporary files "
            "can need several times the SRA size."
        )
        ctx.console.warning(result.warnings[-1], sample=sample.name)

    log = ctx.layout.command_log(sample, "download")
    cmd = build_command(sample, ctx)
    attempts = int(ctx.config["download.retries"]) + 1
    wait = float(ctx.config["download.retry_wait_seconds"])
    timeout = float(ctx.config["download.timeout_hours"]) * 3600
    last_error: ToolExecutionError | None = None
    for attempt in range(1, attempts + 1):
        try:
            run_command(
                cmd,
                log_path=log,
                timeout=timeout,
                description=f"prefetch {sample.accession} (attempt {attempt}/{attempts})",
            )
            last_error = None
            break
        except ToolExecutionError as exc:
            last_error = exc
            kind = classify_failure(exc.output_tail)
            if kind == "permanent":
                raise ToolExecutionError(
                    f"prefetch could not download {sample.accession}: the accession appears to be "
                    "invalid, unavailable, too large for --max-size, or access-controlled",
                    hint=f"Check the accession on https://www.ncbi.nlm.nih.gov/sra/{sample.accession} "
                    f"and the log {log}. Raise download.max_size if the run is larger.",
                    returncode=exc.returncode,
                    log_path=str(log),
                    output_tail=exc.output_tail,
                ) from exc
            if attempt < attempts:
                delay = wait * (2 ** (attempt - 1))
                why = "network problem" if kind == "network" else exc.message
                ctx.console.warning(
                    f"download attempt {attempt}/{attempts} failed ({why}); "
                    f"retrying in {delay:.0f} s",
                    sample=sample.name,
                )
                if STOP_EVENT.wait(delay):
                    raise InterruptedRun("interrupted while waiting to retry a download") from exc
    if last_error is not None:
        kind = classify_failure(last_error.output_tail)
        reason = (
            "NCBI could not be reached (network problem)"
            if kind == "network"
            else "the cause could not be determined automatically"
        )
        raise ToolExecutionError(
            f"prefetch failed for {sample.accession} after {attempts} attempt(s); {reason}",
            hint=_NETWORK_HINT
            + f" If the network is fine, check that {sample.accession} exists and is public: "
            f"https://www.ncbi.nlm.nih.gov/sra/{sample.accession}. Details: {log}",
            returncode=last_error.returncode,
            log_path=str(log),
            output_tail=last_error.output_tail,
        )

    sra = find_sra_file(sample, ctx)
    if sra is None:
        raise IntegrityError(
            f"prefetch reported success but no {sample.accession}.sra/.sralite file was found in "
            f"{ctx.layout.sra_accession_dir(sample)}",
            hint=f"Inspect {log}. Some runs are only available from cloud providers or require "
            "authorisation.",
        )
    size = sra.stat().st_size
    if size == 0:
        raise IntegrityError(f"downloaded SRA file is empty: {sra}")
    result.outputs["sra"] = sra
    result.metrics.update(
        {"sra_file": sra.name, "sra_bytes": size, "sra_lite": sra.suffix == ".sralite"}
    )
    result.notes.append(f"Downloaded {sra.name} ({human_bytes(size)}); prefetch verified it.")
    if sra.suffix == ".sralite":
        result.warnings.append(
            "NCBI delivered an SRA-lite file: base qualities are simplified (not the original "
            "per-base values), so quality-based trimming and quality plots are less informative."
        )
    return result
