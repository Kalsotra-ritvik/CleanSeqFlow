"""Command-line interface: ``transcriptomics-preprocess``."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from . import __version__
from .config import build_config, cpu_count
from .console import Console
from .errors import ExitCode, InputValidationError, PreprocessorError
from .explain import render_explanation
from .validation import Sample, read_metadata_csv, samples_from_accessions

PROG = "transcriptomics-preprocess"

EPILOG = f"""\
examples:
  {PROG} doctor
  {PROG} check-input --metadata metadata.csv
  {PROG} run --accessions SRR26053372 SRR26053373 --layout paired --outdir RNAseq_results --threads 8
  {PROG} run --metadata metadata.csv --outdir RNAseq_results --threads 8
  {PROG} run --metadata metadata.csv --outdir RNAseq_results --dry-run
  {PROG} run --metadata metadata.csv --outdir RNAseq_results --resume
  {PROG} explain qc

exit codes: 0 success, 2 invalid input/usage, 3 missing dependency, 4 partial failure,
            5 failure, 130 interrupted.
This tool performs preprocessing and QC only (no alignment or quantification).
"""


class _Formatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def _positive_int(text: str) -> int:
    try:
        v = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a whole number, got {text!r}") from None
    if v < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1, got {v}")
    return v


def _nonneg_int(text: str) -> int:
    try:
        v = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a whole number, got {text!r}") from None
    if v < 0:
        raise argparse.ArgumentTypeError(f"must be 0 or more, got {v}")
    return v


def _add_input_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("input (choose one)")
    src = g.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--metadata",
        "-m",
        type=Path,
        metavar="CSV",
        help="sample sheet with columns sample, accession, layout (+ optional extra columns)",
    )
    src.add_argument(
        "--accessions",
        "-a",
        nargs="+",
        metavar="SRR",
        help="one or more SRA run accessions (SRR/ERR/DRR); each becomes a sample of the same name",
    )
    g.add_argument(
        "--layout",
        choices=["paired", "single"],
        help="sequencing layout; required with --accessions, default for empty CSV cells",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROG,
        description="Automated SRA-to-clean-FASTQ preprocessing and quality control for bulk RNA-seq "
        "(SRA Toolkit → FastQC → MultiQC → fastp → FastQC → MultiQC).",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    # run -------------------------------------------------------------------
    run = sub.add_parser(
        "run",
        help="download, convert, QC and clean reads",
        formatter_class=_Formatter,
        epilog=EPILOG,
        description="Run the preprocessing pipeline. Without --resume/--overwrite "
        "the output directory must be new or empty.",
    )
    _add_input_args(run)
    o = run.add_argument_group("output and execution")
    o.add_argument("--outdir", "-o", type=Path, required=True, help="output directory")
    o.add_argument(
        "--config", "-c", type=Path, help="YAML configuration file (see config/default.yaml)"
    )
    o.add_argument(
        "--threads", "-t", type=_positive_int, help="total CPU threads (default from config: 4)"
    )
    o.add_argument(
        "--jobs", "-j", type=_positive_int, help="samples processed at the same time (default 1)"
    )
    o.add_argument(
        "--tmpdir", type=Path, help="directory for large temporary files (default <outdir>/tmp)"
    )
    o.add_argument("--project", help="project name shown in reports (default: output folder name)")
    o.add_argument(
        "--on-error",
        choices=["continue", "stop"],
        help="after a sample fails: continue with others (default) or stop",
    )
    mode = o.add_mutually_exclusive_group()
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="validate input and show planned actions; download and process nothing",
    )
    mode.add_argument(
        "--resume",
        action="store_true",
        help="continue a previous run in --outdir, reusing valid completed steps",
    )
    mode.add_argument(
        "--overwrite",
        action="store_true",
        help="recompute every step in an existing pipeline output directory",
    )
    o.add_argument(
        "--remove-sra",
        action="store_true",
        help="delete each downloaded .sra file after its sample completed successfully "
        "(raw FASTQ files are kept)",
    )
    o.add_argument(
        "--no-checksums", action="store_true", help="skip SHA-256 checksums of final files"
    )
    o.add_argument(
        "--validation",
        choices=["full", "quick"],
        help="FASTQ validation depth: full (every record, exact counts) or quick",
    )

    d = run.add_argument_group("download")
    d.add_argument(
        "--max-download-size", metavar="SIZE", help="largest SRA file to download, e.g. 50G"
    )
    d.add_argument("--retries", type=_nonneg_int, help="download retries after network failures")

    f = run.add_argument_group("fastp cleaning (defaults: fastp's own defaults)")
    f.add_argument(
        "--min-length", type=_positive_int, help="discard reads shorter than this (fastp -l)"
    )
    f.add_argument(
        "--qualified-quality",
        type=_nonneg_int,
        help="Phred quality for a base to count as qualified (fastp -q)",
    )
    f.add_argument(
        "--unqualified-percent",
        type=_nonneg_int,
        help="max %% of unqualified bases per read (fastp -u)",
    )
    f.add_argument("--n-base-limit", type=_nonneg_int, help="max N bases per read (fastp -n)")
    f.add_argument(
        "--average-quality", type=_nonneg_int, help="min mean read quality, 0 = off (fastp -e)"
    )
    f.add_argument("--adapter-r1", metavar="SEQ", help="adapter sequence for R1 / single-end reads")
    f.add_argument("--adapter-r2", metavar="SEQ", help="adapter sequence for R2")
    f.add_argument("--no-adapter-trimming", action="store_true", help="disable adapter trimming")
    f.add_argument("--no-quality-filter", action="store_true", help="disable quality filtering")
    f.add_argument(
        "--cut-front", action="store_true", help="sliding-window quality trimming from 5' end"
    )
    f.add_argument(
        "--cut-tail", action="store_true", help="sliding-window quality trimming from 3' end"
    )
    f.add_argument(
        "--cut-window-size", type=_positive_int, help="window size for --cut-front/--cut-tail"
    )
    f.add_argument(
        "--cut-mean-quality", type=_positive_int, help="mean quality for --cut-front/--cut-tail"
    )
    f.add_argument("--trim-poly-g", choices=["auto", "on", "off"], help="polyG tail trimming")

    out = run.add_argument_group("terminal output")
    verb = out.add_mutually_exclusive_group()
    verb.add_argument(
        "--quiet", "-q", action="store_true", help="only warnings, errors and the final summary"
    )
    verb.add_argument(
        "--verbose", "-v", action="store_true", help="also show command lines and details"
    )
    out.add_argument("--no-color", action="store_true", help="disable coloured output")

    # check-input -------------------------------------------------------------
    chk = sub.add_parser(
        "check-input",
        help="validate a sample sheet or accession list and exit",
        formatter_class=_Formatter,
    )
    _add_input_args(chk)
    chk.add_argument("--config", "-c", type=Path, help="also validate this configuration file")

    # doctor ------------------------------------------------------------------
    doc = sub.add_parser(
        "doctor",
        help="check the installation and environment",
        formatter_class=_Formatter,
        description="Check external programs, versions, Python environment, write "
        "permissions and free space.",
    )
    doc.add_argument("--outdir", type=Path, help="also check that this output location is usable")
    doc.add_argument("--tmpdir", type=Path, help="also check this temporary directory")
    doc.add_argument(
        "--smoke-test",
        action="store_true",
        help="run FastQC, fastp and MultiQC on a tiny synthetic dataset",
    )
    doc.add_argument("--json", action="store_true", help="print results as JSON")

    # explain -----------------------------------------------------------------
    exp = sub.add_parser("explain", help="explain the pipeline steps and QC metrics")
    exp.add_argument(
        "topic", nargs="?", choices=["steps", "qc", "limitations", "all"], default="all"
    )
    return parser


def _load_samples(args: argparse.Namespace) -> tuple[list[Sample], str]:
    if args.metadata is not None:
        return read_metadata_csv(args.metadata, args.layout), f"metadata file {args.metadata}"
    return samples_from_accessions(args.accessions, args.layout), "accessions " + " ".join(
        args.accessions
    )


def _overrides(args: argparse.Namespace) -> dict[str, Any]:
    o: dict[str, Any] = {
        "threads": args.threads,
        "jobs": args.jobs,
        "tmpdir": str(args.tmpdir) if args.tmpdir is not None else None,
        "on_error": args.on_error,
        "validation.mode": args.validation,
        "download.max_size": args.max_download_size,
        "download.retries": args.retries,
        "fastp.length_required": args.min_length,
        "fastp.qualified_quality_phred": args.qualified_quality,
        "fastp.unqualified_percent_limit": args.unqualified_percent,
        "fastp.n_base_limit": args.n_base_limit,
        "fastp.average_qual": args.average_quality,
        "fastp.adapter_sequence": args.adapter_r1,
        "fastp.adapter_sequence_r2": args.adapter_r2,
        "fastp.cut_window_size": args.cut_window_size,
        "fastp.cut_mean_quality": args.cut_mean_quality,
        "fastp.trim_poly_g": args.trim_poly_g,
    }
    if args.remove_sra:
        o["download.keep_sra"] = False
    if args.no_checksums:
        o["output.checksums"] = False
    if args.no_adapter_trimming:
        o["fastp.adapter_trimming"] = False
    if args.no_quality_filter:
        o["fastp.quality_filtering"] = False
    if args.cut_front:
        o["fastp.cut_front"] = True
    if args.cut_tail:
        o["fastp.cut_tail"] = True
    return o


def _report_error(console: Console, exc: PreprocessorError) -> int:
    console.error(exc.message, hint=exc.hint)
    for problem in getattr(exc, "problems", []):
        console.error(f"  - {problem}")
    return exc.exit_code


def cmd_run(args: argparse.Namespace, console: Console, argv: Sequence[str]) -> int:
    from .workflow import Pipeline, RunOptions

    samples, desc = _load_samples(args)
    config = build_config(args.config, _overrides(args))
    if config.threads > cpu_count():
        console.warning(
            f"--threads {config.threads} exceeds the {cpu_count()} CPUs available; "
            "this oversubscribes the machine"
        )
    mode = "resume" if args.resume else "overwrite" if args.overwrite else "new"
    opts = RunOptions(
        outdir=args.outdir,
        samples=samples,
        config=config,
        mode=mode,
        dry_run=args.dry_run,
        command_line=[PROG, *argv],
        input_description=desc,
        project=args.project,
    )
    return Pipeline(opts, console).run()


def cmd_check_input(args: argparse.Namespace, console: Console) -> int:
    samples, desc = _load_samples(args)
    if args.config is not None:
        build_config(args.config)
        console.summary(f"Configuration file {args.config}: valid")
    console.summary(f"Input ({desc}): {len(samples)} valid sample(s)")
    paired = sum(s.paired for s in samples)
    console.summary(f"  paired-end: {paired}, single-end: {len(samples) - paired}")
    for s in samples:
        meta = ", ".join(f"{k}={v}" for k, v in s.metadata.items())
        console.summary(f"  {s.name:20s} {s.accession:14s} {s.layout:6s} {meta}")
    console.summary(
        "Note: this checks the format of the input only; whether each accession exists "
        "at NCBI is checked when it is downloaded."
    )
    return ExitCode.OK


def _doctor_checks(args: argparse.Namespace) -> dict[str, Any]:
    from .dependencies import detect_all, environment_info
    from .utils import check_writable_dir, free_bytes, human_bytes

    tools = detect_all()
    checks: list[dict[str, Any]] = []
    for key, st in tools.items():
        if not st.ok:
            checks.append(
                {"check": f"tool {key}", "ok": False, "detail": st.error, "severity": "error"}
            )
        else:
            same = st.version_matches
            checks.append(
                {
                    "check": f"tool {key}",
                    "ok": True,
                    "detail": f"{st.version} at {st.path}"
                    + ("" if same else f" (tested version {st.expected})"),
                    "severity": "ok" if same else "warning",
                }
            )
    env = environment_info()
    py_ok = sys.version_info >= (3, 10)
    checks.append(
        {
            "check": "python",
            "ok": py_ok,
            "detail": f"{env['python']} ({env['python_executable']})",
            "severity": "ok" if py_ok else "error",
        }
    )
    try:
        import yaml

        checks.append(
            {
                "check": "python package PyYAML",
                "ok": True,
                "detail": yaml.__version__,
                "severity": "ok",
            }
        )
    except ImportError:  # pragma: no cover - import already needed by config
        checks.append(
            {
                "check": "python package PyYAML",
                "ok": False,
                "detail": "missing",
                "severity": "error",
            }
        )
    from .utils import is_within

    outside = sorted(
        st.path
        for st in tools.values()
        if st.path and not is_within(Path(st.path), Path(sys.prefix))
    )
    if outside:
        checks.append(
            {
                "check": "tool locations",
                "ok": True,
                "severity": "warning",
                "detail": f"{', '.join(outside)} not inside the active environment {sys.prefix}; "
                "use the launcher installed by install.sh for a consistent environment",
            }
        )
    for label, path in (("output directory", args.outdir), ("temporary directory", args.tmpdir)):
        if path is None:
            continue
        p = Path(path).expanduser()
        target = p
        while not target.exists() and target.parent != target:
            target = target.parent
        writable = os.access(target, os.W_OK)
        free = free_bytes(target)
        checks.append(
            {
                "check": f"{label} {p}",
                "ok": writable,
                "severity": "ok" if writable else "error",
                "detail": f"{'writable' if writable else 'NOT writable'} ({target}); free {human_bytes(free)}",
            }
        )
    tmp = Path(tempfile.gettempdir())
    err = check_writable_dir(tmp)
    checks.append(
        {
            "check": f"system temp {tmp}",
            "ok": err is None,
            "severity": "ok" if err is None else "warning",
            "detail": err or f"writable; free {human_bytes(free_bytes(tmp))}",
        }
    )
    return {
        "environment": env,
        "tools": {k: v.to_dict() for k, v in tools.items()},
        "checks": checks,
    }


def _smoke_test(console: Console) -> list[dict[str, Any]]:
    """Run the real QC tools on tiny synthetic paired and single-end data."""
    from .config import build_config
    from .context import StageContext
    from .layout import OutputLayout
    from .metrics import read_fastp_json, read_fastqc_zip
    from .runner import run_command
    from .synthetic import write_synthetic_fastq

    results: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="tpp-smoke-") as d:
        root = Path(d)
        r1, r2 = root / "smoke_R1.fastq.gz", root / "smoke_R2.fastq.gz"
        write_synthetic_fastq(r1, r2, n_reads=2000, accession="SMOKE1")
        cfg = build_config(overrides={"threads": 2})
        sample = Sample("smoke", "SRR0000001", "paired")
        ctx = StageContext(OutputLayout(root), cfg, root / "tmp", console, {}, 2)
        (root / "tmp").mkdir()

        def step(name: str, fn: Any) -> None:
            try:
                detail = fn()
                results.append(
                    {"check": f"smoke test: {name}", "ok": True, "detail": detail, "severity": "ok"}
                )
            except Exception as exc:
                results.append(
                    {
                        "check": f"smoke test: {name}",
                        "ok": False,
                        "severity": "error",
                        "detail": f"{exc}",
                    }
                )

        from . import fastp as fastp_mod
        from . import fastqc as fastqc_mod
        from . import multiqc as multiqc_mod

        def t_fastqc() -> str:
            res = fastqc_mod.run(sample, [r1, r2], root / "fqc", "fastqc_raw", ctx)
            m = read_fastqc_zip(res.outputs["R1_zip"])
            assert m.total_sequences == 2000, (
                f"FastQC counted {m.total_sequences} reads, expected 2000"
            )
            return "FastQC report parsed; 2000 reads counted"

        def t_fastp() -> str:
            res = fastp_mod.run(sample, [r1, r2], ctx)
            m = read_fastp_json(res.outputs["json"], True)
            assert m.reads_before == 4000, f"fastp counted {m.reads_before} reads, expected 4000"
            assert m.adapter_trimmed_reads, (
                "fastp trimmed no adapters from adapter-containing reads"
            )
            return f"fastp kept {m.reads_after} of 4000 reads; trimmed adapters from {m.adapter_trimmed_reads}"

        def t_multiqc() -> str:
            zips = sorted((root / "fqc").glob("*_fastqc.zip"))
            out = multiqc_mod.run(
                zips, root / "mqc", "smoke_multiqc", "smoke test", root / "mqc.log", 600, root
            )
            return f"report written ({out['html'].stat().st_size} bytes)"

        def t_sra() -> str:
            versions = []
            for exe in ("prefetch", "fasterq-dump"):
                r = run_command(
                    [exe, "--version"],
                    log_path=root / "sra.log",
                    timeout=60,
                    description=f"{exe} --version",
                )
                versions.append(f"{exe} exit {r.returncode}")
            return ", ".join(versions) + " (no network download attempted)"

        step("FastQC", t_fastqc)
        step("fastp", t_fastp)
        step("MultiQC", t_multiqc)
        step("SRA Toolkit executables", t_sra)
    return results


def cmd_doctor(args: argparse.Namespace, console: Console) -> int:
    report = _doctor_checks(args)
    if args.smoke_test:
        report["checks"] += _smoke_test(console)
    errors = [c for c in report["checks"] if not c["ok"]]
    warnings = [c for c in report["checks"] if c["ok"] and c["severity"] == "warning"]
    if args.json:
        print(json.dumps({**report, "ok": not errors}, indent=2, default=str))
    else:
        console.summary(f"{PROG} {__version__} — installation check")
        console.summary(f"Environment: {report['environment']['environment_prefix']}")
        for c in report["checks"]:
            mark = {
                "ok": console.style("OK  ", "green"),
                "warning": console.style("WARN", "yellow"),
                "error": console.style("FAIL", "red"),
            }[c["severity"] if c["ok"] else "error"]
            console.summary(f"  [{mark}] {c['check']}: {c['detail']}")
        if errors:
            console.summary(
                console.style(f"{len(errors)} problem(s) found.", "red")
                + " Re-run ./install.sh (or ./install.sh --repair) to fix the environment."
            )
        else:
            console.summary(
                console.style("All required checks passed.", "green")
                + (f" {len(warnings)} warning(s)." if warnings else "")
            )
    return ExitCode.DEPENDENCY if errors else ExitCode.OK


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    if not argv:
        parser.print_help()
        return ExitCode.USAGE
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:  # argparse already printed the message
        return int(exc.code or 0) if exc.code in (0, None) else ExitCode.USAGE
    if args.command is None:
        parser.print_help()
        return ExitCode.USAGE
    if args.command == "explain":
        sys.stdout.write(render_explanation(args.topic))
        return ExitCode.OK
    console = Console(
        quiet=getattr(args, "quiet", False),
        verbose=getattr(args, "verbose", False),
        color=False if getattr(args, "no_color", False) else None,
    )
    try:
        if args.command == "run":
            return cmd_run(args, console, argv)
        if args.command == "check-input":
            return cmd_check_input(args, console)
        if args.command == "doctor":
            return cmd_doctor(args, console)
    except InputValidationError as exc:
        return _report_error(console, exc)
    except PreprocessorError as exc:
        return _report_error(console, exc)
    except KeyboardInterrupt:
        console.error("interrupted")
        return ExitCode.INTERRUPTED
    raise AssertionError(f"unhandled command {args.command}")  # pragma: no cover


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
