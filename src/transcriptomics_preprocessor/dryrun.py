"""Dry-run mode: validate everything and show the plan without touching any data.

A dry run never downloads, never runs FastQC/fastp/MultiQC, never creates the
output directory and never writes the run journal.
"""

from __future__ import annotations

import os
import shlex
from pathlib import Path
from typing import TYPE_CHECKING

from . import conversion, download, fastp, fastqc, multiqc
from .context import StageContext
from .dependencies import detect_all
from .errors import ExitCode, StateError
from .explain import STAGE_LABELS
from .state import SAMPLE_STAGES, RunState
from .utils import free_bytes, human_bytes
from .validation import check_output_path

if TYPE_CHECKING:
    from .workflow import Pipeline


def _nearest_existing(path: Path) -> Path:
    p = path
    while not p.exists() and p.parent != p:
        p = p.parent
    return p


def dry_run(p: Pipeline) -> int:
    c = p.console
    L = p.layout
    cfg = p.config
    out: list[str] = []
    add = out.append
    code = ExitCode.OK

    add(c.style("DRY RUN — nothing will be downloaded, processed, created or modified.", "bold"))
    add("")
    add("Validated configuration:")
    add(f"  Output directory : {L.root}")
    add(f"  Temporary dir    : {p.tmp_root}")
    add(
        f"  Threads          : {cfg.threads} total; {cfg.jobs} sample(s) at a time; "
        f"{cfg.threads_per_job} thread(s) per sample"
    )
    add(f"  Parallel downloads: up to {min(cfg.jobs, int(cfg['download.max_parallel']))}")
    add(f"  On sample failure: {cfg['on_error']}")
    add(f"  FASTQ validation : {cfg['validation.mode']}")
    add(
        f"  Keep SRA files   : {'yes' if cfg['download.keep_sra'] else 'no (removed after success)'}"
    )
    add(f"  Execution mode   : {p.opts.mode}")

    # Output directory checks (read-only).
    problems = check_output_path(p.opts.outdir)
    root = L.root
    if root.exists() and any(root.iterdir()) and not L.marker_file.exists():
        problems.append(f"{root} is not empty and was not created by this tool")
    has_state = L.state_file.exists()
    if has_state and p.opts.mode == "new":
        problems.append(
            f"{root} already contains a previous run; a real run needs --resume or --overwrite"
        )
    parent = _nearest_existing(root)
    if not os.access(parent, os.W_OK):
        problems.append(f"no write permission in {parent}")
    if problems:
        code = ExitCode.USAGE
        add("")
        add(c.style("A real run would stop with these problem(s):", "red"))
        add("".join(f"  - {m}\n" for m in problems).rstrip())

    add("")
    add("Disk space:")
    for label, path in (("output", root), ("temporary", p.tmp_root)):
        add(f"  {label:9s}: {human_bytes(free_bytes(path))} free at {_nearest_existing(path)}")
    add(
        "  Note: SRA files, decompressed FASTQ during conversion, compressed raw FASTQ and "
        "cleaned FASTQ all need space; plan for several times the SRA size per sample."
    )

    add("")
    add("External tools:")
    tools = detect_all()
    p.tools = tools
    missing = []
    for key, st in tools.items():
        if st.ok:
            pin = "" if st.version_matches else f" (tested: {st.expected})"
            add(f"  {key:13s} {st.version}{pin}  {st.path}")
        else:
            missing.append(key)
            add(f"  {key:13s} MISSING — {st.error}")
    if missing:
        add(
            c.style(
                f"  A real run would stop: missing {', '.join(missing)}. Run ./install.sh.", "red"
            )
        )
        code = code or ExitCode.DEPENDENCY

    add("")
    add(f"Samples ({len(p.opts.samples)}):")
    for s in p.opts.samples:
        meta = ", ".join(f"{k}={v}" for k, v in s.metadata.items())
        add(f"  {s.name:20s} {s.accession:14s} {s.layout:6s} {meta}")

    state = None
    if has_state and p.opts.mode == "resume":
        try:
            state = RunState.load(L.state_file)
        except StateError as exc:
            add(c.style(f"Run journal unreadable: {exc.message}", "red"))
            code = code or ExitCode.FAILED
    p.state = state if state is not None else RunState.new(L.state_file)

    ctx = StageContext(L, cfg, p.tmp_root, c, p.tool_versions, cfg.threads_per_job)
    add("")
    add("Planned actions:")
    for s in p.opts.samples:
        add(f"  Sample {s.name}:")
        for stage in SAMPLE_STAGES:
            if state is not None:
                ok, why = p._reusable(s, stage)
                verdict = f"reuse ({why})" if ok else f"run ({why})"
            else:
                verdict = "run"
            add(f"    - {STAGE_LABELS[stage]}: {verdict}")
        sra = L.sra_accession_dir(s) / f"{s.accession}.sra"
        raw, clean = L.raw_fastq(s), L.clean_fastq(s)
        cmds = [
            download.build_command(s, ctx),
            conversion.build_command(
                sra,
                L.conversion_work_dir(s),
                L.tmp_for(p.tmp_root, s, "fasterq-dump"),
                cfg.threads_per_job,
            ),
            fastqc.build_command(
                raw,
                L.raw_fastqc_dir(s),
                L.tmp_for(p.tmp_root, s, "fastqc_raw"),
                cfg.threads_per_job,
                int(cfg["fastqc.memory_mb"]),
            ),
            fastp.build_command(
                s, raw, clean, L.fastp_json(s), L.fastp_html(s), cfg, cfg.threads_per_job
            ),
            fastqc.build_command(
                clean,
                L.clean_fastqc_dir(s),
                L.tmp_for(p.tmp_root, s, "fastqc_clean"),
                cfg.threads_per_job,
                int(cfg["fastqc.memory_mb"]),
            ),
        ]
        add(
            "    Commands that would be executed (paths are final names; temporary names may differ):"
        )
        add("".join(f"      $ {shlex.join(cmd)}\n" for cmd in cmds).rstrip())
    mq = multiqc.build_command(
        L.raw_qc_dir / "multiqc" / "raw_multiqc_report_inputs.txt",
        L.raw_qc_dir / "multiqc",
        "raw_multiqc_report",
        "Raw reads QC (before fastp)",
    )
    add("  Then, across all samples:")
    add(f"      $ {shlex.join(mq)}")
    add("      (and the same for the cleaned reads + fastp reports), then the summary report.")
    add("")
    add("fastp parameters:")
    add(
        "".join(
            f"  - {line}\n"
            for line in fastp.describe_parameters(cfg, any(s.paired for s in p.opts.samples))
        ).rstrip()
    )
    add("")
    add(
        "Dry run finished: input and configuration are valid."
        if code == ExitCode.OK
        else "Dry run finished: problems were found (see above)."
    )
    c.summary("\n".join(out))
    return int(code)
