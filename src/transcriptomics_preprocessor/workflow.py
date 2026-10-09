"""Pipeline orchestration: stage ordering, checkpoints, concurrency and failure policy."""

from __future__ import annotations

import os
import shutil
import signal
import threading
import time
import traceback
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from types import FrameType
from typing import Any

from . import __version__, conversion, download, fastp, fastqc, multiqc, reporting
from .comparison import build_comparison
from .config import PipelineConfig
from .console import Console
from .context import StageContext, StageResult
from .dependencies import TOOLS, ToolStatus, detect_all, environment_info, require
from .errors import (
    ExitCode,
    IntegrityError,
    InterruptedRun,
    OutputDirectoryError,
    PreprocessorError,
)
from .fastq_validation import validate_fastq, validate_fastq_pair
from .layout import OutputLayout
from .runner import STOP_EVENT, terminate_all
from .state import (
    AGGREGATE_STAGES,
    SAMPLE_STAGES,
    RunState,
    Status,
    can_reuse,
    exclusive_lock,
    outputs_intact,
    outputs_record,
)
from .utils import (
    atomic_write_json,
    atomic_write_text,
    check_writable_dir,
    file_signature,
    free_bytes,
    human_bytes,
    sha256_file,
    stable_hash,
    utc_now_iso,
)
from .validation import Sample, check_output_path

_READ_LABELS = ("R1", "R2", "reads")

#: Upstream dependencies of each per-sample stage: (stage, use outputs?, use fingerprint?).
_UPSTREAM: dict[str, list[tuple[str, bool, bool]]] = {
    "download": [],
    "convert": [("download", True, False)],
    "validate_raw": [("convert", True, False)],
    "fastqc_raw": [("convert", True, False)],
    "fastp": [("convert", True, False)],
    "validate_clean": [("fastp", True, False)],
    "fastqc_clean": [("fastp", True, False)],
    "publish": [("fastp", True, False), ("validate_clean", False, True)],
}
_STAGE_TOOLS: dict[str, list[str]] = {
    "download": ["prefetch"],
    "convert": ["fasterq-dump", "pigz"],
    "validate_raw": [],
    "fastqc_raw": ["fastqc"],
    "fastp": ["fastp"],
    "validate_clean": [],
    "fastqc_clean": ["fastqc"],
    "publish": [],
    "multiqc_raw": ["multiqc"],
    "multiqc_clean": ["multiqc"],
}


@dataclass
class RunOptions:
    """Everything the CLI collects for a run."""

    outdir: Path
    samples: list[Sample]
    config: PipelineConfig
    mode: str = "new"  # new | resume | overwrite
    dry_run: bool = False
    command_line: list[str] = field(default_factory=list)
    input_description: str = ""
    project: str | None = None


def _labels(outputs: dict[str, Any]) -> list[str]:
    return [k for k in _READ_LABELS if k in outputs]


class Pipeline:
    def __init__(self, opts: RunOptions, console: Console) -> None:
        self.opts = opts
        self.console = console
        self.config = opts.config
        self.layout = OutputLayout(opts.outdir.expanduser().absolute())
        tmp = self.config["tmpdir"]
        self.tmp_root = Path(tmp).expanduser().absolute() if tmp else self.layout.default_tmp_dir
        self.tools: dict[str, ToolStatus] = {}
        self.state: RunState | None = None
        self._abort = threading.Event()
        self._download_slots = threading.Semaphore(
            max(1, min(self.config.jobs, int(self.config["download.max_parallel"])))
        )
        self.started = utc_now_iso()

    # ------------------------------------------------------------------ setup
    @property
    def tool_versions(self) -> dict[str, str | None]:
        return {k: (s.version if s else None) for k, s in self.tools.items()}

    def _prepare_outdir(self) -> None:
        problems = check_output_path(self.opts.outdir)
        if problems:
            raise OutputDirectoryError("; ".join(problems))
        root = self.layout.root
        if root.exists() and any(root.iterdir()) and not self.layout.marker_file.exists():
            raise OutputDirectoryError(
                f"{root} is not empty and was not created by this tool",
                hint="Choose a new or empty --outdir. Existing files are never mixed with "
                "pipeline output, and nothing outside pipeline-managed folders is deleted.",
            )
        has_state = self.layout.state_file.exists()
        if has_state and self.opts.mode == "new":
            raise OutputDirectoryError(
                f"{root} already contains results from a previous run",
                hint="Use --resume to continue it (completed, still-valid steps are reused), "
                "--overwrite to recompute every step, or a different --outdir.",
            )
        err = check_writable_dir(root)
        if err:
            raise OutputDirectoryError(err)
        self.layout.marker_file.write_text(
            f"Output directory of transcriptomics-preprocess {__version__}. Do not delete this "
            "file: it marks the folder as safe for the pipeline to manage.\n",
            encoding="utf-8",
        )
        for d in (self.layout.logs_dir, self.layout.info_dir, self.layout.reports_dir):
            d.mkdir(parents=True, exist_ok=True)
        err = check_writable_dir(self.tmp_root)
        if err:
            raise OutputDirectoryError(f"temporary directory unusable: {err}")

    def _load_state(self) -> RunState:
        path = self.layout.state_file
        if self.opts.mode == "overwrite" and path.exists():
            backup = path.with_name(f"state.before-overwrite.{int(time.time())}.json")
            os.replace(path, backup)
            self.console.info(
                f"Previous run journal moved to {self.layout.rel(backup)}; "
                "every step will be recomputed."
            )
        if path.exists():
            state = RunState.load(path)
            self.console.info(
                "Resuming: completed steps are reused only if their parameters, "
                "tool versions, inputs and outputs are unchanged."
            )
        else:
            if self.opts.mode == "resume":
                self.console.info("No previous run found in this directory; starting a new run.")
            state = RunState.new(path)
        for s in self.opts.samples:
            entry = state.sample_entry(s.name)
            entry.update({"accession": s.accession, "layout": s.layout, "metadata": s.metadata})
        state.data["runs"].append(
            {
                "started": self.started,
                "mode": self.opts.mode,
                "command": self.opts.command_line,
                "version": __version__,
            }
        )
        state.save()
        return state

    def _write_run_info(self) -> None:
        info = self.layout.info_dir
        atomic_write_text(info / "config_used.yaml", self.config.to_yaml())
        atomic_write_json(
            info / "software_versions.json",
            {
                "environment": environment_info(),
                "tools": {k: v.to_dict() for k, v in self.tools.items()},
            },
        )
        atomic_write_text(info / "command_line.txt", " ".join(self.opts.command_line) + "\n")
        rows = ["sample,accession,layout"]
        rows += [f"{s.name},{s.accession},{s.layout}" for s in self.opts.samples]
        atomic_write_text(info / "samples_input.csv", "\n".join(rows) + "\n")

    # ------------------------------------------------------------ fingerprints
    def _stage_params(self, sample: Sample, stage: str) -> dict[str, Any]:
        cfg = self.config
        if stage == "download":
            return {"accession": sample.accession}
        if stage == "convert":
            return {
                "layout": sample.layout,
                "keep_unpaired": cfg["conversion.keep_unpaired"],
                "fasterq_dump_mode": ["--split-3", "--skip-technical"],
            }
        if stage in ("validate_raw", "validate_clean"):
            mode = cfg["validation.mode"]
            return {
                "mode": mode,
                "quick_records": cfg["validation.quick_records"] if mode == "quick" else None,
            }
        if stage == "fastp":
            p = fastp.stage_parameters(cfg, sample.paired)
            p.pop("compression_level", None)
            return p
        if stage == "publish":
            return {"checksums": cfg["output.checksums"], "keep_sra": cfg["download.keep_sra"]}
        return {}

    def _fingerprint(self, sample: Sample, stage: str) -> str:
        assert self.state is not None
        upstream = {}
        for up, use_outputs, use_fp in _UPSTREAM[stage]:
            rec = self.state.stage(sample.name, up)
            entry: dict[str, Any] = {}
            if use_outputs:
                entry["outputs"] = rec.get("outputs")
            if use_fp:
                entry["fingerprint"] = rec.get("fingerprint")
            upstream[up] = entry
        return stable_hash(
            {
                "stage": stage,
                "params": self._stage_params(sample, stage),
                "tools": {t: self.tool_versions.get(t) for t in _STAGE_TOOLS[stage]},
                "upstream": upstream,
            }
        )

    def _reusable(self, sample: Sample, stage: str) -> tuple[bool, str]:
        assert self.state is not None
        rec = self.state.stage(sample.name, stage)
        fp = self._fingerprint(sample, stage)
        if stage == "download" and rec.get("removed_by_cleanup"):
            # The SRA file was deleted on request after a successful conversion. It is only
            # needed again if the conversion itself has to be redone.
            if rec.get("status") == Status.COMPLETED.value and rec.get("fingerprint") == fp:
                conv_ok, _ = self._reusable(sample, "convert")
                if conv_ok:
                    return (
                        True,
                        "SRA file removed after conversion (cleanup option); conversion still valid",
                    )
            return False, "SRA file was removed and is needed again"
        return can_reuse(rec, fp, self.layout.root)

    # ------------------------------------------------------------- execution
    def _context(self) -> StageContext:
        return StageContext(
            self.layout,
            self.config,
            self.tmp_root,
            self.console,
            self.tool_versions,
            self.config.threads_per_job,
        )

    def _execute(
        self,
        sample: Sample,
        stage: str,
        position: tuple[int, int],
        func: Callable[[StageContext], StageResult],
        describe: Callable[[], tuple[list[str], list[str], list[str]]],
    ) -> bool:
        """Run one stage for one sample with checkpointing. Returns success."""
        assert self.state is not None
        ok, why = self._reusable(sample, stage)
        if ok:
            self.console.step_end(stage, sample.name, "REUSED", reason=f"valid checkpoint ({why})")
            return True
        rec_before = self.state.stage(sample.name, stage)
        if (
            rec_before.get("status") not in (None, Status.NOT_STARTED.value)
            and self.opts.mode == "resume"
        ):
            self.console.detail(f"[{sample.name}] {stage}: re-running because {why}")
        fp = self._fingerprint(sample, stage)
        tools = _STAGE_TOOLS[stage]
        tool = tools[0] if tools else None
        inputs, params, outputs = describe()
        self.console.step_start(
            stage,
            sample.name,
            position,
            inputs,
            tool,
            self.tool_versions.get(tool) if tool else None,
            params,
            outputs,
        )
        self.state.reset(sample.name, stage)
        self.state.update(
            sample.name, stage, status=Status.RUNNING, started=utc_now_iso(), fingerprint=None
        )
        t0 = time.monotonic()
        try:
            result = func(self._context())
        except InterruptedRun:
            self.state.update(
                sample.name,
                stage,
                status=Status.FAILED,
                message="interrupted before completion",
                finished=utc_now_iso(),
            )
            raise
        except PreprocessorError as exc:
            self._record_failure(sample, stage, exc.message, exc.hint, time.monotonic() - t0)
            return False
        except Exception as exc:  # isolate unexpected errors per sample, but record them loudly
            self.console.detail(traceback.format_exc())
            self._record_failure(
                sample,
                stage,
                f"unexpected internal error: {exc!r}",
                "This is probably a bug; the traceback is in logs/pipeline.log. Please report it.",
                time.monotonic() - t0,
            )
            return False
        duration = time.monotonic() - t0
        self.state.update(
            sample.name,
            stage,
            status=Status.COMPLETED,
            fingerprint=fp,
            finished=utc_now_iso(),
            duration_s=round(duration, 2),
            outputs=outputs_record(result.outputs, self.layout.root),
            metrics=result.metrics,
            warnings=result.warnings,
            notes=result.notes,
            tool_version=self.tool_versions.get(tool) if tool else None,
            message=None,
        )
        for w in result.warnings:
            if w not in [x["message"] for x in self.console.warnings]:
                self.console.warning(w, sample=sample.name)
        self.console.step_end(stage, sample.name, "COMPLETED", duration, result.notes)
        return True

    def _record_failure(
        self, sample: Sample, stage: str, message: str, hint: str | None, duration: float
    ) -> None:
        assert self.state is not None
        self.state.update(
            sample.name,
            stage,
            status=Status.FAILED,
            message=message,
            hint=hint,
            finished=utc_now_iso(),
            duration_s=round(duration, 2),
            fingerprint=None,
        )
        self.console.step_end(stage, sample.name, "FAILED", duration)
        self.console.error(message, hint=hint, sample=sample.name)

    def _skip(self, sample: Sample, stage: str, reason: str) -> None:
        assert self.state is not None
        self.state.update(
            sample.name,
            stage,
            status=Status.SKIPPED,
            message=reason,
            fingerprint=None,
            finished=utc_now_iso(),
        )
        self.console.step_end(stage, sample.name, "SKIPPED", reason=reason)

    def _outputs(self, sample: Sample, stage: str) -> dict[str, Path]:
        assert self.state is not None
        rec = self.state.stage(sample.name, stage)
        return {k: self.layout.abs(v["path"]) for k, v in (rec.get("outputs") or {}).items()}

    def _read_files(self, sample: Sample, stage: str) -> list[Path]:
        outs = self._outputs(sample, stage)
        return [outs[k] for k in _labels(outs)]

    # ------------------------------------------------------- stage functions
    def _stage_functions(
        self, sample: Sample
    ) -> dict[
        str,
        tuple[
            Callable[[StageContext], StageResult],
            Callable[[], tuple[list[str], list[str], list[str]]],
        ],
    ]:
        L = self.layout
        rel = L.rel
        cfg = self.config

        def do_download(ctx: StageContext) -> StageResult:
            with self._download_slots:
                return download.run(sample, ctx)

        def desc_download() -> tuple[list[str], list[str], list[str]]:
            return (
                [
                    f"SRA run accession {sample.accession} (sample '{sample.name}', {sample.layout}-end)"
                ],
                [
                    f"Maximum download size: {cfg['download.max_size']}",
                    f"Retries: {cfg['download.retries']}",
                    "Checksum verification: yes",
                ],
                [rel(L.sra_accession_dir(sample)) + "/"],
            )

        def do_convert(ctx: StageContext) -> StageResult:
            sra = self._outputs(sample, "download").get("sra")
            if sra is None or not sra.is_file():
                raise IntegrityError(
                    f"the downloaded SRA file for {sample.accession} is missing",
                    hint="Re-run with --resume to download it again.",
                )
            return conversion.run(sample, sra, ctx)

        def desc_convert() -> tuple[list[str], list[str], list[str]]:
            sra = self._outputs(sample, "download").get("sra")
            return (
                [rel(sra) if sra else "(downloaded SRA file)"],
                [
                    f"Layout: {sample.layout}-end",
                    "fasterq-dump --split-3 --skip-technical",
                    f"Threads: {cfg.threads_per_job}",
                    f"Temporary directory: {self.tmp_root}",
                ],
                [rel(p) for p in L.raw_fastq(sample)],
            )

        def make_validate(stage: str, source: str) -> Callable[[StageContext], StageResult]:
            def do_validate(ctx: StageContext) -> StageResult:
                files = self._read_files(sample, source)
                quick = cfg["validation.mode"] == "quick"
                limit = int(cfg["validation.quick_records"]) if quick else None
                res = StageResult()
                if sample.paired:
                    s1, s2 = validate_fastq_pair(files[0], files[1], max_records=limit)
                    res.metrics = {"R1": s1.to_dict(), "R2": s2.to_dict()}
                    if s1.records is not None:
                        res.notes.append(
                            f"{s1.records:,} read pairs ({2 * s1.records:,} reads); "
                            "R1/R2 counts and read names match."
                        )
                    else:
                        res.notes.append(
                            f"Quick check of the first {limit:,} pairs passed "
                            "(counts not determined in quick mode)."
                        )
                    for st in (s1, s2):
                        res.notes.append(
                            f"{Path(st.path).name}: lengths {st.min_length}-{st.max_length} bp, "
                            f"{human_bytes(st.file_size)}"
                        )
                else:
                    s = validate_fastq(files[0], max_records=limit)
                    res.metrics = {"reads": s.to_dict()}
                    res.notes.append(
                        f"{s.records:,} reads; lengths {s.min_length}-{s.max_length} bp"
                        if s.records is not None
                        else f"Quick check of the first {limit:,} reads passed."
                    )
                return res

            return do_validate

        def make_validate_desc(source: str) -> Callable[[], tuple[list[str], list[str], list[str]]]:
            def desc() -> tuple[list[str], list[str], list[str]]:
                files = self._read_files(sample, source)
                return (
                    [rel(f) for f in files],
                    [
                        f"Mode: {cfg['validation.mode']}"
                        + (" (all records)" if cfg["validation.mode"] == "full" else ""),
                        "Checks: FASTQ structure, sequence/quality characters, gzip integrity"
                        + (", R1/R2 pairing" if sample.paired else ""),
                    ],
                    [],
                )

            return desc

        def make_fastqc(
            stage: str, source: str, outdir: Path
        ) -> Callable[[StageContext], StageResult]:
            def do(ctx: StageContext) -> StageResult:
                return fastqc.run(sample, self._read_files(sample, source), outdir, stage, ctx)

            return do

        def make_fastqc_desc(
            source: str, outdir: Path
        ) -> Callable[[], tuple[list[str], list[str], list[str]]]:
            def desc() -> tuple[list[str], list[str], list[str]]:
                return (
                    [rel(f) for f in self._read_files(sample, source)],
                    [
                        f"Threads: {min(cfg.threads_per_job, 2 if sample.paired else 1)}",
                        f"Memory per file: {cfg['fastqc.memory_mb']} MB",
                    ],
                    [rel(outdir) + "/"],
                )

            return desc

        def do_fastp(ctx: StageContext) -> StageResult:
            return fastp.run(sample, self._read_files(sample, "convert"), ctx)

        def desc_fastp() -> tuple[list[str], list[str], list[str]]:
            params = [f"Threads: {min(cfg.threads_per_job, fastp.FASTP_MAX_THREADS)}"]
            params += fastp.describe_parameters(cfg, sample.paired)
            outs = [rel(p) for p in L.clean_fastq(sample)] + [
                rel(L.fastp_html(sample)),
                rel(L.fastp_json(sample)),
            ]
            return [rel(f) for f in self._read_files(sample, "convert")], params, outs

        def do_publish(ctx: StageContext) -> StageResult:
            return self._publish(sample)

        def desc_publish() -> tuple[list[str], list[str], list[str]]:
            params = [
                f"SHA-256 checksums: {'yes' if cfg['output.checksums'] else 'no'}",
                f"Keep SRA files: {'yes' if cfg['download.keep_sra'] else 'no (removed after success)'}",
            ]
            return (
                [rel(f) for f in self._read_files(sample, "fastp")],
                params,
                [rel(p) for p in L.final_fastq(sample)],
            )

        return {
            "download": (do_download, desc_download),
            "convert": (do_convert, desc_convert),
            "validate_raw": (
                make_validate("validate_raw", "convert"),
                make_validate_desc("convert"),
            ),
            "fastqc_raw": (
                make_fastqc("fastqc_raw", "convert", L.raw_fastqc_dir(sample)),
                make_fastqc_desc("convert", L.raw_fastqc_dir(sample)),
            ),
            "fastp": (do_fastp, desc_fastp),
            "validate_clean": (
                make_validate("validate_clean", "fastp"),
                make_validate_desc("fastp"),
            ),
            "fastqc_clean": (
                make_fastqc("fastqc_clean", "fastp", L.clean_fastqc_dir(sample)),
                make_fastqc_desc("fastp", L.clean_fastqc_dir(sample)),
            ),
            "publish": (do_publish, desc_publish),
        }

    def _publish(self, sample: Sample) -> StageResult:
        """Place validated cleaned files in the final directory (hard link, else copy)."""
        assert self.state is not None
        res = StageResult()
        L = self.layout
        L.final_dir.mkdir(parents=True, exist_ok=True)
        sources = self._read_files(sample, "fastp")
        clean_metrics = self.state.stage(sample.name, "validate_clean").get("metrics") or {}
        labels = ["R1", "R2"] if sample.paired else ["reads"]
        linked = 0
        for label, src, dest in zip(labels, sources, L.final_fastq(sample), strict=True):
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            try:
                os.link(src, dest)
                linked += 1
            except OSError:
                shutil.copy2(src, dest)
            res.outputs[label] = dest
            info: dict[str, Any] = {
                "file": L.rel(dest),
                "bytes": dest.stat().st_size,
                "records": (clean_metrics.get(label) or {}).get("records"),
            }
            if self.config["output.checksums"]:
                info["sha256"] = sha256_file(dest)
            res.metrics[label] = info
        res.notes.append(
            f"Final files: {', '.join(L.rel(p) for p in res.outputs.values())} "
            + (
                "(hard links to the fastp output; no extra disk space used)"
                if linked == len(labels)
                else "(copied)"
            )
        )
        if not self.config["download.keep_sra"]:
            sra = self._outputs(sample, "download").get("sra")
            if (
                sra is not None
                and sra.is_file()
                and sra.resolve().is_relative_to(L.sra_dir.resolve())
            ):
                size = sra.stat().st_size
                sra.unlink()
                self.state.update(sample.name, "download", removed_by_cleanup=True)
                res.notes.append(
                    f"Removed {L.rel(sra)} ({human_bytes(size)}) as requested "
                    "(--remove-sra); raw FASTQ files are kept."
                )
        return res

    def _process_sample(self, sample: Sample, index: int, total: int) -> None:
        funcs = self._stage_functions(sample)
        failed_stage: str | None = None
        for stage in SAMPLE_STAGES:
            if STOP_EVENT.is_set():
                raise InterruptedRun("interrupted")
            if failed_stage is not None:
                self._skip(sample, stage, f"not run because '{failed_stage}' failed")
                continue
            if self._abort.is_set():
                self._skip(sample, stage, "run stopped after a failure (on_error: stop)")
                continue
            func, describe = funcs[stage]
            if not self._execute(sample, stage, (index, total), func, describe):
                failed_stage = stage
                if self.config["on_error"] == "stop":
                    self._abort.set()

    # -------------------------------------------------------------- aggregates
    def _aggregate(self, stage: str, inputs: list[Path], filename: str, title: str) -> None:
        assert self.state is not None
        rec = self.state.stage(None, stage)
        if not inputs:
            self.state.update(
                None,
                stage,
                status=Status.SKIPPED,
                message="no sample reached this point",
                fingerprint=None,
                outputs={},
            )
            self.console.step_end(stage, None, "SKIPPED", reason="no input reports available")
            return
        fp = stable_hash(
            {
                "stage": stage,
                "tool": self.tool_versions.get("multiqc"),
                "inputs": [(self.layout.rel(p), file_signature(p)) for p in inputs],
            }
        )
        ok, _why = can_reuse(rec, fp, self.layout.root)
        if ok:
            self.console.step_end(stage, None, "REUSED", reason="inputs unchanged")
            return
        outdir = (
            self.layout.raw_qc_dir if stage == "multiqc_raw" else self.layout.clean_qc_dir
        ) / "multiqc"
        self.console.step_start(
            stage,
            None,
            None,
            [f"{len(inputs)} report file(s)"],
            "multiqc",
            self.tool_versions.get("multiqc"),
            [
                "Explicit input list (only samples of this run)",
                "Network features disabled (version check, AI, MegaQC upload)",
            ],
            [self.layout.rel(outdir / f"{filename}.html")],
        )
        self.state.update(
            None, stage, status=Status.RUNNING, started=utc_now_iso(), fingerprint=None
        )
        t0 = time.monotonic()
        try:
            outs = multiqc.run(
                inputs,
                outdir,
                filename,
                title,
                self.layout.command_log(None, stage),
                float(self.config["multiqc.timeout_hours"]) * 3600,
                self.layout.root,
            )
        except InterruptedRun:
            self.state.update(None, stage, status=Status.FAILED, message="interrupted")
            raise
        except PreprocessorError as exc:
            self.state.update(
                None,
                stage,
                status=Status.FAILED,
                message=exc.message,
                hint=exc.hint,
                finished=utc_now_iso(),
            )
            self.console.step_end(stage, None, "FAILED", time.monotonic() - t0)
            self.console.error(exc.message, hint=exc.hint)
            return
        dur = time.monotonic() - t0
        self.state.update(
            None,
            stage,
            status=Status.COMPLETED,
            fingerprint=fp,
            finished=utc_now_iso(),
            duration_s=round(dur, 2),
            outputs=outputs_record(outs, self.layout.root),
            message=None,
        )
        self.console.step_end(
            stage, None, "COMPLETED", dur, [f"Report: {self.layout.rel(outs['html'])}"]
        )

    def _completed(self, sample: Sample, stage: str) -> bool:
        assert self.state is not None
        return self.state.status(sample.name, stage) == Status.COMPLETED

    def _run_aggregates(self) -> None:
        raw_inputs: list[Path] = []
        clean_inputs: list[Path] = []
        for s in self.opts.samples:
            if self._completed(s, "fastqc_raw"):
                raw_inputs += [
                    p for k, p in self._outputs(s, "fastqc_raw").items() if k.endswith("_zip")
                ]
            if self._completed(s, "fastqc_clean"):
                clean_inputs += [
                    p for k, p in self._outputs(s, "fastqc_clean").items() if k.endswith("_zip")
                ]
            if self._completed(s, "fastp"):
                clean_inputs.append(self._outputs(s, "fastp")["json"])
        self._aggregate(
            "multiqc_raw", raw_inputs, "raw_multiqc_report", "Raw reads QC (before fastp)"
        )
        self._aggregate(
            "multiqc_clean", clean_inputs, "clean_multiqc_report", "Cleaned reads QC (after fastp)"
        )

    # ---------------------------------------------------------------- summary
    def build_summary(self, interrupted: bool) -> dict[str, Any]:
        assert self.state is not None
        samples_out = []
        counts = {"attempted": len(self.opts.samples), "completed": 0, "failed": 0, "skipped": 0}
        for s in self.opts.samples:
            entry = self.state.sample_entry(s.name)
            stages = {st: dict(self.state.stage(s.name, st)) for st in SAMPLE_STAGES}
            statuses = [rec.get("status") for rec in stages.values()]
            if stages["publish"].get("status") == Status.COMPLETED.value:
                overall = "COMPLETED"
                counts["completed"] += 1
            elif Status.FAILED.value in statuses:
                overall = "FAILED"
                counts["failed"] += 1
            else:
                overall = "INCOMPLETE" if interrupted else "SKIPPED"
                counts["skipped"] += 1
            cmp = None
            if stages["fastp"].get("status") == Status.COMPLETED.value:
                cmp = build_comparison(
                    s.name,
                    s.paired,
                    stages["validate_raw"].get("metrics"),
                    stages["validate_clean"].get("metrics")
                    if stages["validate_clean"].get("status") == "COMPLETED"
                    else None,
                    stages["fastp"].get("metrics"),
                    stages["fastqc_raw"].get("metrics"),
                    stages["fastqc_clean"].get("metrics")
                    if stages["fastqc_clean"].get("status") == "COMPLETED"
                    else None,
                    self.config,
                ).to_dict()
            failure = next(
                (
                    {"stage": st, "message": rec.get("message"), "hint": rec.get("hint")}
                    for st, rec in stages.items()
                    if rec.get("status") == Status.FAILED.value
                ),
                None,
            )
            samples_out.append(
                {
                    "sample": s.name,
                    "accession": s.accession,
                    "layout": s.layout,
                    "metadata": s.metadata,
                    "status": overall,
                    "failure": failure,
                    "stages": stages,
                    "comparison": cmp,
                    "final_files": stages["publish"].get("metrics")
                    if overall == "COMPLETED"
                    else None,
                    "removed_by_cleanup": entry.get("stages", {})
                    .get("download", {})
                    .get("removed_by_cleanup", False),
                }
            )
        aggregates = {
            st: dict(self.state.stage(None, st)) for st in AGGREGATE_STAGES if st != "report"
        }
        stage_warnings = [
            {"sample": s["sample"], "stage": st, "message": w}
            for s in samples_out
            for st, rec in s["stages"].items()
            for w in (rec.get("warnings") or [])
        ]
        return {
            "tool": "transcriptomics-preprocessor",
            "version": __version__,
            "project": self.opts.project or self.layout.root.name,
            "started": self.started,
            "finished": utc_now_iso(),
            "interrupted": interrupted,
            "mode": self.opts.mode,
            "command_line": self.opts.command_line,
            "input": self.opts.input_description,
            "outdir": str(self.layout.root),
            "tmpdir": str(self.tmp_root),
            "environment": environment_info(),
            "tools": {k: v.to_dict() for k, v in self.tools.items()},
            "config": self.config.as_nested(),
            "counts": counts,
            "samples": samples_out,
            "aggregates": aggregates,
            "warnings": stage_warnings,
            "console_warnings": self.console.warnings,
        }

    # ------------------------------------------------------------------- main
    def _install_signal_handlers(self) -> dict[int, Any]:
        previous: dict[int, Any] = {}

        def handler(signum: int, frame: FrameType | None) -> None:
            if STOP_EVENT.is_set():
                raise KeyboardInterrupt
            STOP_EVENT.set()
            self.console.warning(
                f"received {signal.Signals(signum).name}: stopping running tools and saving the "
                "run journal (press Ctrl-C again to force quit)"
            )

        if threading.current_thread() is threading.main_thread():
            for sig in (signal.SIGINT, signal.SIGTERM):
                previous[sig] = signal.signal(sig, handler)
        return previous

    def run(self) -> int:
        if self.opts.dry_run:
            from .dryrun import dry_run

            return dry_run(self)
        self._prepare_outdir()
        with exclusive_lock(self.layout.lock_file):
            self.console.attach_logs(
                self.layout.logs_dir / "pipeline.log", self.layout.logs_dir / "events.jsonl"
            )
            try:
                return self._run_locked()
            finally:
                self.console.close()

    def _run_locked(self) -> int:
        self.console.summary(
            f"transcriptomics-preprocess {__version__} — {len(self.opts.samples)} "
            f"sample(s) — output: {self.layout.root}"
        )
        self.tools = detect_all()
        needed = sorted({t for st in _STAGE_TOOLS.values() for t in st})
        require(self.tools, needed)
        for key in needed:
            st = self.tools[key]
            if not st.version_matches:
                self.console.warning(
                    f"{key} version {st.version} differs from the tested version "
                    f"{st.expected}; results may differ from the validated setup"
                )
        self.state = self._load_state()
        self._write_run_info()
        self._space_report()

        previous = self._install_signal_handlers()
        interrupted = False
        try:
            jobs = min(self.config.jobs, len(self.opts.samples))
            total = len(self.opts.samples)
            if jobs <= 1:
                for i, s in enumerate(self.opts.samples, 1):
                    self._process_sample(s, i, total)
            else:
                with ThreadPoolExecutor(max_workers=jobs, thread_name_prefix="sample") as pool:
                    futures = [
                        pool.submit(self._process_sample, s, i, total)
                        for i, s in enumerate(self.opts.samples, 1)
                    ]
                    errors = []
                    for f in futures:
                        try:
                            f.result()
                        except InterruptedRun as exc:
                            errors.append(exc)
                    if errors:
                        raise errors[0]
            if STOP_EVENT.is_set():
                raise InterruptedRun("interrupted")
            self._run_aggregates()
        except (InterruptedRun, KeyboardInterrupt):
            interrupted = True
            STOP_EVENT.set()
            terminate_all()
            self.state.mark_interrupted()
        finally:
            for sig, h in previous.items():
                signal.signal(sig, h)

        summary = self.build_summary(interrupted)
        try:
            written = reporting.write_reports(summary, self.layout)
            self.state.update(
                None,
                "report",
                status=Status.COMPLETED,
                finished=utc_now_iso(),
                outputs=outputs_record(written, self.layout.root),
            )
        except Exception as exc:  # report failure must not hide processing results
            self.console.detail(traceback.format_exc())
            self.console.error(
                f"report generation failed: {exc}",
                hint="Processing results are intact; see logs/pipeline.log.",
            )
            self.state.update(None, "report", status=Status.FAILED, message=str(exc))
            written = {}
        code = self._exit_code(summary, interrupted, report_ok=bool(written))
        self.state.data["runs"][-1].update({"finished": utc_now_iso(), "exit_code": code})
        self.state.save()
        self._final_message(summary, written, code)
        return code

    def _space_report(self) -> None:
        for label, path in (
            ("output directory", self.layout.root),
            ("temporary directory", self.tmp_root),
        ):
            free = free_bytes(path)
            self.console.info(f"Free space in {label} ({path}): {human_bytes(free)}")

    def _exit_code(self, summary: dict[str, Any], interrupted: bool, report_ok: bool) -> int:
        if interrupted:
            return ExitCode.INTERRUPTED
        c = summary["counts"]
        if c["completed"] == 0:
            return ExitCode.FAILED
        agg_failed = any(
            a.get("status") == Status.FAILED.value for a in summary["aggregates"].values()
        )
        if c["failed"] or c["skipped"] or agg_failed or not report_ok:
            return ExitCode.PARTIAL
        return ExitCode.OK

    def _final_message(self, summary: dict[str, Any], written: dict[str, Path], code: int) -> None:
        c = summary["counts"]
        rel = self.layout.rel
        lines = [
            "",
            self.console.style("RUN SUMMARY", "bold"),
            f"Samples: {c['attempted']} attempted, {c['completed']} completed, "
            f"{c['failed']} failed, {c['skipped']} skipped/incomplete",
        ]
        for s in summary["samples"]:
            extra = ""
            if s["failure"]:
                extra = f" — failed at {s['failure']['stage']}: {s['failure']['message']}"
            elif s["comparison"] and s["comparison"]["screening_flags"]:
                extra = f" — {len(s['comparison']['screening_flags'])} screening flag(s) for review"
            lines.append(f"  {s['sample']}: {s['status']}{extra}")
        if c["completed"]:
            lines.append(f"Cleaned FASTQ files: {self.layout.final_dir}")
        for key, label in (("html", "Summary report"), ("json", "Machine-readable summary")):
            if key in written:
                lines.append(f"{label}: {rel(written[key])}")
        for st, label in (("multiqc_raw", "Raw MultiQC"), ("multiqc_clean", "Cleaned MultiQC")):
            rec = summary["aggregates"].get(st, {})
            if rec.get("status") == "COMPLETED":
                html = (rec.get("outputs") or {}).get("html", {}).get("path")
                if html:
                    lines.append(f"{label}: {html}")
        if code == ExitCode.INTERRUPTED:
            lines.append(
                "The run was interrupted. Continue later with the same command plus --resume."
            )
        elif code != ExitCode.OK:
            lines.append(
                "Some steps failed. Fix the reported problem(s) and re-run with --resume; "
                "completed steps will be reused."
            )
        lines.append(
            "Next: open the summary report, then the MultiQC reports, and review any "
            "flagged samples before alignment."
        )
        self.console.summary("\n".join(lines))


def plan_reuse(pipeline: Pipeline, sample: Sample, stage: str) -> tuple[bool, str]:
    """Read-only reuse decision (used by dry-run)."""
    return pipeline._reusable(sample, stage)


__all__ = ["TOOLS", "Pipeline", "RunOptions", "detect_all", "outputs_intact"]
