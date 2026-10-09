"""Stage 3: adapter trimming and quality filtering with fastp."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from .config import PipelineConfig
from .context import StageContext, StageResult
from .errors import IntegrityError
from .metrics import MetricsError, read_fastp_json
from .runner import run_command
from .utils import remove_if_exists
from .validation import Sample

#: fastp uses at most 16 worker threads.
FASTP_MAX_THREADS = 16


def stage_parameters(config: PipelineConfig, paired: bool) -> dict[str, Any]:
    """The fastp settings that influence results (used for fingerprints and reports)."""
    p = config.section("fastp")
    p.pop("timeout_hours", None)
    if not paired:
        for key in ("adapter_sequence_r2", "detect_adapter_for_pe", "trim_front2", "trim_tail2"):
            p.pop(key, None)
    return p


def build_command(
    sample: Sample,
    inputs: list[Path],
    outputs: list[Path],
    json_path: Path,
    html_path: Path,
    config: PipelineConfig,
    threads: int,
) -> list[str]:
    """Translate the configuration into a fastp argument list."""
    p = config.section("fastp")
    cmd = ["fastp", "--in1", str(inputs[0]), "--out1", str(outputs[0])]
    if sample.paired:
        cmd += ["--in2", str(inputs[1]), "--out2", str(outputs[1])]
    cmd += [
        "--json",
        str(json_path),
        "--html",
        str(html_path),
        "--report_title",
        f"fastp report: {sample.name}",
        "--thread",
        str(max(1, min(threads, FASTP_MAX_THREADS))),
        "--compression",
        str(p["compression_level"]),
    ]
    if not p["adapter_trimming"]:
        cmd.append("--disable_adapter_trimming")
    else:
        if p["adapter_sequence"]:
            cmd += ["--adapter_sequence", p["adapter_sequence"]]
        if sample.paired:
            if p["adapter_sequence_r2"]:
                cmd += ["--adapter_sequence_r2", p["adapter_sequence_r2"]]
            if p["detect_adapter_for_pe"]:
                cmd.append("--detect_adapter_for_pe")
    if not p["quality_filtering"]:
        cmd.append("--disable_quality_filtering")
    else:
        cmd += [
            "--qualified_quality_phred",
            str(p["qualified_quality_phred"]),
            "--unqualified_percent_limit",
            str(p["unqualified_percent_limit"]),
            "--n_base_limit",
            str(p["n_base_limit"]),
        ]
        if p["average_qual"]:
            cmd += ["--average_qual", str(p["average_qual"])]
    if not p["length_filtering"]:
        cmd.append("--disable_length_filtering")
    else:
        cmd += ["--length_required", str(p["length_required"])]
        if p["length_limit"]:
            cmd += ["--length_limit", str(p["length_limit"])]
    if p["cut_front"] or p["cut_tail"]:
        if p["cut_front"]:
            cmd.append("--cut_front")
        if p["cut_tail"]:
            cmd.append("--cut_tail")
        cmd += [
            "--cut_window_size",
            str(p["cut_window_size"]),
            "--cut_mean_quality",
            str(p["cut_mean_quality"]),
        ]
    # Fixed trimming. For paired data R2 values are always passed explicitly,
    # because fastp otherwise silently copies the R1 values to R2.
    if p["trim_front1"]:
        cmd += ["--trim_front1", str(p["trim_front1"])]
    if p["trim_tail1"]:
        cmd += ["--trim_tail1", str(p["trim_tail1"])]
    if sample.paired and (
        p["trim_front1"] or p["trim_tail1"] or p["trim_front2"] or p["trim_tail2"]
    ):
        cmd += ["--trim_front2", str(p["trim_front2"]), "--trim_tail2", str(p["trim_tail2"])]
    if p["trim_poly_g"] == "on":
        cmd.append("--trim_poly_g")
    elif p["trim_poly_g"] == "off":
        cmd.append("--disable_trim_poly_g")
    if p["trim_poly_x"]:
        cmd.append("--trim_poly_x")
    if p["low_complexity_filter"]:
        cmd += ["--low_complexity_filter", "--complexity_threshold", str(p["complexity_threshold"])]
    if p["overrepresentation_analysis"]:
        cmd.append("--overrepresentation_analysis")
    return cmd


def describe_parameters(config: PipelineConfig, paired: bool) -> list[str]:
    """Plain-language description of every enabled fastp operation."""
    p = config.section("fastp")
    lines: list[str] = []
    if p["adapter_trimming"]:
        if p["adapter_sequence"]:
            lines.append(
                f"Adapter trimming with the given adapter {p['adapter_sequence']}"
                + (
                    f" / R2 {p['adapter_sequence_r2']}"
                    if paired and p["adapter_sequence_r2"]
                    else ""
                )
            )
        elif paired:
            lines.append(
                "Adapter trimming: R1/R2 overlap analysis"
                + (" plus adapter sequence auto-detection" if p["detect_adapter_for_pe"] else "")
                + " (removes adapter read-through from short inserts)"
            )
        else:
            lines.append("Adapter trimming with automatic adapter detection")
    else:
        lines.append("Adapter trimming: DISABLED")
    if p["quality_filtering"]:
        lines.append(
            f"Quality filter: discard a read if more than {p['unqualified_percent_limit']}% of its "
            f"bases have Phred quality < {p['qualified_quality_phred']}, or it has more than "
            f"{p['n_base_limit']} N bases"
            + (f", or its mean quality is < {p['average_qual']}" if p["average_qual"] else "")
        )
    else:
        lines.append("Quality filtering: DISABLED")
    if p["length_filtering"]:
        lines.append(
            f"Length filter: discard reads shorter than {p['length_required']} bases after trimming"
            + (f" or longer than {p['length_limit']}" if p["length_limit"] else "")
        )
    else:
        lines.append("Length filtering: DISABLED")
    if p["cut_front"] or p["cut_tail"]:
        ends = " and ".join(e for e, on in (("5'", p["cut_front"]), ("3'", p["cut_tail"])) if on)
        lines.append(
            f"Sliding-window trimming at the {ends} end(s): window {p['cut_window_size']}, "
            f"mean quality >= {p['cut_mean_quality']}"
        )
    trims = [(k, p[k]) for k in ("trim_front1", "trim_tail1", "trim_front2", "trim_tail2") if p[k]]
    if trims:
        lines.append("Fixed trimming: " + ", ".join(f"{k}={v}" for k, v in trims))
    lines.append(
        {
            "auto": "polyG trimming: automatic (fastp enables it only when the read header identifies "
            "NextSeq/NovaSeq data, which SRA-converted headers usually do not)",
            "on": "polyG tail trimming: forced on",
            "off": "polyG tail trimming: disabled",
        }[p["trim_poly_g"]]
    )
    if p["trim_poly_x"]:
        lines.append("polyX tail trimming: on")
    if p["low_complexity_filter"]:
        lines.append(f"Low-complexity filter: on (threshold {p['complexity_threshold']}%)")
    lines.append("Deduplication: OFF (duplicates are expected in RNA-seq and are kept)")
    return lines


def run(sample: Sample, raw: list[Path], ctx: StageContext) -> StageResult:
    result = StageResult()
    layout = ctx.layout
    outdir = layout.fastp_sample_dir(sample)
    outdir.mkdir(parents=True, exist_ok=True)
    finals = layout.clean_fastq(sample)
    partials = [f.with_name(f.name.replace(".fastq.gz", ".partial.fastq.gz")) for f in finals]
    json_path, html_path = layout.fastp_json(sample), layout.fastp_html(sample)
    for p in (*finals, *partials, json_path, html_path):
        remove_if_exists(p)
    cmd = build_command(sample, raw, partials, json_path, html_path, ctx.config, ctx.threads)
    run_command(
        cmd,
        log_path=layout.command_log(sample, "fastp"),
        timeout=float(ctx.config["fastp.timeout_hours"]) * 3600,
        description=f"fastp {sample.name}",
    )
    for p in (*partials, json_path, html_path):
        if not p.is_file() or p.stat().st_size == 0:
            raise IntegrityError(
                f"fastp finished but {p.name} is missing or empty",
                hint=f"See {layout.command_log(sample, 'fastp')}",
            )
    try:
        m = read_fastp_json(json_path, sample.paired)
    except MetricsError as exc:
        raise IntegrityError(f"fastp JSON report cannot be used: {exc}") from exc
    for part, final in zip(partials, finals, strict=True):
        os.replace(part, final)
    labels = ["R1", "R2"] if sample.paired else ["reads"]
    for label, f in zip(labels, finals, strict=True):
        result.outputs[label] = f
    result.outputs["json"] = json_path
    result.outputs["html"] = html_path
    result.metrics = m.to_dict()
    result.metrics["command_line"] = cmd
    if m.reads_after == 0:
        result.warnings.append("fastp removed ALL reads; check the parameters and the raw QC.")
    return result
