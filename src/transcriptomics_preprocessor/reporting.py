"""Final reports: HTML summary, JSON summary, CSV tables, final manifest and checksums.

The HTML report is self-contained (no external scripts, fonts or images) and
separates five kinds of information:

1. observed measurements (counts from the FASTQ validator),
2. tool-generated statistics (fastp, FastQC values),
3. software quality flags (FastQC PASS/WARN/FAIL),
4. automated screening recommendations (configurable heuristics),
5. items that need expert review.
"""

from __future__ import annotations

import csv
import html
import io
import os
from pathlib import Path
from typing import Any

from .explain import FLAG_DISCLAIMER, LIMITATIONS, QC_MODULES, STAGE_LABELS
from .layout import OutputLayout
from .state import SAMPLE_STAGES
from .utils import atomic_write_json, atomic_write_text, human_bytes

NA = "not available"


def _fmt(v: Any, unit: str = "") -> str:
    if v is None:
        return NA
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, int):
        return f"{v:,}{unit}"
    if isinstance(v, float):
        return f"{v:,.2f}{unit}"
    return f"{v}{unit}"


def _csv(rows: list[list[Any]]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    for r in rows:
        w.writerow(["" if v is None else v for v in r])
    return buf.getvalue()


# ----------------------------------------------------------------- tables
def status_rows(summary: dict[str, Any]) -> list[list[Any]]:
    rows: list[list[Any]] = [
        [
            "sample",
            "accession",
            "layout",
            "overall_status",
            *SAMPLE_STAGES,
            "failed_stage",
            "failure_message",
        ]
    ]
    for s in summary["samples"]:
        f = s["failure"] or {}
        rows.append(
            [
                s["sample"],
                s["accession"],
                s["layout"],
                s["status"],
                *[s["stages"][st].get("status", "NOT_STARTED") for st in SAMPLE_STAGES],
                f.get("stage"),
                f.get("message"),
            ]
        )
    return rows


def comparison_rows(summary: dict[str, Any]) -> list[list[Any]]:
    rows: list[list[Any]] = [
        ["sample", "metric", "before", "after", "change", "unit", "definition", "source"]
    ]
    for s in summary["samples"]:
        for r in (s["comparison"] or {}).get("rows", []):
            rows.append(
                [
                    s["sample"],
                    r["metric"],
                    r["before"],
                    r["after"],
                    r["change"],
                    r["unit"].strip(),
                    r["definition"],
                    r["source"],
                ]
            )
    return rows


def manifest_rows(summary: dict[str, Any]) -> list[list[Any]]:
    meta_keys: list[str] = []
    for s in summary["samples"]:
        for k in s["metadata"]:
            if k not in meta_keys:
                meta_keys.append(k)
    rows: list[list[Any]] = [
        [
            "sample",
            "accession",
            "layout",
            *meta_keys,
            "status",
            "screening_status",
            "r1_or_single_fastq",
            "r2_fastq",
            "reads",
            "read_pairs",
        ]
    ]
    for s in summary["samples"]:
        ff = s["final_files"] or {}
        first = ff.get("R1") or ff.get("reads") or {}
        second = ff.get("R2") or {}
        reads = None
        if first.get("records") is not None:
            reads = (
                first["records"] + (second.get("records") or 0)
                if s["layout"] == "paired"
                else first["records"]
            )
        rows.append(
            [
                s["sample"],
                s["accession"],
                s["layout"],
                *[s["metadata"].get(k, "") for k in meta_keys],
                s["status"],
                (s["comparison"] or {}).get("screening_status"),
                first.get("file"),
                second.get("file"),
                reads,
                first.get("records") if s["layout"] == "paired" else None,
            ]
        )
    return rows


def inventory_rows(layout: OutputLayout) -> list[list[Any]]:
    """Every file under the output directory except temporary files, with size."""
    rows: list[list[Any]] = [["path", "category", "bytes", "size"]]
    root = layout.root
    categories = {
        "01_sra": "SRA archive (raw download)",
        "02_raw_fastq": "raw FASTQ",
        "03_raw_qc": "raw QC report",
        "04_fastp": "fastp output (intermediate)",
        "05_clean_qc": "post-cleaning QC report",
        "06_final_fastq": "FINAL cleaned FASTQ / manifest",
        "07_reports": "summary report",
        "logs": "log",
        "pipeline_info": "run information",
    }
    for dirpath, dirnames, filenames in os.walk(root):
        rel_dir = Path(dirpath).relative_to(root)
        top = rel_dir.parts[0] if rel_dir.parts else ""
        if top == "tmp" or any(part.startswith(".work_") for part in rel_dir.parts):
            dirnames[:] = []
            continue
        dirnames.sort()
        for name in sorted(filenames):
            p = Path(dirpath) / name
            if p.is_symlink() or name.endswith((".tmp", ".partial")):
                continue
            size = p.stat().st_size
            rows.append(
                [str(p.relative_to(root)), categories.get(top, "other"), size, human_bytes(size)]
            )
    return rows


# ------------------------------------------------------------------- html
_CSS = """
:root{--fg:#1d2329;--muted:#5b6670;--bg:#fff;--panel:#f5f7f9;--line:#d9dee3;--ok:#1a7f37;--warn:#9a6700;
--bad:#cf222e;--info:#0969da}
@media (prefers-color-scheme:dark){:root{--fg:#e6edf3;--muted:#9da7b1;--bg:#0d1117;--panel:#161b22;
--line:#30363d;--ok:#3fb950;--warn:#d29922;--bad:#f85149;--info:#58a6ff}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.5 system-ui,-apple-system,Segoe UI,Roboto,sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 64px}
h1{font-size:1.6rem;margin:.2em 0}h2{font-size:1.25rem;margin-top:2em;border-bottom:1px solid var(--line);
padding-bottom:.3em}h3{font-size:1.05rem;margin-top:1.4em}
.muted{color:var(--muted)}.panel{background:var(--panel);border:1px solid var(--line);border-radius:8px;
padding:12px 16px;margin:12px 0}
.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;margin:8px 0;font-size:.92rem}
th,td{border:1px solid var(--line);padding:5px 8px;text-align:left;vertical-align:top}
th{background:var(--panel)}td.num{text-align:right;font-variant-numeric:tabular-nums}
.badge{display:inline-block;padding:1px 8px;border-radius:10px;font-size:.8rem;font-weight:600;
border:1px solid currentColor}
.COMPLETED,.PASS,.NO_FLAGS{color:var(--ok)}.FAILED,.FAIL,.INCONSISTENT{color:var(--bad)}
.WARN,.REVIEW,.SKIPPED,.INCOMPLETE{color:var(--warn)}.NOT_STARTED,.RUNNING{color:var(--muted)}
.tag{font-size:.75rem;text-transform:uppercase;letter-spacing:.04em;color:var(--info);font-weight:700}
code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.88em;word-break:break-all}
ul{padding-left:1.3em}details{margin:8px 0}summary{cursor:pointer;font-weight:600}
"""


def _e(x: Any) -> str:
    return html.escape(str(x), quote=True)


def _badge(status: str | None) -> str:
    s = status or "NOT_STARTED"
    return f'<span class="badge {_e(s)}">{_e(s)}</span>'


def _table(header: list[str], rows: list[list[str]], num_cols: tuple[int, ...] = ()) -> str:
    h = "".join(f"<th>{_e(c)}</th>" for c in header)
    body = []
    for r in rows:
        cells = "".join(
            f'<td class="num">{c}</td>' if i in num_cols else f"<td>{c}</td>"
            for i, c in enumerate(r)
        )
        body.append(f"<tr>{cells}</tr>")
    return f'<div class="scroll"><table><thead><tr>{h}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def _link(path: str | None, label: str | None = None) -> str:
    if not path:
        return NA
    # Links are relative to 07_reports/.
    return f'<a href="../{_e(path)}">{_e(label or path)}</a>'


def render_html(summary: dict[str, Any]) -> str:
    c = summary["counts"]
    parts: list[str] = []
    add = parts.append
    title = f"RNA-seq preprocessing report — {summary['project']}"
    add(
        f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
        f"content='width=device-width,initial-scale=1'><title>{_e(title)}</title><style>{_CSS}</style>"
        "</head><body><main>"
    )
    add(f"<h1>{_e(title)}</h1>")
    add(
        f"<p class='muted'>transcriptomics-preprocess {_e(summary['version'])} · started "
        f"{_e(summary['started'])} · finished {_e(summary['finished'])} (UTC) · mode {_e(summary['mode'])}</p>"
    )

    # At a glance
    if summary["interrupted"]:
        verdict = "The run was INTERRUPTED; results below are partial. Re-run with --resume."
    elif c["completed"] == c["attempted"] and c["attempted"]:
        verdict = "All samples completed every processing step."
    elif c["completed"]:
        verdict = "Some samples did not complete; see failures below."
    else:
        verdict = "No sample completed processing."
    flagged = [
        s["sample"]
        for s in summary["samples"]
        if s["comparison"] and s["comparison"]["screening_status"] != "NO_FLAGS"
    ]
    add("<h2>At a glance</h2><div class='panel'>")
    add(f"<p><strong>{_e(verdict)}</strong></p>")
    add(
        f"<p>Samples attempted: {c['attempted']} · completed: {c['completed']} · failed: {c['failed']} · "
        f"skipped/incomplete: {c['skipped']}</p>"
    )
    add("<p>What to inspect next:</p><ol>")
    if c["failed"]:
        add(
            "<li>Read the failure messages in <em>Processing status</em> and fix the cause, then re-run "
            "with <code>--resume</code>.</li>"
        )
    if flagged:
        add(f"<li>Review samples with screening flags: {_e(', '.join(flagged))}.</li>")
    add(
        "<li>Open the raw and cleaned MultiQC reports (links below) and compare the plots.</li>"
        "<li>Check that read counts per sample are what you expect from the SRA records.</li>"
        "<li>Only then use the files in <code>06_final_fastq/</code> for alignment.</li></ol></div>"
    )

    # Status
    add("<h2>Processing status</h2>")
    rows = []
    for s in summary["samples"]:
        cond = ", ".join(f"{k}={v}" for k, v in s["metadata"].items())
        scr = (s["comparison"] or {}).get("screening_status")
        rows.append(
            [
                _e(s["sample"]),
                _e(s["accession"]),
                _e(s["layout"]),
                _e(cond) or "—",
                _badge(s["status"]),
                _badge(scr) if scr else NA,
            ]
        )
    add(_table(["Sample", "Accession", "Layout", "Metadata", "Status", "Screening"], rows))
    add("<details><summary>Status of every step</summary>")
    rows = []
    for s in summary["samples"]:
        rows.append(
            [_e(s["sample"])] + [_badge(s["stages"][st].get("status")) for st in SAMPLE_STAGES]
        )
    add(_table(["Sample"] + [STAGE_LABELS[st] for st in SAMPLE_STAGES], rows))
    add("</details>")
    for s in summary["samples"]:
        if s["failure"]:
            f = s["failure"]
            add(
                f"<div class='panel'><span class='tag'>Failure</span><p><strong>{_e(s['sample'])}</strong> "
                f"failed at <em>{_e(STAGE_LABELS.get(f['stage'], f['stage']))}</em>: {_e(f['message'])}</p>"
                + (f"<p>What to do: {_e(f['hint'])}</p>" if f.get("hint") else "")
                + "</div>"
            )

    # Outputs
    add("<h2>Outputs</h2>")
    agg = summary["aggregates"]

    def agg_link(st: str) -> str:
        rec = agg.get(st) or {}
        if rec.get("status") != "COMPLETED":
            return f"{_badge(rec.get('status'))} {_e(rec.get('message') or '')}"
        return _link(((rec.get("outputs") or {}).get("html") or {}).get("path"))

    add(
        _table(
            ["Report", "Location"],
            [
                ["Raw reads MultiQC (before fastp)", agg_link("multiqc_raw")],
                [
                    "Cleaned reads MultiQC (after fastp, with fastp statistics)",
                    agg_link("multiqc_clean"),
                ],
                ["Machine-readable summary", _link("07_reports/summary.json")],
                ["Before/after comparison table", _link("07_reports/qc_comparison.csv")],
                ["Final FASTQ manifest", _link("06_final_fastq/MANIFEST.csv")],
                ["File inventory", _link("07_reports/file_inventory.csv")],
            ],
        )
    )
    rows = []
    for s in summary["samples"]:
        for label, info in (s["final_files"] or {}).items():
            rows.append(
                [
                    _e(s["sample"]),
                    _e(label),
                    _e(info.get("file")),
                    _e(_fmt(info.get("records"))),
                    _e(human_bytes(info.get("bytes"))),
                    f"<code>{_e(info.get('sha256') or NA)}</code>",
                ]
            )
    if rows:
        add("<h3>Final cleaned FASTQ files</h3>")
        add(_table(["Sample", "File", "Path", "Reads in file", "Size", "SHA-256"], rows, (3, 4)))

    # Per-sample comparisons
    add("<h2>Before-versus-after comparison</h2>")
    add(
        "<p class='muted'>Each row compares the same statistic from the same source before and after "
        "fastp. <em>Reads</em> are individual reads; for paired-end data one read pair = two reads. "
        "Values that could not be obtained are shown as “not available”, never estimated.</p>"
    )
    for s in summary["samples"]:
        cmpd = s["comparison"]
        if not cmpd:
            continue
        add(
            f"<h3>{_e(s['sample'])} <span class='muted'>({_e(s['accession'])}, {_e(s['layout'])}-end)</span> "
            f"{_badge(cmpd['screening_status'])}</h3>"
        )
        add("<span class='tag'>Observed measurements and tool-generated statistics</span>")
        rows = [
            [
                _e(r["metric"]),
                _e(_fmt(r["before"])),
                _e(_fmt(r["after"])),
                _e(r["change"]),
                _e(r["source"]),
                f"<span class='muted'>{_e(r['definition'])}</span>",
            ]
            for r in cmpd["rows"]
        ]
        add(_table(["Metric", "Before", "After", "Change", "Source", "Definition"], rows, (1, 2)))
        removed = cmpd.get("reads_removed_by_reason") or {}
        if removed:
            add(
                "<p>Reads removed by fastp, by reason (individual reads): "
                + ", ".join(f"{_e(k)}: {_e(_fmt(v))}" for k, v in removed.items())
                + "</p>"
            )
        if cmpd["observations"]:
            add("<ul>" + "".join(f"<li>{_e(o)}</li>" for o in cmpd["observations"]) + "</ul>")
        if cmpd["consistency_issues"]:
            add(
                "<div class='panel'><span class='tag'>Consistency check failed — needs review</span><ul>"
                + "".join(f"<li class='FAIL'>{_e(i)}</li>" for i in cmpd["consistency_issues"])
                + "</ul></div>"
            )
        add("<span class='tag'>Automated screening (heuristics; configurable)</span>")
        if cmpd["screening_flags"]:
            add(
                "<ul>"
                + "".join(f"<li class='WARN'>{_e(f)}</li>" for f in cmpd["screening_flags"])
                + "</ul>"
            )
        else:
            add(
                "<p>No screening rule was triggered. This does not prove the sample is suitable; "
                "it only means none of the configured checks fired.</p>"
            )

    # FastQC flags
    add("<h2>FastQC module flags</h2>")
    add(
        f"<div class='panel'><span class='tag'>Software quality flags</span><p>{_e(FLAG_DISCLAIMER)}</p></div>"
    )
    mods = list(QC_MODULES)
    rows = []
    for s in summary["samples"]:
        for when, st in (("raw", "fastqc_raw"), ("cleaned", "fastqc_clean")):
            m = s["stages"][st].get("metrics") or {}
            for mate, vals in m.items():
                flags = vals.get("module_status") or {}
                rows.append(
                    [_e(s["sample"]), _e(mate), _e(when)]
                    + [_badge(flags[k]) if k in flags else "—" for k in mods]
                )
    if rows:
        add(_table(["Sample", "File", "Reads"] + [QC_MODULES[k][0] for k in mods], rows))
    add(
        "<details><summary>What each FastQC module assesses</summary><ul>"
        + "".join(f"<li><strong>{_e(t)}</strong>: {_e(d)}</li>" for t, d in QC_MODULES.values())
        + "</ul></details>"
    )

    # Warnings
    add("<h2>Warnings</h2>")
    warns = summary["warnings"] + [
        {"sample": w.get("sample"), "stage": "", "message": w["message"]}
        for w in summary["console_warnings"]
        if w["message"] not in {x["message"] for x in summary["warnings"]}
    ]
    if warns:
        add(
            "<ul>"
            + "".join(
                f"<li>{_e(w['sample'] or 'run')}{(' / ' + _e(w['stage'])) if w['stage'] else ''}: {_e(w['message'])}</li>"
                for w in warns
            )
            + "</ul>"
        )
    else:
        add("<p>No warnings were raised.</p>")

    # Parameters & versions
    add("<h2>Software and parameters</h2>")
    rows = []
    for k, t in summary["tools"].items():
        rows.append(
            [_e(k), _e(t.get("version") or NA), _e(t.get("expected")), _e(t.get("path") or NA)]
        )
    add(_table(["Tool", "Version used", "Tested version", "Path"], rows))
    env = summary["environment"]
    add(
        f"<p class='muted'>Python {_e(env['python'])} · {_e(env['platform'])} · environment "
        f"{_e(env['environment_prefix'])}</p>"
    )
    fp = summary["config"].get("fastp", {})
    add(_table(["fastp parameter", "Value"], [[_e(k), _e(v)] for k, v in fp.items()]))
    add(
        "<details><summary>Full configuration</summary>"
        + _table(
            ["Section", "Key", "Value"],
            [
                [_e(sec), _e(k), _e(v)]
                for sec, vals in summary["config"].items()
                for k, v in (vals.items() if isinstance(vals, dict) else [("", vals)])
            ],
        )
        + "</details>"
    )
    add(f"<p>Command: <code>{_e(' '.join(summary['command_line']))}</code></p>")
    add(f"<p>Input: {_e(summary['input'])}</p>")

    add("<h2>Interpretation and limitations</h2>")
    add(
        "<div class='panel'><span class='tag'>Conclusions requiring expert review</span><p>Whether a "
        "sample is suitable for a particular analysis is a scientific judgement. The values above "
        "describe technical read properties; they do not by themselves establish sample identity, "
        "biological validity, or adequacy for differential expression.</p></div><ul>"
        + "".join(f"<li>{_e(x)}</li>" for x in LIMITATIONS)
        + "</ul>"
    )
    add("</main></body></html>")
    return "\n".join(parts)


def _final_manifest(summary: dict[str, Any], layout: OutputLayout) -> dict[str, Path]:
    written: dict[str, Path] = {}
    rows: list[list[Any]] = [
        ["sample", "accession", "layout", "file_role", "path", "reads", "bytes", "sha256"]
    ]
    sums: list[str] = []
    for s in summary["samples"]:
        for label, info in (s["final_files"] or {}).items():
            name = Path(info["file"]).name
            rows.append(
                [
                    s["sample"],
                    s["accession"],
                    s["layout"],
                    label,
                    name,
                    info.get("records"),
                    info.get("bytes"),
                    info.get("sha256"),
                ]
            )
            if info.get("sha256"):
                sums.append(f"{info['sha256']}  {name}\n")
    if len(rows) > 1:
        layout.final_dir.mkdir(parents=True, exist_ok=True)
        p = layout.final_dir / "MANIFEST.csv"
        atomic_write_text(p, _csv(rows))
        written["final_manifest"] = p
        if sums:
            q = layout.final_dir / "SHA256SUMS"
            atomic_write_text(q, "".join(sums))
            written["checksums"] = q
    return written


def write_reports(summary: dict[str, Any], layout: OutputLayout) -> dict[str, Path]:
    """Write all report files atomically and return their paths."""
    d = layout.reports_dir
    d.mkdir(parents=True, exist_ok=True)
    written = _final_manifest(summary, layout)
    targets = {
        "json": d / "summary.json",
        "status_csv": d / "processing_status.csv",
        "comparison_csv": d / "qc_comparison.csv",
        "manifest_csv": d / "sample_manifest.csv",
        "html": d / "summary_report.html",
        "inventory_csv": d / "file_inventory.csv",
    }
    atomic_write_json(targets["json"], summary)
    atomic_write_text(targets["status_csv"], _csv(status_rows(summary)))
    atomic_write_text(targets["comparison_csv"], _csv(comparison_rows(summary)))
    atomic_write_text(targets["manifest_csv"], _csv(manifest_rows(summary)))
    atomic_write_text(targets["html"], render_html(summary))
    # Inventory last so it lists the other reports too (it cannot list itself).
    atomic_write_text(targets["inventory_csv"], _csv(inventory_rows(layout)))
    written.update(targets)
    return written
