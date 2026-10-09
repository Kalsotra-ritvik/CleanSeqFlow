from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from transcriptomics_preprocessor.comparison import build_comparison
from transcriptomics_preprocessor.config import build_config
from transcriptomics_preprocessor.layout import OutputLayout
from transcriptomics_preprocessor.metrics import read_fastp_json
from transcriptomics_preprocessor.reporting import render_html, write_reports

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "tool_outputs"
CFG = build_config()


def _val(records: int) -> dict[str, Any]:
    return {"records": records}


def _rows(cmp: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {r["metric"]: r for r in cmp["rows"]}


def test_paired_reads_vs_pairs_and_consistency() -> None:
    fp = read_fastp_json(FIX / "fastp_1.3.6_paired.json", True).to_dict()
    cmp = build_comparison(
        "s",
        True,
        {"R1": _val(3000), "R2": _val(3000)},
        {"R1": _val(2835), "R2": _val(2835)},
        fp,
        {"R1": {"max_adapter_percent": 19.9}, "R2": {"max_adapter_percent": 20.0}},
        {"R1": {"max_adapter_percent": 0.0}, "R2": {"max_adapter_percent": 0.0}},
        CFG,
    ).to_dict()
    rows = _rows(cmp)
    assert rows["Reads"]["before"] == 6000 and rows["Reads"]["after"] == 5670
    assert rows["Read pairs"]["before"] == 3000 and rows["Read pairs"]["after"] == 2835
    assert rows["Reads retained"]["after"] == 94.5
    assert cmp["consistency_issues"] == []
    assert any("decreased from 19.90% to 0.00%" in o for o in cmp["observations"])
    # 2,835 pairs < default screening threshold of 10 million -> flagged for review, not "failed"
    assert cmp["screening_status"] == "REVIEW"


def test_count_inconsistency_is_reported() -> None:
    fp = read_fastp_json(FIX / "fastp_1.3.6_paired.json", True).to_dict()
    cmp = build_comparison(
        "s",
        True,
        {"R1": _val(2999), "R2": _val(2999)},
        {"R1": _val(2835), "R2": _val(2835)},
        fp,
        None,
        None,
        CFG,
    ).to_dict()
    assert cmp["screening_status"] == "INCONSISTENT"
    assert any("raw files contain 5,998" in i for i in cmp["consistency_issues"])


def test_missing_metrics_remain_unavailable() -> None:
    cmp = build_comparison("s", False, None, None, None, None, None, CFG).to_dict()
    for r in cmp["rows"]:
        assert r["before"] in (None,) and r["after"] is None
        assert r["change"] == "not available"
    assert any("not available" in o for o in cmp["observations"])
    assert cmp["screening_flags"] == []  # no flag can be raised without data


def test_no_invented_improvement_and_negligible_change() -> None:
    fp = read_fastp_json(FIX / "fastp_1.3.6_single.json", False).to_dict()
    cmp = build_comparison(
        "s",
        False,
        {"reads": _val(1500)},
        {"reads": _val(1434)},
        fp,
        {"reads": {"max_adapter_percent": 19.4}},
        {"reads": {"max_adapter_percent": 19.386}},
        CFG,
    ).to_dict()
    obs = " ".join(cmp["observations"])
    assert "essentially unchanged" in obs and "decreased" not in obs
    assert "automatic adapter detection probably failed" in obs
    assert any("Residual adapter content" in f for f in cmp["screening_flags"])


def test_increase_is_reported_honestly() -> None:
    cmp = build_comparison(
        "s",
        False,
        None,
        None,
        None,
        {"reads": {"max_adapter_percent": 1.0}},
        {"reads": {"max_adapter_percent": 3.0}},
        CFG,
    ).to_dict()
    assert any("INCREASED" in o for o in cmp["observations"])


def test_screening_thresholds_are_configurable() -> None:
    fp = read_fastp_json(FIX / "fastp_1.3.6_paired.json", True).to_dict()
    lenient = build_config(overrides={"screening.min_clean_reads": 100})
    cmp = build_comparison(
        "s",
        True,
        {"R1": _val(3000), "R2": _val(3000)},
        {"R1": _val(2835), "R2": _val(2835)},
        fp,
        None,
        None,
        lenient,
    ).to_dict()
    assert cmp["screening_flags"] == [] and cmp["screening_status"] == "NO_FLAGS"
    strict = build_config(
        overrides={"screening.min_retained_percent": 99, "screening.min_clean_q30_rate": 0.9}
    )
    cmp = build_comparison(
        "s",
        True,
        {"R1": _val(3000), "R2": _val(3000)},
        {"R1": _val(2835), "R2": _val(2835)},
        fp,
        None,
        None,
        strict,
    ).to_dict()
    text = " ".join(cmp["screening_flags"])
    assert "94.50% of reads were retained" in text and "Q30 rate" in text


def _summary(tmp_path: Path, sample_extra: dict[str, Any] | None = None) -> dict[str, Any]:
    stages = {
        st: {"status": "NOT_STARTED"}
        for st in (
            "download",
            "convert",
            "validate_raw",
            "fastqc_raw",
            "fastp",
            "validate_clean",
            "fastqc_clean",
            "publish",
        )
    }
    stages["download"] = {
        "status": "FAILED",
        "message": "prefetch failed <script>",
        "hint": "check network",
    }
    s = {
        "sample": "s<1>",
        "accession": "SRR1234567",
        "layout": "paired",
        "metadata": {"condition": "a&b"},
        "status": "FAILED",
        "failure": {"stage": "download", "message": "prefetch failed <script>", "hint": "x"},
        "stages": stages,
        "comparison": None,
        "final_files": None,
        "removed_by_cleanup": False,
    }
    s.update(sample_extra or {})
    return {
        "tool": "t",
        "version": "0.1.0",
        "project": "p",
        "started": "a",
        "finished": "b",
        "interrupted": False,
        "mode": "new",
        "command_line": ["x"],
        "input": "i",
        "outdir": str(tmp_path),
        "tmpdir": "t",
        "environment": {"python": "3", "platform": "linux", "environment_prefix": "/e"},
        "tools": {"fastp": {"version": None, "expected": "1.3.6", "path": None}},
        "config": CFG.as_nested(),
        "counts": {"attempted": 1, "completed": 0, "failed": 1, "skipped": 0},
        "samples": [s],
        "aggregates": {
            "multiqc_raw": {"status": "SKIPPED"},
            "multiqc_clean": {"status": "SKIPPED"},
        },
        "warnings": [],
        "console_warnings": [],
    }


def test_report_with_failures_and_missing_metrics_is_escaped(tmp_path: Path) -> None:
    html = render_html(_summary(tmp_path))
    assert "<script>" not in html and "&lt;script&gt;" in html
    assert "No sample completed processing" in html
    assert "not available" in html  # missing tool version shown as such


def test_write_reports_files(tmp_path: Path) -> None:
    layout = OutputLayout(tmp_path)
    written = write_reports(_summary(tmp_path), layout)
    for key in ("json", "status_csv", "comparison_csv", "manifest_csv", "html", "inventory_csv"):
        assert written[key].is_file() and written[key].stat().st_size > 0
    assert json.loads(written["json"].read_text())["counts"]["failed"] == 1
    rows = list(csv.reader(written["status_csv"].open()))
    assert rows[1][3] == "FAILED" and rows[1][-2] == "download"
    assert "final_manifest" not in written  # no completed sample -> no final manifest


def test_change_wording_uses_percentage_points() -> None:
    from transcriptomics_preprocessor.comparison import ComparisonRow

    assert (
        ComparisonRow("q", "d", "s", 73.13, 79.84, "%").change
        == "increased by 6.71 percentage points"
    )
    assert ComparisonRow("r", "d", "s", 3000, 2834).change == "decreased by 166"
    assert ComparisonRow("l", "d", "s", 100, 92, " bp").change == "decreased by 8 bp"
    assert ComparisonRow("x", "d", "s", 1.0, 1.0).change == "unchanged"
    assert ComparisonRow("x", "d", "s", None, 1.0).change == "not available"
