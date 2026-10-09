"""Parsers tested against REAL FastQC 0.12.1 / fastp 1.3.6 outputs (tests/fixtures/tool_outputs).

Expected values were read from the fixture files independently (awk / json
module), not produced by the parsers under test.
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pytest

from transcriptomics_preprocessor.metrics import (
    MetricsError,
    parse_fastp_json,
    parse_fastqc_data,
    read_fastp_json,
    read_fastqc_zip,
)

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "tool_outputs"


def test_fastqc_data_real_output() -> None:
    m = parse_fastqc_data((FIX / "fastqc_0.12.1_raw_R1_fastqc_data.txt").read_text())
    assert m.total_sequences == 3000
    assert m.gc_percent == 49
    assert m.sequence_length == "100"
    assert m.poor_quality_sequences == 0
    assert m.max_adapter_percent == pytest.approx(19.9)
    assert m.adapter_type_at_max == "Illumina Universal Adapter"
    assert m.min_position_mean_quality == pytest.approx(24.069)
    assert m.deduplicated_percent == pytest.approx(88.86666666666667)
    assert m.overrepresented_count == 2
    assert m.module_status["adapter_content"] == "FAIL"
    assert m.module_status["basic_statistics"] == "PASS"


def test_fastqc_zip_uses_summary_flags(tmp_path: Path) -> None:
    z = tmp_path / "x_fastqc.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr(
            "x_fastqc/fastqc_data.txt", (FIX / "fastqc_0.12.1_raw_R1_fastqc_data.txt").read_text()
        )
        zf.writestr("x_fastqc/summary.txt", (FIX / "fastqc_0.12.1_raw_R1_summary.txt").read_text())
    m = read_fastqc_zip(z)
    assert m.module_status["per_sequence_gc"] == "FAIL"
    assert m.module_status["per_base_content"] == "WARN"
    assert len(m.module_status) == 10


def test_fastqc_bad_inputs(tmp_path: Path) -> None:
    with pytest.raises(MetricsError):
        parse_fastqc_data("not fastqc")
    bad = tmp_path / "bad.zip"
    bad.write_bytes(b"PK\x03\x04garbage")
    with pytest.raises(MetricsError):
        read_fastqc_zip(bad)
    empty = tmp_path / "e.zip"
    with zipfile.ZipFile(empty, "w") as zf:
        zf.writestr("x/other.txt", "x")
    with pytest.raises(MetricsError, match=r"fastqc_data\.txt"):
        read_fastqc_zip(empty)


def test_fastqc_missing_modules_stay_none() -> None:
    m = parse_fastqc_data(
        "##FastQC\t0.12.1\n>>Basic Statistics\tpass\n#Measure\tValue\nTotal Sequences\t5\n>>END_MODULE\n"
    )
    assert m.total_sequences == 5
    assert m.max_adapter_percent is None and m.deduplicated_percent is None
    assert m.min_position_mean_quality is None


def test_fastp_paired_real_output() -> None:
    m = read_fastp_json(FIX / "fastp_1.3.6_paired.json", paired=True)
    assert (m.reads_before, m.reads_after) == (6000, 5670)
    assert (m.pairs_before, m.pairs_after) == (3000, 2835)  # reads vs read pairs
    assert m.low_quality_reads == 330 and m.passed_filter_reads == 5670
    assert m.adapter_trimmed_reads == 1204
    assert m.q30_rate_after == pytest.approx(0.798549)
    assert m.insert_size_peak == 57
    assert m.retained_fraction == pytest.approx(5670 / 6000)
    assert m.fastp_version == "1.3.6"


def test_fastp_single_real_output_without_adapter_section() -> None:
    m = read_fastp_json(FIX / "fastp_1.3.6_single.json", paired=False)
    assert (m.reads_before, m.reads_after) == (1500, 1434)
    assert m.pairs_before is None  # single-end: no pairs
    assert m.adapter_trimmed_reads is None  # absent in the report -> not available, not 0


def test_fastp_layout_mismatch_detected() -> None:
    data = json.loads((FIX / "fastp_1.3.6_paired.json").read_text())
    with pytest.raises(MetricsError, match="paired"):
        parse_fastp_json(data, paired=False)


def test_fastp_missing_values_not_invented() -> None:
    m = parse_fastp_json({"summary": {"before_filtering": {"total_reads": 10}}}, paired=False)
    assert m.reads_before == 10 and m.reads_after is None and m.retained_fraction is None
    assert m.q30_rate_before is None and m.duplication_rate is None
    with pytest.raises(MetricsError):
        parse_fastp_json({"nothing": 1}, paired=False)


def test_fastp_unreadable(tmp_path: Path) -> None:
    p = tmp_path / "x.json"
    p.write_text("{truncated")
    with pytest.raises(MetricsError):
        read_fastp_json(p, paired=True)
