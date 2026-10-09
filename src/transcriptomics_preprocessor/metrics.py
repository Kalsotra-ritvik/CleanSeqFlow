"""Extraction of QC metrics from the tools' own machine-readable outputs.

Sources and definitions
-----------------------
* **FastQC** (``*_fastqc.zip``): ``summary.txt`` (PASS/WARN/FAIL per module)
  and ``fastqc_data.txt`` (module tables). Metrics are per *file*, i.e. per
  read mate for paired-end data.
* **fastp** (``*.fastp.json``): ``summary``, ``filtering_result``,
  ``adapter_cutting`` and ``duplication`` sections. fastp counts *reads*: for
  paired-end data each pair contributes two reads, so read pairs =
  ``total_reads / 2``.

A metric that cannot be found is returned as ``None`` and must be shown as
"not available"; values are never estimated or filled in.
"""

from __future__ import annotations

import json
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

#: FastQC modules in report order, with the short keys used in this package.
FASTQC_MODULES: dict[str, str] = {
    "Basic Statistics": "basic_statistics",
    "Per base sequence quality": "per_base_quality",
    "Per tile sequence quality": "per_tile_quality",
    "Per sequence quality scores": "per_sequence_quality",
    "Per base sequence content": "per_base_content",
    "Per sequence GC content": "per_sequence_gc",
    "Per base N content": "per_base_n",
    "Sequence Length Distribution": "length_distribution",
    "Sequence Duplication Levels": "duplication",
    "Overrepresented sequences": "overrepresented",
    "Adapter Content": "adapter_content",
    "Kmer Content": "kmer_content",
}


class MetricsError(Exception):
    """A report file exists but cannot be parsed."""


@dataclass
class FastqcMetrics:
    """Selected values from one FastQC report (one FASTQ file)."""

    file: str
    total_sequences: int | None = None
    poor_quality_sequences: int | None = None
    sequence_length: str | None = None
    gc_percent: float | None = None
    module_status: dict[str, str] = field(default_factory=dict)
    #: Lowest per-position *mean* quality (Per base sequence quality module).
    min_position_mean_quality: float | None = None
    #: Mean of the per-position mean qualities (unweighted across positions).
    mean_position_mean_quality: float | None = None
    #: Most frequent per-read mean quality (Per sequence quality scores).
    modal_read_quality: float | None = None
    #: Highest cumulative adapter percentage at any position, over all adapter types.
    max_adapter_percent: float | None = None
    adapter_type_at_max: str | None = None
    #: FastQC "Total Deduplicated Percentage" (estimate on the first 100,000 sequences).
    deduplicated_percent: float | None = None
    overrepresented_count: int | None = None
    max_n_percent: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return dict(self.__dict__)


def _modules(text: str) -> dict[str, tuple[str, list[list[str]]]]:
    """Split fastqc_data.txt into ``{module name: (status, rows)}``."""
    modules: dict[str, tuple[str, list[list[str]]]] = {}
    current: str | None = None
    status = ""
    rows: list[list[str]] = []
    for line in text.splitlines():
        if line.startswith(">>END_MODULE"):
            if current is not None:
                modules[current] = (status, rows)
            current, rows = None, []
        elif line.startswith(">>"):
            parts = line[2:].split("\t")
            current = parts[0]
            status = parts[1] if len(parts) > 1 else ""
            rows = []
        elif current is not None:
            rows.append(line.split("\t"))
    return modules


def _float(x: str) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _int(x: str) -> int | None:
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def parse_fastqc_data(text: str, file_label: str = "") -> FastqcMetrics:
    """Parse the contents of ``fastqc_data.txt``."""
    if not text.startswith("##FastQC"):
        raise MetricsError("not a FastQC data file (missing '##FastQC' header)")
    mods = _modules(text)
    m = FastqcMetrics(file=file_label)
    for name, (status, _rows) in mods.items():
        key = FASTQC_MODULES.get(name)
        if key:
            m.module_status[key] = status.upper()

    if "Basic Statistics" in mods:
        stats = {
            r[0]: r[1]
            for r in mods["Basic Statistics"][1]
            if len(r) >= 2 and not r[0].startswith("#")
        }
        m.total_sequences = _int(stats.get("Total Sequences", ""))
        m.poor_quality_sequences = _int(stats.get("Sequences flagged as poor quality", ""))
        m.sequence_length = stats.get("Sequence length")
        m.gc_percent = _float(stats.get("%GC", ""))
        if not m.file:
            m.file = stats.get("Filename", "")

    if "Per base sequence quality" in mods:
        means = [
            v
            for r in mods["Per base sequence quality"][1]
            if len(r) >= 2 and not r[0].startswith("#") and (v := _float(r[1])) is not None
        ]
        if means:
            m.min_position_mean_quality = min(means)
            m.mean_position_mean_quality = sum(means) / len(means)

    if "Per sequence quality scores" in mods:
        best: tuple[float, float] | None = None
        for r in mods["Per sequence quality scores"][1]:
            if len(r) >= 2 and not r[0].startswith("#"):
                q, c = _float(r[0]), _float(r[1])
                if q is not None and c is not None and (best is None or c > best[1]):
                    best = (q, c)
        if best is not None:
            m.modal_read_quality = best[0]

    if "Adapter Content" in mods:
        rows = mods["Adapter Content"][1]
        header = next((r for r in rows if r and r[0].startswith("#")), None)
        best_val: float | None = None
        best_type: str | None = None
        for r in rows:
            if not r or r[0].startswith("#"):
                continue
            for i, cell in enumerate(r[1:], start=1):
                v = _float(cell)
                if v is not None and (best_val is None or v > best_val):
                    best_val = v
                    best_type = header[i] if header and i < len(header) else None
        m.max_adapter_percent = best_val
        m.adapter_type_at_max = best_type

    if "Sequence Duplication Levels" in mods:
        for r in mods["Sequence Duplication Levels"][1]:
            if r and r[0] == "#Total Deduplicated Percentage" and len(r) >= 2:
                m.deduplicated_percent = _float(r[1])

    if "Overrepresented sequences" in mods:
        m.overrepresented_count = sum(
            1
            for r in mods["Overrepresented sequences"][1]
            if r and r[0] and not r[0].startswith("#")
        )

    if "Per base N content" in mods:
        ns = [
            v
            for r in mods["Per base N content"][1]
            if len(r) >= 2 and not r[0].startswith("#") and (v := _float(r[1])) is not None
        ]
        m.max_n_percent = max(ns) if ns else None
    return m


def read_fastqc_zip(zip_path: Path) -> FastqcMetrics:
    """Read ``fastqc_data.txt`` and ``summary.txt`` from a FastQC zip archive."""
    try:
        with zipfile.ZipFile(zip_path) as zf:
            data_name = next((n for n in zf.namelist() if n.endswith("/fastqc_data.txt")), None)
            if data_name is None:
                raise MetricsError(f"{zip_path.name}: fastqc_data.txt not found in archive")
            text = zf.read(data_name).decode("utf-8", errors="replace")
            summary_name = data_name.replace("fastqc_data.txt", "summary.txt")
            summary = (
                zf.read(summary_name).decode("utf-8", errors="replace")
                if summary_name in zf.namelist()
                else ""
            )
    except zipfile.BadZipFile as exc:
        raise MetricsError(f"{zip_path.name}: not a valid zip archive ({exc})") from exc
    metrics = parse_fastqc_data(text)
    # summary.txt holds the authoritative PASS/WARN/FAIL flags shown in the HTML report.
    for line in summary.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1] in FASTQC_MODULES:
            metrics.module_status[FASTQC_MODULES[parts[1]]] = parts[0].strip().upper()
    return metrics


@dataclass
class FastpMetrics:
    """Selected values from one fastp JSON report."""

    paired: bool
    fastp_version: str | None
    reads_before: int | None
    reads_after: int | None
    bases_before: int | None
    bases_after: int | None
    q20_rate_before: float | None
    q20_rate_after: float | None
    q30_rate_before: float | None
    q30_rate_after: float | None
    gc_before: float | None
    gc_after: float | None
    r1_mean_length_before: float | None
    r1_mean_length_after: float | None
    r2_mean_length_before: float | None
    r2_mean_length_after: float | None
    passed_filter_reads: int | None
    low_quality_reads: int | None
    too_many_n_reads: int | None
    too_short_reads: int | None
    too_long_reads: int | None
    adapter_dimer_reads: int | None
    adapter_trimmed_reads: int | None
    adapter_trimmed_bases: int | None
    duplication_rate: float | None
    insert_size_peak: int | None
    command: str | None

    @property
    def pairs_before(self) -> int | None:
        return self.reads_before // 2 if self.paired and self.reads_before is not None else None

    @property
    def pairs_after(self) -> int | None:
        return self.reads_after // 2 if self.paired and self.reads_after is not None else None

    @property
    def retained_fraction(self) -> float | None:
        if not self.reads_before or self.reads_after is None:
            return None
        return self.reads_after / self.reads_before

    def to_dict(self) -> dict[str, Any]:
        d = dict(self.__dict__)
        d["pairs_before"] = self.pairs_before
        d["pairs_after"] = self.pairs_after
        d["retained_fraction"] = self.retained_fraction
        return d


def _get(d: Any, *keys: str) -> Any:
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def _num(x: Any) -> Any:
    return x if isinstance(x, (int, float)) and not isinstance(x, bool) else None


def parse_fastp_json(data: dict[str, Any], paired: bool) -> FastpMetrics:
    """Extract metrics from a parsed fastp JSON report."""
    if not isinstance(data, dict) or "summary" not in data:
        raise MetricsError("not a fastp JSON report (no 'summary' section)")
    seq_desc = str(_get(data, "summary", "sequencing") or "")
    if seq_desc:
        reported_paired = seq_desc.startswith("paired")
        if reported_paired != paired:
            raise MetricsError(
                f"fastp report describes '{seq_desc}' but the sample is "
                f"{'paired' if paired else 'single'}-end"
            )
    b = _get(data, "summary", "before_filtering") or {}
    a = _get(data, "summary", "after_filtering") or {}
    f = data.get("filtering_result") or {}
    ad = data.get("adapter_cutting") or {}
    return FastpMetrics(
        paired=paired,
        fastp_version=_get(data, "summary", "fastp_version"),
        reads_before=_num(b.get("total_reads")),
        reads_after=_num(a.get("total_reads")),
        bases_before=_num(b.get("total_bases")),
        bases_after=_num(a.get("total_bases")),
        q20_rate_before=_num(b.get("q20_rate")),
        q20_rate_after=_num(a.get("q20_rate")),
        q30_rate_before=_num(b.get("q30_rate")),
        q30_rate_after=_num(a.get("q30_rate")),
        gc_before=_num(b.get("gc_content")),
        gc_after=_num(a.get("gc_content")),
        r1_mean_length_before=_num(b.get("read1_mean_length")),
        r1_mean_length_after=_num(a.get("read1_mean_length")),
        r2_mean_length_before=_num(b.get("read2_mean_length")),
        r2_mean_length_after=_num(a.get("read2_mean_length")),
        passed_filter_reads=_num(f.get("passed_filter_reads")),
        low_quality_reads=_num(f.get("low_quality_reads")),
        too_many_n_reads=_num(f.get("too_many_N_reads")),
        too_short_reads=_num(f.get("too_short_reads")),
        too_long_reads=_num(f.get("too_long_reads")),
        adapter_dimer_reads=_num(f.get("adapter_dimer_reads")),
        adapter_trimmed_reads=_num(ad.get("adapter_trimmed_reads")),
        adapter_trimmed_bases=_num(ad.get("adapter_trimmed_bases")),
        duplication_rate=_num(_get(data, "duplication", "rate")),
        insert_size_peak=_num(_get(data, "insert_size", "peak")),
        command=data.get("command") if isinstance(data.get("command"), str) else None,
    )


def read_fastp_json(path: Path, paired: bool) -> FastpMetrics:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise MetricsError(f"cannot read fastp JSON {path}: {exc}") from exc
    return parse_fastp_json(data, paired)
