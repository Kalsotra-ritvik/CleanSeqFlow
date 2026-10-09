"""Before-versus-after comparison and automated screening.

Rules followed here:

* Only compare like with like: each row states its definition and source, and
  the same source is used for the "before" and "after" value.
* Reads and read pairs are reported separately. For paired-end data, "reads"
  counts R1 + R2 records, "read pairs" counts R1 records.
* Missing values stay ``None`` ("not available"); no value is estimated.
* The direction of a change is stated from the numbers ("decreased from A to
  B"); no improvement is claimed that the numbers do not show.
* Screening flags are heuristics with stated, configurable thresholds. They
  mark samples for human review; they are not biological verdicts.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from .config import PipelineConfig

#: FastQC 0.12.1 default WARN level for Adapter Content (Configuration/limits.txt).
FASTQC_ADAPTER_WARN = 5.0
#: Differences smaller than this (percentage points) are described as "essentially unchanged".
NEGLIGIBLE_PERCENT_POINTS = 0.5


@dataclass
class ComparisonRow:
    metric: str
    definition: str
    source: str
    before: float | int | str | None
    after: float | int | str | None
    unit: str = ""

    @property
    def change(self) -> str:
        b, a = self.before, self.after
        if b is None or a is None:
            return "not available"
        if not isinstance(b, (int, float)) or not isinstance(a, (int, float)):
            return "unchanged" if a == b else f"{b} → {a}"
        if a == b:
            return "unchanged"
        direction = "increased" if a > b else "decreased"
        return f"{direction} by {abs(a - b):,.2f}{self.unit}".replace(".00", "")

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["change"] = self.change
        return d


@dataclass
class SampleComparison:
    sample: str
    paired: bool
    rows: list[ComparisonRow] = field(default_factory=list)
    removed: dict[str, int | None] = field(default_factory=dict)
    flags: list[str] = field(default_factory=list)
    observations: list[str] = field(default_factory=list)
    consistency_issues: list[str] = field(default_factory=list)

    @property
    def screening_status(self) -> str:
        if self.consistency_issues:
            return "INCONSISTENT"
        return "REVIEW" if self.flags else "NO_FLAGS"

    def to_dict(self) -> dict[str, Any]:
        return {
            "sample": self.sample,
            "paired": self.paired,
            "screening_status": self.screening_status,
            "rows": [r.to_dict() for r in self.rows],
            "reads_removed_by_reason": self.removed,
            "screening_flags": self.flags,
            "observations": self.observations,
            "consistency_issues": self.consistency_issues,
        }


def _pct(x: float | None) -> float | None:
    return None if x is None else round(100.0 * x, 2)


def _sum(values: list[int | None]) -> int | None:
    if not values or any(v is None for v in values):
        return None
    return sum(v for v in values if v is not None)


def _fastqc(m: dict[str, Any] | None, key: str) -> Any:
    return None if not m else m.get(key)


def build_comparison(
    sample: str,
    paired: bool,
    raw_validation: dict[str, Any] | None,
    clean_validation: dict[str, Any] | None,
    fastp: dict[str, Any] | None,
    fastqc_raw: dict[str, Any] | None,
    fastqc_clean: dict[str, Any] | None,
    config: PipelineConfig,
) -> SampleComparison:
    """Assemble the comparison for one sample from recorded stage metrics.

    ``*_validation`` map mate labels (``R1``/``R2``/``reads``) to FASTQ
    validator statistics; ``fastqc_*`` map mate labels to FastQC metrics;
    ``fastp`` is the dictionary from :class:`~.metrics.FastpMetrics`.
    """
    cmp = SampleComparison(sample=sample, paired=paired)
    mates = ["R1", "R2"] if paired else ["reads"]
    rv, cv = raw_validation or {}, clean_validation or {}
    fp = fastp or {}

    reads_before = _sum([rv.get(m, {}).get("records") for m in mates]) if rv else None
    reads_after = _sum([cv.get(m, {}).get("records") for m in mates]) if cv else None
    cmp.rows.append(
        ComparisonRow(
            "Reads",
            "Individual reads in all files of the sample (R1 + R2 for paired-end)",
            "FASTQ validator (exact count)",
            reads_before,
            reads_after,
        )
    )
    if paired:
        pairs_b = rv.get("R1", {}).get("records") if rv else None
        pairs_a = cv.get("R1", {}).get("records") if cv else None
        cmp.rows.append(
            ComparisonRow(
                "Read pairs",
                "Read pairs (records in R1; R2 must match)",
                "FASTQ validator (exact count)",
                pairs_b,
                pairs_a,
            )
        )
    retained = (
        round(100.0 * reads_after / reads_before, 2)
        if reads_before and reads_after is not None
        else None
    )
    cmp.rows.append(
        ComparisonRow(
            "Reads retained",
            "Reads after cleaning as a percentage of reads before",
            "derived from the counts above",
            100.0 if reads_before else None,
            retained,
            "%",
        )
    )
    cmp.rows.append(
        ComparisonRow(
            "Bases",
            "Total sequenced bases over all reads",
            "fastp JSON summary",
            fp.get("bases_before"),
            fp.get("bases_after"),
        )
    )
    cmp.rows.append(
        ComparisonRow(
            "Q20 bases",
            "Percentage of bases with Phred quality >= 20",
            "fastp JSON summary",
            _pct(fp.get("q20_rate_before")),
            _pct(fp.get("q20_rate_after")),
            "%",
        )
    )
    cmp.rows.append(
        ComparisonRow(
            "Q30 bases",
            "Percentage of bases with Phred quality >= 30",
            "fastp JSON summary",
            _pct(fp.get("q30_rate_before")),
            _pct(fp.get("q30_rate_after")),
            "%",
        )
    )
    cmp.rows.append(
        ComparisonRow(
            "GC content",
            "Percentage of G and C among all bases",
            "fastp JSON summary",
            _pct(fp.get("gc_before")),
            _pct(fp.get("gc_after")),
            "%",
        )
    )
    for mate, key in (("R1", "r1"), ("R2", "r2")):
        if mate == "R2" and not paired:
            continue
        cmp.rows.append(
            ComparisonRow(
                f"Mean read length ({mate if paired else 'reads'})",
                "Mean length of reads",
                "fastp JSON summary",
                fp.get(f"{key}_mean_length_before"),
                fp.get(f"{key}_mean_length_after"),
                " bp",
            )
        )

    for mate in mates:
        fr = (fastqc_raw or {}).get(mate)
        fc = (fastqc_clean or {}).get(mate)
        label = mate if paired else "reads"
        cmp.rows.append(
            ComparisonRow(
                f"Max adapter content ({label})",
                "Highest cumulative % of reads with a known adapter at any position, any adapter type",
                "FastQC Adapter Content",
                _fastqc(fr, "max_adapter_percent"),
                _fastqc(fc, "max_adapter_percent"),
                "%",
            )
        )
        cmp.rows.append(
            ComparisonRow(
                f"Lowest per-position mean quality ({label})",
                "Minimum over read positions of the mean Phred quality",
                "FastQC Per base sequence quality",
                _fastqc(fr, "min_position_mean_quality"),
                _fastqc(fc, "min_position_mean_quality"),
            )
        )
        cmp.rows.append(
            ComparisonRow(
                f"Estimated % remaining after deduplication ({label})",
                "FastQC 'Total Deduplicated Percentage' (estimate; lower = more duplication)",
                "FastQC Sequence Duplication Levels",
                _fastqc(fr, "deduplicated_percent"),
                _fastqc(fc, "deduplicated_percent"),
                "%",
            )
        )
        cmp.rows.append(
            ComparisonRow(
                f"FastQC FAIL flags ({label})",
                "Number of FastQC modules flagged FAIL",
                "FastQC summary.txt",
                _count_flags(fr, "FAIL"),
                _count_flags(fc, "FAIL"),
            )
        )

    if fp:
        cmp.removed = {
            "low_quality": fp.get("low_quality_reads"),
            "too_many_N": fp.get("too_many_n_reads"),
            "too_short": fp.get("too_short_reads"),
            "too_long": fp.get("too_long_reads"),
            "adapter_dimer": fp.get("adapter_dimer_reads"),
        }
        if fp.get("adapter_trimmed_reads") is not None:
            cmp.observations.append(
                f"fastp trimmed adapter sequence from {fp['adapter_trimmed_reads']:,} reads "
                f"({fp.get('adapter_trimmed_bases') or 0:,} bases)."
            )
        elif config["fastp.adapter_trimming"]:
            raw_max = max(
                (
                    v
                    for m in mates
                    if (v := _fastqc((fastqc_raw or {}).get(m), "max_adapter_percent")) is not None
                ),
                default=None,
            )
            msg = (
                "fastp reported no adapter trimming for this sample (its JSON report has no "
                "adapter-cutting section)."
            )
            if raw_max is not None and raw_max > FASTQC_ADAPTER_WARN:
                msg += (
                    f" FastQC found up to {raw_max:.2f}% adapter content in the raw reads, so "
                    "automatic adapter detection probably failed (it is less reliable for "
                    "single-end data and small files). Set the adapter explicitly, e.g. "
                    "--adapter-r1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA for Illumina TruSeq."
                )
            cmp.observations.append(msg)
        _check_consistency(cmp, fp, reads_before, reads_after)

    _observe_adapters(cmp, mates, fastqc_raw, fastqc_clean, paired)
    _screen(cmp, config, retained, paired, mates, fp, cv, fastqc_clean)
    return cmp


def _count_flags(m: dict[str, Any] | None, flag: str) -> int | None:
    if not m or not m.get("module_status"):
        return None
    return sum(1 for v in m["module_status"].values() if v == flag)


def _check_consistency(
    cmp: SampleComparison, fp: dict[str, Any], reads_before: int | None, reads_after: int | None
) -> None:
    if (
        reads_before is not None
        and fp.get("reads_before") is not None
        and reads_before != fp["reads_before"]
    ):
        cmp.consistency_issues.append(
            f"fastp counted {fp['reads_before']:,} input reads but the raw files contain {reads_before:,}"
        )
    if (
        reads_after is not None
        and fp.get("reads_after") is not None
        and reads_after != fp["reads_after"]
    ):
        cmp.consistency_issues.append(
            f"fastp reports {fp['reads_after']:,} output reads but the cleaned files contain {reads_after:,}"
        )
    parts = [
        fp.get(k)
        for k in (
            "passed_filter_reads",
            "low_quality_reads",
            "too_many_n_reads",
            "too_short_reads",
            "too_long_reads",
            "adapter_dimer_reads",
        )
    ]
    if fp.get("reads_before") is not None and all(p is not None for p in parts):
        total = sum(p for p in parts if p is not None)
        if total != fp["reads_before"]:
            cmp.consistency_issues.append(
                f"fastp filtering categories add up to {total:,}, not {fp['reads_before']:,} input reads"
            )


def _observe_adapters(
    cmp: SampleComparison,
    mates: list[str],
    raw: dict[str, Any] | None,
    clean: dict[str, Any] | None,
    paired: bool,
) -> None:
    for mate in mates:
        b = _fastqc((raw or {}).get(mate), "max_adapter_percent")
        a = _fastqc((clean or {}).get(mate), "max_adapter_percent")
        label = mate if paired else "reads"
        if b is None or a is None:
            cmp.observations.append(f"Adapter content ({label}): not available for comparison.")
        elif abs(a - b) < NEGLIGIBLE_PERCENT_POINTS and a != b:
            cmp.observations.append(
                f"Adapter content ({label}) is essentially unchanged ({b:.2f}% before, {a:.2f}% after; "
                "FastQC, max over positions)."
            )
        elif a < b:
            cmp.observations.append(
                f"Adapter content ({label}) decreased from {b:.2f}% to {a:.2f}% (FastQC, max over positions)."
            )
        elif a > b:
            cmp.observations.append(
                f"Adapter content ({label}) INCREASED from {b:.2f}% to {a:.2f}%; inspect the reports."
            )
        else:
            cmp.observations.append(f"Adapter content ({label}) unchanged at {a:.2f}%.")


def _screen(
    cmp: SampleComparison,
    config: PipelineConfig,
    retained: float | None,
    paired: bool,
    mates: list[str],
    fp: dict[str, Any],
    cv: dict[str, Any],
    fastqc_clean: dict[str, Any] | None,
) -> None:
    min_ret = float(config["screening.min_retained_percent"])
    if retained is not None and retained < min_ret:
        cmp.flags.append(
            f"Only {retained:.2f}% of reads were retained (screening threshold {min_ret:g}%). "
            "Check the fastp report for the main removal reason and whether the raw data or the "
            "parameters explain it."
        )
    min_reads = int(config["screening.min_clean_reads"])
    clean_units = cv.get("R1", {}).get("records") if paired else cv.get("reads", {}).get("records")
    unit = "read pairs" if paired else "reads"
    if clean_units is not None and clean_units < min_reads:
        cmp.flags.append(
            f"{clean_units:,} {unit} after cleaning is below the screening threshold of "
            f"{min_reads:,}. Whether this depth is sufficient depends on the study design."
        )
    max_adapt = float(config["screening.max_clean_adapter_percent"])
    for mate in mates:
        a = _fastqc((fastqc_clean or {}).get(mate), "max_adapter_percent")
        if a is not None and a > max_adapt:
            cmp.flags.append(
                f"Residual adapter content in cleaned {mate} is {a:.2f}% (screening threshold "
                f"{max_adapt:g}%). Consider specifying the adapter sequence explicitly."
            )
    q30 = fp.get("q30_rate_after")
    min_q30 = float(config["screening.min_clean_q30_rate"])
    if q30 is not None and q30 < min_q30:
        cmp.flags.append(
            f"Q30 rate after cleaning is {100 * q30:.2f}% (screening threshold {100 * min_q30:g}%)."
        )
