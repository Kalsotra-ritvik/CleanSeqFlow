"""Plain-language explanations of pipeline stages and QC metrics.

This text is shown in the terminal while the pipeline runs, by the
``explain`` command, and in the HTML report.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class StepInfo:
    number: int
    title: str
    purpose: str
    tool: str
    details: str


#: The six user-facing steps. Internal stages map onto them via STAGE_TO_STEP.
STEPS: dict[int, StepInfo] = {
    1: StepInfo(
        1,
        "DOWNLOAD SRA DATA",
        "Fetch each sequencing run from the NCBI Sequence Read Archive (SRA).",
        "prefetch (SRA Toolkit)",
        "SRA stores runs in its own compressed archive format (.sra), not as FASTQ. prefetch "
        "downloads the archive, resumes interrupted transfers and verifies the file against "
        "NCBI's checksum. The .sra files are kept unless you ask for their removal.",
    ),
    2: StepInfo(
        2,
        "CONVERT SRA TO FASTQ",
        "Extract the reads from the SRA archive into FASTQ text files.",
        "fasterq-dump (SRA Toolkit) + pigz",
        "fasterq-dump --split-3 writes mate 1 and mate 2 of paired-end runs to separate files in "
        "the same order (so R1 line N pairs with R2 line N), and puts reads whose mate is missing "
        "into a separate file. The FASTQ files are then gzip-compressed. Conversion is a separate "
        "step from download: a finished download does not imply a successful conversion.",
    ),
    3: StepInfo(
        3,
        "RAW READ VALIDATION AND QC",
        "Confirm the FASTQ files are intact, then measure raw read quality before any cleaning.",
        "built-in FASTQ validator, FastQC, MultiQC",
        "Every record is checked for correct FASTQ structure; paired files must have the same "
        "number of reads with matching names. FastQC then summarises base quality, GC content, "
        "adapter content, duplication and more for each file, and MultiQC combines all samples "
        "into one report.",
    ),
    4: StepInfo(
        4,
        "FASTP PREPROCESSING",
        "Remove adapter sequence and filter reads according to the selected quality thresholds "
        "before downstream alignment.",
        "fastp",
        "fastp trims adapter read-through (common when the RNA fragment is shorter than the read) "
        "and discards reads that fail the quality, N-content or length filters. Paired reads are "
        "kept or discarded together so R1/R2 stay synchronised. Every removed read is counted "
        "by reason in the fastp report.",
    ),
    5: StepInfo(
        5,
        "POST-CLEANING VALIDATION AND QC",
        "Check the cleaned files and measure quality again so it can be compared with the raw data.",
        "built-in FASTQ validator, FastQC, MultiQC",
        "The same checks as step 3 are repeated on fastp's output. A second, separate MultiQC "
        "report combines the post-cleaning FastQC reports with the fastp statistics.",
    ),
    6: StepInfo(
        6,
        "FINAL OUTPUT AND REPORTS",
        "Publish validated cleaned FASTQ files and summarise what happened to every sample.",
        "built-in reporting",
        "Cleaned files are placed in 06_final_fastq (as hard links when possible, so no extra "
        "space is used), with a manifest, file inventory and checksums. The summary report "
        "compares measured before/after metrics and lists warnings and failures.",
    ),
}

STAGE_TO_STEP: dict[str, int] = {
    "download": 1,
    "convert": 2,
    "validate_raw": 3,
    "fastqc_raw": 3,
    "multiqc_raw": 3,
    "fastp": 4,
    "validate_clean": 5,
    "fastqc_clean": 5,
    "multiqc_clean": 5,
    "publish": 6,
    "report": 6,
}

STAGE_LABELS: dict[str, str] = {
    "download": "SRA download",
    "convert": "SRA-to-FASTQ conversion",
    "validate_raw": "raw FASTQ validation",
    "fastqc_raw": "FastQC (raw reads)",
    "multiqc_raw": "MultiQC (raw reads)",
    "fastp": "fastp preprocessing",
    "validate_clean": "cleaned FASTQ validation",
    "fastqc_clean": "FastQC (cleaned reads)",
    "multiqc_clean": "MultiQC (cleaned reads + fastp)",
    "publish": "final FASTQ publication",
    "report": "summary report",
}

#: What each FastQC module assesses, and how to read it for bulk RNA-seq.
QC_MODULES: dict[str, tuple[str, str]] = {
    "basic_statistics": (
        "Basic Statistics",
        "Number of reads, read length range and overall GC%. Check that the read count is what "
        "you expect for the run.",
    ),
    "per_base_quality": (
        "Per base sequence quality",
        "Distribution of Phred quality scores at each position along the read. Quality usually "
        "declines towards the 3' end on Illumina instruments; a gradual decline is normal.",
    ),
    "per_tile_quality": (
        "Per tile sequence quality",
        "Quality deviations by flow-cell tile. Usually not available for SRA-derived reads, "
        "because the tile information is not kept in the read names.",
    ),
    "per_sequence_quality": (
        "Per sequence quality scores",
        "Distribution of mean quality per read. A subpopulation of uniformly poor reads suggests "
        "a technical problem in part of the run.",
    ),
    "per_base_content": (
        "Per base sequence content",
        "Proportion of A/C/G/T at each position. RNA-seq libraries commonly FAIL this module at "
        "the first ~10-15 bases because random-hexamer priming is not perfectly random; this is "
        "a known library-preparation effect and normally not a reason to discard data.",
    ),
    "per_sequence_gc": (
        "Per sequence GC content",
        "GC% distribution compared with a theoretical normal distribution. Transcriptomes are not "
        "expected to follow it exactly; sharp extra peaks can indicate contamination or dominant "
        "transcripts (e.g. rRNA).",
    ),
    "per_base_n": (
        "Per base N content",
        "Percentage of uncalled bases (N) per position. Elevated values point to sequencing "
        "problems at specific cycles.",
    ),
    "length_distribution": (
        "Sequence length distribution",
        "Read lengths. Raw Illumina reads usually have one length; after trimming a range of "
        "lengths is expected, so FastQC's WARN here is normal for cleaned data.",
    ),
    "duplication": (
        "Sequence duplication levels",
        "How often identical sequences occur. FastQC tracks only sequences first seen among the "
        "first 100,000 reads of each file, so the value is an estimate. In RNA-seq "
        "high duplication is EXPECTED: highly expressed transcripts are sampled many times. PCR "
        "over-amplification and low-input libraries also raise it. Duplicates are therefore not "
        "removed by this pipeline; FastQC cannot tell biological from technical duplicates.",
    ),
    "overrepresented": (
        "Overrepresented sequences",
        "Sequences making up > 0.1% of reads. In RNA-seq these are often very abundant "
        "transcripts (mitochondrial, ribosomal RNA); adapter or primer hits indicate "
        "technical contamination.",
    ),
    "adapter_content": (
        "Adapter content",
        "Cumulative percentage of reads containing known adapter sequences by position. Rising "
        "curves mean that fragments were shorter than the read length. This should drop after "
        "fastp trimming.",
    ),
    "kmer_content": (
        "Kmer content",
        "Position-biased short sequences. Disabled by default in FastQC 0.12.",
    ),
}

FLAG_DISCLAIMER = (
    "FastQC PASS/WARN/FAIL flags compare each file with thresholds designed for generic, "
    "whole-genome DNA libraries. Several modules routinely WARN or FAIL for good RNA-seq data "
    "(per-base sequence content, duplication, GC content, overrepresented sequences). A flag is "
    "a prompt to look at the plot, not a decision about whether a sample is usable."
)

LIMITATIONS = [
    "Automated QC cannot establish that a sample is biologically valid; it summarises technical "
    "properties of the reads. Sample identity, contamination by other organisms, rRNA content, "
    "strandedness, RNA degradation and library complexity need further checks (e.g. alignment-"
    "based QC).",
    "Screening rules in the report are configurable heuristics chosen for this tool, not "
    "universal thresholds; adapt them to the study design.",
    "FastQC duplication estimates are based on sequences first seen in the first 100,000 reads "
    "of each file.",
    "fastp's duplication rate is computed on the input reads and is an estimate.",
    "Read names from fasterq-dump do not contain the original Illumina instrument information, "
    "so fastp's automatic polyG trimming (for NextSeq/NovaSeq two-colour chemistry) is usually "
    "not triggered; set fastp.trim_poly_g: on if your data come from such instruments.",
    "The pipeline does not perform alignment, quantification or differential expression.",
]


def stage_header(stage: str) -> StepInfo:
    return STEPS[STAGE_TO_STEP[stage]]


def render_explanation(topic: str | None = None) -> str:
    """Text for the ``explain`` command."""
    lines: list[str] = []
    if topic in (None, "steps", "all"):
        lines.append("PIPELINE STEPS")
        lines.append("==============")
        for info in STEPS.values():
            lines.append(f"\n[STEP {info.number}/{len(STEPS)}] {info.title}")
            lines.append(f"  Purpose : {info.purpose}")
            lines.append(f"  Tool    : {info.tool}")
            lines.append(f"  Details : {info.details}")
    if topic in (None, "qc", "all"):
        lines.append("\nFASTQC MODULES")
        lines.append("==============")
        for title, text in QC_MODULES.values():
            lines.append(f"\n* {title}\n  {text}")
        lines.append(f"\nAbout PASS/WARN/FAIL:\n  {FLAG_DISCLAIMER}")
    if topic in (None, "limitations", "all"):
        lines.append("\nLIMITATIONS")
        lines.append("===========")
        lines.extend(f"- {item}" for item in LIMITATIONS)
    if not lines:
        raise ValueError(f"unknown topic: {topic}")
    return "\n".join(lines) + "\n"
