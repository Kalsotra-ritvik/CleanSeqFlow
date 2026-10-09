"""Output directory structure. All paths the pipeline writes are defined here."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .validation import Sample

#: Marker identifying a directory as a Transcriptomics Preprocessor output directory.
MARKER_NAME = ".transcriptomics-preprocessor-output"


@dataclass(frozen=True)
class OutputLayout:
    """Locations of every output, relative to the run's output directory.

    ::

        <outdir>/
          01_sra/<accession>/            downloaded SRA archives (prefetch)
          02_raw_fastq/                  raw FASTQ converted from SRA (gzip)
            unpaired/                    mate-less reads from paired runs, if any
          03_raw_qc/fastqc/<sample>/     FastQC on raw reads
          03_raw_qc/multiqc/             MultiQC across raw FastQC reports
          04_fastp/<sample>/             fastp output FASTQ + HTML/JSON reports
          05_clean_qc/fastqc/<sample>/   FastQC on cleaned reads
          05_clean_qc/multiqc/           MultiQC across cleaned FastQC + fastp
          06_final_fastq/                cleaned FASTQ for downstream analysis
          07_reports/                    summary report, tables, comparisons
          logs/                          pipeline logs and per-command logs
          pipeline_info/                 run journal, configuration, versions
          tmp/                           default temporary directory
    """

    root: Path

    # Top-level directories -------------------------------------------------
    @property
    def sra_dir(self) -> Path:
        return self.root / "01_sra"

    @property
    def raw_fastq_dir(self) -> Path:
        return self.root / "02_raw_fastq"

    @property
    def unpaired_dir(self) -> Path:
        return self.raw_fastq_dir / "unpaired"

    @property
    def raw_qc_dir(self) -> Path:
        return self.root / "03_raw_qc"

    @property
    def fastp_dir(self) -> Path:
        return self.root / "04_fastp"

    @property
    def clean_qc_dir(self) -> Path:
        return self.root / "05_clean_qc"

    @property
    def final_dir(self) -> Path:
        return self.root / "06_final_fastq"

    @property
    def reports_dir(self) -> Path:
        return self.root / "07_reports"

    @property
    def logs_dir(self) -> Path:
        return self.root / "logs"

    @property
    def info_dir(self) -> Path:
        return self.root / "pipeline_info"

    @property
    def default_tmp_dir(self) -> Path:
        return self.root / "tmp"

    @property
    def state_file(self) -> Path:
        return self.info_dir / "state.json"

    @property
    def lock_file(self) -> Path:
        return self.info_dir / "run.lock"

    @property
    def marker_file(self) -> Path:
        return self.root / MARKER_NAME

    # Per-sample paths ----------------------------------------------------------
    def sra_accession_dir(self, sample: Sample) -> Path:
        return self.sra_dir / sample.accession

    def conversion_work_dir(self, sample: Sample) -> Path:
        return self.raw_fastq_dir / f".work_{sample.name}"

    def raw_fastq(self, sample: Sample) -> list[Path]:
        d = self.raw_fastq_dir
        if sample.paired:
            return [d / f"{sample.name}_R1.fastq.gz", d / f"{sample.name}_R2.fastq.gz"]
        return [d / f"{sample.name}.fastq.gz"]

    def unpaired_fastq(self, sample: Sample) -> Path:
        return self.unpaired_dir / f"{sample.name}_unpaired.fastq.gz"

    def raw_fastqc_dir(self, sample: Sample) -> Path:
        return self.raw_qc_dir / "fastqc" / sample.name

    def clean_fastqc_dir(self, sample: Sample) -> Path:
        return self.clean_qc_dir / "fastqc" / sample.name

    def fastp_sample_dir(self, sample: Sample) -> Path:
        return self.fastp_dir / sample.name

    def clean_fastq(self, sample: Sample) -> list[Path]:
        d = self.fastp_sample_dir(sample)
        if sample.paired:
            return [d / f"{sample.name}_R1.fastq.gz", d / f"{sample.name}_R2.fastq.gz"]
        return [d / f"{sample.name}.fastq.gz"]

    def fastp_json(self, sample: Sample) -> Path:
        return self.fastp_sample_dir(sample) / f"{sample.name}.fastp.json"

    def fastp_html(self, sample: Sample) -> Path:
        return self.fastp_sample_dir(sample) / f"{sample.name}.fastp.html"

    def final_fastq(self, sample: Sample) -> list[Path]:
        d = self.final_dir
        if sample.paired:
            return [d / f"{sample.name}_R1.fastq.gz", d / f"{sample.name}_R2.fastq.gz"]
        return [d / f"{sample.name}.fastq.gz"]

    def command_log(self, sample: Sample | None, stage: str) -> Path:
        sub = sample.name if sample is not None else "_aggregate"
        return self.logs_dir / "commands" / sub / f"{stage}.log"

    def tmp_for(self, tmp_root: Path, sample: Sample, stage: str) -> Path:
        return tmp_root / f"{sample.name}.{stage}"

    @staticmethod
    def fastqc_stem(fastq: Path) -> str:
        """FastQC names its outputs after the input file minus FASTQ/gzip extensions."""
        name = fastq.name
        for ext in (".gz", ".bz2"):
            name = name.removesuffix(ext)
        for ext in (".fastq", ".fq"):
            name = name.removesuffix(ext)
        return name

    def rel(self, path: Path) -> str:
        """Path relative to the output root (for portable records)."""
        try:
            return str(Path(path).resolve().relative_to(self.root.resolve()))
        except ValueError:
            return str(path)

    def abs(self, rel: str) -> Path:
        p = Path(rel)
        return p if p.is_absolute() else self.root / p
