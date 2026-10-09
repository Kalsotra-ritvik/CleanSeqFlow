"""Transcriptomics Preprocessor.

Automated, resumable preprocessing and quality control of bulk RNA-seq reads,
from NCBI SRA run accessions to quality-controlled, cleaned FASTQ files.

The package orchestrates established tools (SRA Toolkit, FastQC, fastp,
MultiQC); it does not perform alignment, quantification or differential
expression.
"""

__version__ = "0.1.0"
__all__ = ["__version__"]
