# How the pipeline works

This text is identical to the output of `transcriptomics-preprocess explain all`
(a test keeps them in sync).

```text
PIPELINE STEPS
==============

[STEP 1/6] DOWNLOAD SRA DATA
  Purpose : Fetch each sequencing run from the NCBI Sequence Read Archive (SRA).
  Tool    : prefetch (SRA Toolkit)
  Details : SRA stores runs in its own compressed archive format (.sra), not as FASTQ. prefetch downloads the archive, resumes interrupted transfers and verifies the file against NCBI's checksum. The .sra files are kept unless you ask for their removal.

[STEP 2/6] CONVERT SRA TO FASTQ
  Purpose : Extract the reads from the SRA archive into FASTQ text files.
  Tool    : fasterq-dump (SRA Toolkit) + pigz
  Details : fasterq-dump --split-3 writes mate 1 and mate 2 of paired-end runs to separate files in the same order (so R1 line N pairs with R2 line N), and puts reads whose mate is missing into a separate file. The FASTQ files are then gzip-compressed. Conversion is a separate step from download: a finished download does not imply a successful conversion.

[STEP 3/6] RAW READ VALIDATION AND QC
  Purpose : Confirm the FASTQ files are intact, then measure raw read quality before any cleaning.
  Tool    : built-in FASTQ validator, FastQC, MultiQC
  Details : Every record is checked for correct FASTQ structure; paired files must have the same number of reads with matching names. FastQC then summarises base quality, GC content, adapter content, duplication and more for each file, and MultiQC combines all samples into one report.

[STEP 4/6] FASTP PREPROCESSING
  Purpose : Remove adapter sequence and filter reads according to the selected quality thresholds before downstream alignment.
  Tool    : fastp
  Details : fastp trims adapter read-through (common when the RNA fragment is shorter than the read) and discards reads that fail the quality, N-content or length filters. Paired reads are kept or discarded together so R1/R2 stay synchronised. Every removed read is counted by reason in the fastp report.

[STEP 5/6] POST-CLEANING VALIDATION AND QC
  Purpose : Check the cleaned files and measure quality again so it can be compared with the raw data.
  Tool    : built-in FASTQ validator, FastQC, MultiQC
  Details : The same checks as step 3 are repeated on fastp's output. A second, separate MultiQC report combines the post-cleaning FastQC reports with the fastp statistics.

[STEP 6/6] FINAL OUTPUT AND REPORTS
  Purpose : Publish validated cleaned FASTQ files and summarise what happened to every sample.
  Tool    : built-in reporting
  Details : Cleaned files are placed in 06_final_fastq (as hard links when possible, so no extra space is used), with a manifest, file inventory and checksums. The summary report compares measured before/after metrics and lists warnings and failures.

FASTQC MODULES
==============

* Basic Statistics
  Number of reads, read length range and overall GC%. Check that the read count is what you expect for the run.

* Per base sequence quality
  Distribution of Phred quality scores at each position along the read. Quality usually declines towards the 3' end on Illumina instruments; a gradual decline is normal.

* Per tile sequence quality
  Quality deviations by flow-cell tile. Usually not available for SRA-derived reads, because the tile information is not kept in the read names.

* Per sequence quality scores
  Distribution of mean quality per read. A subpopulation of uniformly poor reads suggests a technical problem in part of the run.

* Per base sequence content
  Proportion of A/C/G/T at each position. RNA-seq libraries commonly FAIL this module at the first ~10-15 bases because random-hexamer priming is not perfectly random; this is a known library-preparation effect and normally not a reason to discard data.

* Per sequence GC content
  GC% distribution compared with a theoretical normal distribution. Transcriptomes are not expected to follow it exactly; sharp extra peaks can indicate contamination or dominant transcripts (e.g. rRNA).

* Per base N content
  Percentage of uncalled bases (N) per position. Elevated values point to sequencing problems at specific cycles.

* Sequence length distribution
  Read lengths. Raw Illumina reads usually have one length; after trimming a range of lengths is expected, so FastQC's WARN here is normal for cleaned data.

* Sequence duplication levels
  How often identical sequences occur. FastQC tracks only sequences first seen among the first 100,000 reads of each file, so the value is an estimate. In RNA-seq high duplication is EXPECTED: highly expressed transcripts are sampled many times. PCR over-amplification and low-input libraries also raise it. Duplicates are therefore not removed by this pipeline; FastQC cannot tell biological from technical duplicates.

* Overrepresented sequences
  Sequences making up > 0.1% of reads. In RNA-seq these are often very abundant transcripts (mitochondrial, ribosomal RNA); adapter or primer hits indicate technical contamination.

* Adapter content
  Cumulative percentage of reads containing known adapter sequences by position. Rising curves mean that fragments were shorter than the read length. This should drop after fastp trimming.

* Kmer content
  Position-biased short sequences. Disabled by default in FastQC 0.12.

About PASS/WARN/FAIL:
  FastQC PASS/WARN/FAIL flags compare each file with thresholds designed for generic, whole-genome DNA libraries. Several modules routinely WARN or FAIL for good RNA-seq data (per-base sequence content, duplication, GC content, overrepresented sequences). A flag is a prompt to look at the plot, not a decision about whether a sample is usable.

LIMITATIONS
===========
- Automated QC cannot establish that a sample is biologically valid; it summarises technical properties of the reads. Sample identity, contamination by other organisms, rRNA content, strandedness, RNA degradation and library complexity need further checks (e.g. alignment-based QC).
- Screening rules in the report are configurable heuristics chosen for this tool, not universal thresholds; adapt them to the study design.
- FastQC duplication estimates are based on sequences first seen in the first 100,000 reads of each file.
- fastp's duplication rate is computed on the input reads and is an estimate.
- Read names from fasterq-dump do not contain the original Illumina instrument information, so fastp's automatic polyG trimming (for NextSeq/NovaSeq two-colour chemistry) is usually not triggered; set fastp.trim_poly_g: on if your data come from such instruments.
- The pipeline does not perform alignment, quantification or differential expression.
```
