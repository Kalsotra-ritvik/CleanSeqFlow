# User guide

## 1. Prepare the input

Find the SRA **run** accessions for your study. For a GEO series (e.g. GSE243105):
open the NCBI SRA Run Selector (<https://www.ncbi.nlm.nih.gov/Traces/study/>), search for
the GEO accession, and download the metadata table. Copy the `Run` column into
`accession` and `LibraryLayout` (PAIRED/SINGLE) into `layout`, and choose a short,
meaningful `sample` name for each run.

```csv
sample,accession,layout,condition
siOGT_rep1,SRR26053372,paired,siOGT
siOGT_rep2,SRR26053373,paired,siOGT
siCtrl_rep1,SRR26053375,paired,siControl
```

Check it:

```bash
transcriptomics-preprocess check-input --metadata metadata.csv
```

This checks format only; whether each accession exists at NCBI is checked when it is
downloaded.

## 2. Preview the run

```bash
transcriptomics-preprocess run --metadata metadata.csv --outdir RNAseq_results --threads 8 --dry-run
```

The dry run prints the validated settings, free disk space, detected tool versions and the
exact commands that would be executed. Nothing is downloaded or created.

## 3. Run

```bash
transcriptomics-preprocess run --metadata metadata.csv --outdir RNAseq_results --threads 8
```

Each step prints a block like this (values come from the actual run):

```text
[STEP 4/6] FASTP PREPROCESSING — fastp preprocessing
Purpose: Remove adapter sequence and filter reads according to the selected quality thresholds before downstream alignment.
Sample: ctrl_1 (1 of 3)
Tool: fastp 1.3.6
Input:
  02_raw_fastq/ctrl_1_R1.fastq.gz
  02_raw_fastq/ctrl_1_R2.fastq.gz
Parameters:
  Threads: 4
  Adapter trimming: R1/R2 overlap analysis plus adapter sequence auto-detection (...)
  Quality filter: discard a read if more than 40% of its bases have Phred quality < 15, ...
  Length filter: discard reads shorter than 15 bases after trimming
  ...
  Deduplication: OFF (duplicates are expected in RNA-seq and are kept)
Output:
  04_fastp/ctrl_1/ctrl_1_R1.fastq.gz
  ...
Status: RUNNING
fastp preprocessing [ctrl_1]: COMPLETED (3.1 s)
```

The purpose text is shown the first time each step appears (every time with `--verbose`).
`--quiet` shows only warnings, errors and the final summary; everything is still written to
`logs/pipeline.log`.

### Threads and parallel samples

`--threads N` is the total budget. `--jobs J` processes J samples at once, each with
`N // J` threads (`--jobs` may not exceed `--threads`). At most
`download.max_parallel` (default 2) downloads run at the same time. Start with
`--jobs 1`; increase it on machines with many cores and enough disk space, since every
concurrent sample needs its own temporary space.

### Temporary files

`fasterq-dump` writes large temporary files. They go to `<outdir>/tmp` unless you choose
`--tmpdir`, e.g. a fast local disk. Paths with spaces are fine.

## 4. Read the results

1. `07_reports/summary_report.html`: status of every sample, failures with advice,
   screening flags, before/after tables, FastQC flags, software versions, parameters.
2. `03_raw_qc/multiqc/raw_multiqc_report.html` and
   `05_clean_qc/multiqc/clean_multiqc_report.html`: compare the plots.
3. `06_final_fastq/`: the cleaned files for alignment, with `MANIFEST.csv` and `SHA256SUMS`.

## 5. When something fails

- A failing sample is marked FAILED, its later steps SKIPPED; other samples continue
  (`--on-error stop` stops after the first failure instead).
- The exit code is 4 if some samples completed, 5 if none did.
- Fix the cause, then run the same command with `--resume`; completed, unchanged steps are
  reused and only what is needed is recomputed.

## 6. Changing parameters later

Run again with `--resume` and the new options. Steps whose inputs and settings are
unchanged are reused (downloads, conversion, raw QC); fastp and everything after it are
recomputed. Results are never a mix of old and new settings.

## 7. Configuration file

Every setting can be given in YAML (`--config my.yaml`); command-line options override
the file. See `config/default.yaml` for all keys with explanations and
`examples/config.yaml` for an example. The configuration actually used is saved as
`pipeline_info/config_used.yaml`.

## 8. Next steps (outside this tool)

The cleaned FASTQ files can be aligned (e.g. HISAT2, STAR) and quantified (e.g.
featureCounts). Use the `condition` and other metadata columns from
`07_reports/sample_manifest.csv` to set up the downstream design.
