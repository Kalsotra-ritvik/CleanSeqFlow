# Changelog

All notable changes are documented here. Versions follow [Semantic Versioning](https://semver.org/).

## [0.1.0] - unreleased (release candidate)

First release candidate.

### Added
- `transcriptomics-preprocess` command with `run`, `check-input`, `doctor` and `explain`.
- SRA run accession and CSV sample-sheet input with strict validation.
- Download (`prefetch`, retries, network-error classification), conversion
  (`fasterq-dump --split-3`, layout verification, unpaired reads kept separately),
  record-level FASTQ validation, FastQC and MultiQC on raw and cleaned reads,
  fastp cleaning, before/after comparison, screening flags, HTML/JSON/CSV reports,
  final manifest with SHA-256 checksums.
- Run journal with fingerprinted checkpoints: `--resume`, `--overwrite`, `--dry-run`,
  per-sample failure isolation, `--jobs` concurrency, graceful interruption.
- `install.sh` / `uninstall.sh`: user-level installation with a SHA-256-pinned micromamba
  and an explicit, checksummed conda lock file.

### Software versions (lock file of 2026-10-09)
SRA Toolkit 3.4.1, FastQC 0.12.1, fastp 1.3.6, MultiQC 1.35, pigz 2.8, Python 3.12.15,
micromamba 2.9.0.

### Open items before a public release
- Choose a license (project owner decision) and add `LICENSE`.
- Fill in author names/ORCIDs in `CITATION.cff`.
- Run one real download from NCBI (not possible in the development sandbox).
