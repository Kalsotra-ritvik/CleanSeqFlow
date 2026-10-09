# Troubleshooting

Every error message ends with a "What to do" line. Detailed tool output is in
`<outdir>/logs/commands/<sample>/<step>.log`; the overall log is `<outdir>/logs/pipeline.log`.

## Installation

| Message | Cause and fix |
|---|---|
| `cannot reach https://conda.anaconda.org` | No internet or a proxy is needed. Set `HTTPS_PROXY=http://proxy:port` and re-run. |
| `missing basic system utilities: …` | Ask an administrator to install them (`sudo apt install curl tar bzip2 coreutils`). |
| `checksum mismatch for micromamba…` | The download was corrupted or altered; it was not used. Re-run later. If it persists, open an issue. |
| `not enough free space` | Free about 4 GB or use `--prefix` on another disk. |
| `… exists, is not empty and was not created by this installer` | Choose another `--prefix`. |
| `the smoke test failed` | See `install.log`; try `./install.sh --repair`. |
| `transcriptomics-preprocess: command not found` | `export PATH="$HOME/.local/bin:$PATH"` or log in again. |
| Launcher says the installation is incomplete | `./install.sh --repair`. |

## Input

| Message | Fix |
|---|---|
| `'GSE…' is a GEO Series, not an SRA run accession` | Look up the runs (SRR…) with the SRA Run Selector. |
| `sample name '…' is not allowed` | Use letters, digits, `.`, `_`, `-`; start with a letter or digit. |
| `appears N times` / `listed twice` | Sample names and accessions must be unique. |
| `is not UTF-8 text` | Save as "CSV UTF-8" in your spreadsheet program. |
| `jobs (J) cannot exceed threads (T)` | Lower `--jobs` or raise `--threads`. |
| `unknown configuration key` | Check spelling against `config/default.yaml`. |

## Running

| Message | Cause and fix |
|---|---|
| `already contains results from a previous run` | Add `--resume` (continue) or `--overwrite` (recompute), or choose a new `--outdir`. |
| `is not empty and was not created by this tool` | Use a new, empty output directory. |
| `another run is already using this output directory` | Wait for the other run, or use another `--outdir`. |
| `NCBI could not be reached (network problem)` | Restore connectivity, then `--resume`. Partial downloads continue. |
| `accession appears to be invalid, unavailable, too large … or access-controlled` | Check the accession page on NCBI; for very large runs raise `download.max_size`. Controlled-access (dbGaP) data are not supported. |
| `declared paired-end, but … single-end reads only` (or vice versa) | Correct `layout` for that sample. |
| `fasterq-dump ran out of disk space` | Free space or use `--tmpdir` on a larger disk, then `--resume`. |
| `paired files have different numbers of reads` / `pairing is broken` | The FASTQ files are damaged; delete them and `--resume` to regenerate. |
| `compressed data is corrupt or truncated` | A file was cut off (e.g. disk full); delete it and `--resume`. |
| `unexpected internal error` | A bug. Please open an issue with `logs/pipeline.log`. |
| `the run journal … is unreadable or corrupt` | Restore `pipeline_info/state.json.bak` as `state.json`, or restart with `--overwrite`. |

## Interpreting results

| Observation | Meaning |
|---|---|
| FastQC FAILs *Per base sequence content* on raw RNA-seq | Usually random-hexamer priming bias in the first bases; common and normally harmless. |
| High duplication | Expected in RNA-seq (abundant transcripts). Not removed by design. |
| "Residual adapter content … above the screening threshold" | Adapter auto-detection may have missed the adapter; pass `--adapter-r1`/`--adapter-r2` and `--resume`. |
| "fastp reported no adapter trimming … automatic adapter detection probably failed" | As above. |
| Few reads/pairs after cleaning flagged | The flag uses a configurable screening threshold (default 10 million). Adequate depth depends on the study. |
| `INCONSISTENT` screening status | Read counts from fastp and from the files disagree. Please report it; do not use the sample until it is understood. |
