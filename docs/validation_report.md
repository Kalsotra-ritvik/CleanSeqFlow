# Validation report: release candidate 0.1.0

Date: 2026-10-09. This report records what was tested, how, and with what result. It also
lists what was **not** tested. Test counts come from actual pytest output.

## Test environment

| Item | Value |
|---|---|
| OS | Ubuntu 24.04.5 LTS, x86-64, Linux 6.18 (container), 4 CPUs, 15 GB RAM |
| Runtime environment | created by `./install.sh` from `environment-lock-linux-64.txt` in an empty `HOME` |
| Tools | SRA Toolkit 3.4.1, FastQC 0.12.1, fastp 1.3.6, MultiQC 1.35, pigz 2.8, Python 3.12.15 |
| Network | conda.anaconda.org and PyPI reachable; **all NCBI hosts blocked** (firewall) |
| User | root (affects permission tests, see below) |

## Results

| Suite | Command | Result |
|---|---|---|
| Full suite in the installed environment | `python -m pytest tests -q -rs` (inside `~/.local/share/transcriptomics-preprocessor/env`) | **334 passed, 1 skipped** (5 min 7 s) |
| Unit + CLI tests on Python 3.10 (no tools; CI minimum) | `PYTHONPATH=src python3.10 -m pytest tests` | 307 passed, 28 skipped (tool tests) |
| Lint / format | `ruff check src tests`, `ruff format --check src tests` | passed |
| Types | `mypy` (strict, 26 source files) | no issues |
| Shell scripts | `shellcheck install.sh uninstall.sh scripts/*.sh` | no findings |
| Clean installation | `env -i HOME=<empty dir> ./install.sh` | passed (45 s; 178 packages) |
| Re-installation | same command again | passed, environment reused (14 s) |

The one skipped test (`test_unwritable_output_location`) cannot run as root, because root
ignores file permissions. It runs in GitHub Actions (non-root runner).

## Feature-to-test matrix

Status values: PASSED (assertions passed), BLOCKED (cannot be tested here; reason given),
NOT_TESTED. "fake SRA" means `prefetch`/`fasterq-dump` were test doubles; FastQC, fastp,
MultiQC and pigz were always the real programs.

| # | Feature | Evidence (tests / manual checks) | Status |
|---|---|---|---|
| 1 | Automatic installation | clean `HOME` install; idempotent rerun; launcher tested with `env -i` | PASSED |
| 2 | Dependency discovery, diagnostics | `doctor`, `doctor --json`, `--smoke-test`; missing-tool tests (exit 3) | PASSED |
| 3 | CLI help, invalid commands | `tests/test_cli.py` (usage errors exit 2, help content) | PASSED |
| 4 | Accession input | `test_accessions_*`, `test_accessions_input` | PASSED |
| 5 | Metadata CSV input | `tests/unit/test_validation.py` (delimiters, BOM, comments, extra columns) | PASSED |
| 6 | Duplicate / malformed metadata | `test_all_problems_reported_together`, `test_malformed_csv`, case clash | PASSED |
| 7 | Paired-end processing | `test_full_run_paired_and_single`: R1/R2 counts equal, names checked | PASSED (fake SRA) |
| 8 | Single-end processing | same test; `test_explicit_adapter_fixes_single_end_detection` | PASSED (fake SRA) |
| 9 | SRA download | retry, permanent-error and missing-output tests (fake); real `prefetch` against blocked NCBI reported as network failure | PASSED for failure paths; real download **BLOCKED** (no NCBI access) |
| 10 | SRA-to-FASTQ conversion | layout mismatch, unpaired reads, truncation, no-space, failure tests (fake) | PASSED (fake SRA); real `fasterq-dump` on a real archive **BLOCKED** |
| 11 | FASTQ integrity validation | 23 unit tests incl. truncated/corrupt gzip, pairing, characters | PASSED |
| 12 | Raw FastQC reports | integration tests; parser tested on real FastQC 0.12.1 output | PASSED |
| 13 | Raw MultiQC | input list contains only raw FastQC zips | PASSED |
| 14 | fastp parameter handling | 96 exhaustive switch combinations × PE/SE; boundary tests; real run with `--min-length 40` checked in the recorded command | PASSED |
| 15 | fastp HTML/JSON reports | present, parsed, cross-checked with independent counts | PASSED |
| 16 | Post-cleaning FastQC | integration tests | PASSED |
| 17 | Post-cleaning MultiQC | input list contains cleaned FastQC + fastp JSON only | PASSED |
| 18 | Reads vs read pairs | unit tests on real fastp JSON; integration counts via `zcat` | PASSED |
| 19 | Before/after comparison | values equal their sources (fastp JSON, FastQC zip, independent counts); missing values stay missing; percentage-point wording | PASSED |
| 20 | Final output manifest | SHA-256 recomputed and compared; read counts recounted | PASSED |
| 21 | Log files | per-command logs, `pipeline.log`, `events.jsonl` (event counts checked) | PASSED |
| 22 | Progress and explanations | stdout content checked; `--quiet`/`--verbose` tests | PASSED |
| 23 | HTML/JSON/CSV summaries | content, HTML escaping, failure display, rendered and inspected visually | PASSED |
| 24 | Resume | identical resume (all 18 steps reused, no download) | PASSED |
| 25 | Checkpoint integrity | changed parameter, deleted output, hand-modified output, corrupt journal | PASSED |
| 26 | Dry run | creates nothing; shows commands; reuse plan with `--resume`; journal unchanged | PASSED |
| 27 | Overwrite protection | existing run refused without flag; foreign directory refused even with `--overwrite`/`--resume` | PASSED |
| 28 | Safe interruption | SIGINT (1 job) and SIGTERM (2 parallel jobs): exit 130 quickly, no RUNNING records, resume works; grandchild processes terminated | PASSED |
| 29 | Partial failures | one sample fails (conversion, fastp, FastQC); others complete; exit 4; MultiQC failure after cleaning | PASSED |
| 30 | Temporary directory | custom `--tmpdir` with spaces; emptied after success | PASSED |
| 31 | Disk-space and permission failures | no-space message classified with advice (fake); disk warnings unit-tested indirectly; unwritable directory test | PARTIAL: permission test skipped as root; real disk exhaustion not simulated |
| 32 | Reproducibility / version recording | `software_versions.json`, `config_used.yaml`, versions in journal and report | PASSED |
| 33 | Clean-environment installation | empty `HOME`, `env -i`, no preinstalled tools | PASSED (same machine, fresh user state; not a fresh VM) |
| 34 | Uninstallation safeguards | marker checks, foreign launcher kept, no-tty refusal, dry-run, dangerous prefixes | PASSED |
| 35 | GitHub CI | workflow written and validated locally step by step | NOT_TESTED on GitHub (push to the repository was refused: no GitHub access) |

## Scientific checks

- Read counts reported for each sample equal (a) fastp's `total_reads`, and (b) an
  independent `zcat | wc -l` count of the files, before and after cleaning.
- fastp's filtering categories sum to its input count (checked automatically per sample).
- Adapter content values equal the maximum of the FastQC *Adapter Content* table read
  directly from the zip archives.
- Synthetic data with known properties: fastp's reported low-quality reads (532) equal
  exactly twice the number of low-quality pairs the generator injected (266), and
  adapter-trimmed reads equal twice the short-insert pairs (1,030 → 2,060).
- Finding: fastp's single-end adapter auto-detection found no adapter in a 1,500-read file
  (adapter content 19.40% → 19.39%), but did in a 20,000-read file. The report flags this
  and does not claim an improvement; the recommended fix (explicit adapter) was verified.
- Duplicates are never removed (no `--dedup` in any generated command; tested).

## Defects found and fixed during validation

| Defect | Found by | Fix |
|---|---|---|
| Unreachable NCBI would have been reported as "invalid accession" (real sra-tools 3.4.1 message contains "cannot resolve") | running real `prefetch` offline | network markers checked first; tests use the real message |
| Duplicate sample names missed when the duplicate row had another error | CLI test | duplicates checked over all rows |
| `--dry-run --resume` rejected although supported | integration test | options made combinable |
| `install.sh --prefix /` resolved to the current directory | manual adversarial test | root/empty prefixes rejected before normalisation |
| "increased by 6.71%" for a difference of percentages | report review | reported as percentage points |
| Negligible adapter change described as "decreased" | real run review | "essentially unchanged" below 0.5 points; failed auto-detection explained |
| Machine-specific paths in fixtures | repository audit | sanitised |

## Not tested / limitations

- **Real SRA download and conversion** (network blocked). The code paths for success were
  tested with test doubles that follow `fasterq-dump --split-3` naming; one real run on
  an NCBI accession is required before release.
- Very large datasets (performance measured only on 1 million pairs: the FASTQ validator
  processes about 160,000 read pairs per second, roughly 3 minutes per pass for 30 million pairs).
- Real disk exhaustion, SRA-lite files, runs with technical reads, aarch64, distributions
  other than Ubuntu 24.04, and a fresh virtual machine.
- Not every combination of all parameters was tested: fastp enable switches were tested
  exhaustively (96 combinations); numeric parameters at boundaries; other options
  individually.
