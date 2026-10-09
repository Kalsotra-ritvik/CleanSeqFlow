Real outputs captured from FastQC 0.12.1 and fastp 1.3.6 runs on the synthetic
reads produced by `transcriptomics_preprocessor.synthetic` (2026-10-09).
The fastp JSON files were trimmed to remove bulky per-position curves; every
summary, filtering_result, adapter_cutting, duplication and command section is
unmodified. They are used to test the metric parsers against real formats.
The absolute paths inside the fastp `command` field were replaced by `/path/to/results`.
