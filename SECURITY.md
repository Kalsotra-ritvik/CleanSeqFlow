# Security policy

## Reporting a vulnerability

Please do **not** open a public issue for security problems. Use GitHub's private
vulnerability reporting ("Security" tab → "Report a vulnerability") for this repository.
Include a description, steps to reproduce and the affected version. You should receive a
response within 14 days.

## Security-relevant design

- `install.sh` never uses `sudo`, never pipes downloads into a shell, and verifies the
  micromamba binary against a pinned SHA-256 before executing it. Conda packages come from
  an explicit lock file with checksums.
- External programs are executed with argument lists (no shell), so sample names, paths
  and accessions cannot inject commands. Sample names are restricted to a safe character set.
- Deletion is limited to files and directories the tool created inside its own output or
  installation directory, with guards against symlinks and paths outside those directories.
- MultiQC runs with its version check, AI features and MegaQC upload disabled. The only
  network access during a run is the SRA download by `prefetch`.
- Logs contain commands and paths but no environment variables or credentials.

## Supported versions

Only the latest release receives fixes.

## Known, accepted risks

- Extra sample-sheet columns (e.g. `condition`) are copied unchanged into the output CSV
  files. A value beginning with `=`, `+`, `-` or `@` could be interpreted as a formula if a
  CSV is opened in a spreadsheet program. Values are not altered because they are the
  user's own metadata; do not open CSVs from untrusted sources with formula evaluation on.
- `install.sh` trusts `conda.anaconda.org` over HTTPS for package downloads; package
  integrity is checked against the MD5 checksums in the lock file (the format conda uses),
  and the bootstrap binary against a pinned SHA-256.
