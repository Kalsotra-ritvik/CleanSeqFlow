# Installation

## Summary

```bash
git clone https://github.com/Kalsotra-ritvik/CleanSeqFlow.git
cd CleanSeqFlow
./install.sh
transcriptomics-preprocess doctor
```

`install.sh` is the single, authoritative way to install the software. It needs no
administrator rights and changes nothing outside two locations:

| Location (default) | Contents |
|---|---|
| `~/.local/share/transcriptomics-preprocessor/` | private micromamba, package cache, environment(s), install log |
| `~/.local/bin/transcriptomics-preprocess` | launcher script |

It never runs `sudo`, never modifies the system Python or your own Python packages, never
edits `.bashrc`/`.profile`, and never pipes a download into a shell.

## What the installer does, step by step

1. **Platform check**: Linux required; x86-64 (tested) or aarch64 (untested). Non-Ubuntu
   distributions get a warning but are not blocked.
2. **Utilities and space**: needs `tar`, `bzip2`, `sha256sum`, `curl` or `wget`, and about
   4 GB free. If something is missing it tells you what to ask an administrator to install
   (e.g. `sudo apt install curl bzip2`). It does not run `sudo` itself.
3. **micromamba**: reuses the private copy from an earlier installation, or uses the one
   you pass with `--micromamba PATH`. Otherwise it downloads the conda-forge package
   `micromamba-2.9.0-0.tar.bz2` from `https://conda.anaconda.org` and compares its SHA-256
   with the value hard-coded in `install.sh` (taken from conda-forge's repository index).
   On a mismatch the file is discarded and never executed.
4. **Environment**: creates `envs/tpp-<id>` from `environment-lock-linux-64.txt`, an explicit
   list of packages with MD5 checksums, which micromamba verifies. `<id>` is derived from the
   lock file's content, so a changed lock file produces a new environment next to the old
   one. User configuration files such as `~/.condarc` are ignored (`--no-rc`).
5. **Package**: installs this repository's Python package into the environment with
   `pip --no-deps --no-index` (no access to PyPI; all dependencies come from the lock).
6. **Verification**: runs `transcriptomics-preprocess doctor` inside the environment and
   prints the detected versions.
7. **Smoke test**: FastQC, fastp and MultiQC process a small synthetic dataset, and the
   results are checked (e.g. FastQC must count exactly 2,000 reads; fastp must trim
   adapters that were deliberately inserted). The SRA Toolkit executables are started, but
   nothing is downloaded from NCBI.
8. **Activation**: only after steps 6–7 pass is the environment marked complete and the
   `env` symlink switched to it (atomically). A failed update therefore leaves the
   previous working installation untouched. The launcher is written and tested from an
   empty environment (`env -i`).
9. **Summary**: versions, commands, and PATH advice if `~/.local/bin` is not on `PATH`.

Everything is logged to `~/.local/share/transcriptomics-preprocessor/install.log`.

## Re-running, updating, repairing

| Situation | What happens |
|---|---|
| Re-run with unchanged lock file | Existing complete environment is reused and re-verified; the package is reinstalled from the current source (about 15 s). |
| `git pull` brought a new lock file | A new environment is built and verified; then it becomes active. The old one stays until `./uninstall.sh --prune-old-envs`. |
| Interrupted installation | The incomplete environment has no completion marker; the next run deletes and rebuilds it. |
| `./install.sh --repair` | Rebuilds the environment even if it looks complete. |

## Options

```
--prefix DIR        installation directory
--bin-dir DIR       launcher directory (default ~/.local/bin)
--no-launcher       only create <prefix>/bin/transcriptomics-preprocess
--force-launcher    replace a foreign file named transcriptomics-preprocess in --bin-dir
--micromamba PATH   use an existing micromamba
--from-spec         solve from environment.yml instead of the lock file
--repair            rebuild the environment
--skip-smoke-test   skip the functional test (not recommended)
```

The installer refuses to use `/`, your home directory or system directories as the prefix,
and refuses a non-empty prefix it did not create.

## Without internet access

Installation needs `https://conda.anaconda.org` (conda-forge and bioconda). If it cannot be
reached, the installer stops with a message listing what it needs to download, and reports
that the installation is **not** complete. Behind a proxy, set `HTTPS_PROXY` before running
it. Running the pipeline additionally needs access to NCBI (`*.ncbi.nlm.nih.gov` and its
download hosts).

## Uninstalling

```bash
./uninstall.sh            # asks for confirmation
./uninstall.sh --dry-run  # show what would be removed
./uninstall.sh --yes      # no question (for scripts)
```

It removes only a prefix containing the installer's marker file and only launchers that
contain the installer's marker line. Pipeline output directories are never touched, and
neither is the cloned repository.

## Using the environment directly (developers)

```bash
P=~/.local/share/transcriptomics-preprocessor
eval "$("$P/bin/micromamba" shell activate -s bash -p "$P/env")"
which fastp fastqc multiqc prefetch
```

## Regenerating the lock file (maintainers)

Edit `environment.yml`, run `scripts/update-lock.sh`, reinstall, run the full test suite,
and update the version tables in `README.md` and `CHANGELOG.md`.
