# Contributing

Thank you for helping. Bug reports, documentation fixes and code are all welcome.

## Reporting a problem

Open an issue using the bug-report template and include: the exact command, the
`transcriptomics-preprocess doctor` output, the error message, and the relevant files from
`<outdir>/logs/`. Please do not attach sequencing data or anything confidential.

## Development setup

```bash
./install.sh                                   # creates the runtime environment
P=~/.local/share/transcriptomics-preprocessor
eval "$("$P/bin/micromamba" shell activate -s bash -p "$P/env")"
python -m pip install --no-deps --no-build-isolation -e .   # editable install for development
python -m pip install ruff mypy types-PyYAML                 # developer tools (not runtime deps)
```

## Before opening a pull request

```bash
ruff format --check src tests
ruff check src tests
mypy                      # strict mode, configured in pyproject.toml
python -m pytest tests    # tool-based tests run only inside the environment
shellcheck install.sh uninstall.sh scripts/*.sh
```

Guidelines:

- Keep modules small and testable; external programs are always called with argument
  lists, never through a shell.
- Every bug fix comes with a test that fails without it.
- Do not weaken a test or a scientific check to make a test pass; explain any change to
  an expected value in the pull request.
- Never invent metric values. Missing data stay missing.
- Keep user-facing wording consistent with `explain.py` and the documentation.
- Tests must not download from NCBI; use the doubles in `tests/fixtures/fake_bin`.
- Changing tool versions means regenerating the lock file (`scripts/update-lock.sh`) and
  re-running the full suite.
