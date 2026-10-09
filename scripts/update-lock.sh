#!/usr/bin/env bash
# Regenerate environment-lock-linux-64.txt from environment.yml.
# Usage: scripts/update-lock.sh [path/to/micromamba]
# Afterwards: reinstall (./install.sh), run the full test suite and update the
# version tables in README.md and CHANGELOG.md.
set -Eeuo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
MM="${1:-${XDG_DATA_HOME:-$HOME/.local/share}/transcriptomics-preprocessor/bin/micromamba}"
[ -x "$MM" ] || { echo "micromamba not found at $MM (run ./install.sh first or pass its path)" >&2; exit 1; }
WORK="$(mktemp -d)"
trap 'rm -rf -- "$WORK"' EXIT
export MAMBA_ROOT_PREFIX="$WORK/root"
"$MM" create -y --no-rc -p "$WORK/env" -f "$REPO/environment.yml" --strict-channel-priority
{
  echo "# Explicit, checksummed lock of the runtime environment for linux-64."
  echo "# Solved from environment.yml with micromamba $("$MM" --version) on $(date -u +%Y-%m-%d) (channels: conda-forge, bioconda;"
  echo "# strict channel priority; no Anaconda 'defaults' channel). install.sh uses this file by default."
  echo "# Regenerate with: scripts/update-lock.sh (then run the full test suite)."
  "$MM" env export -p "$WORK/env" --explicit --md5 | grep -v '^# This file may\|^# \$ conda create'
} > "$REPO/environment-lock-linux-64.txt.new"
mv "$REPO/environment-lock-linux-64.txt.new" "$REPO/environment-lock-linux-64.txt"
echo "Updated $REPO/environment-lock-linux-64.txt"
