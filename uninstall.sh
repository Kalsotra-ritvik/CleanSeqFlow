#!/usr/bin/env bash
# uninstall.sh — remove an installation created by install.sh.
#
# Safety rules:
#   * only a prefix containing the installer's marker file is removed;
#   * only launchers containing the installer's marker line are removed;
#   * pipeline OUTPUT directories (your results) are never touched;
#   * nothing is deleted without confirmation (or --yes).
set -Eeuo pipefail

readonly APP="transcriptomics-preprocessor"
readonly CMD="transcriptomics-preprocess"
readonly ROOT_MARK=".tpp-install-root"

PREFIX="${XDG_DATA_HOME:-$HOME/.local/share}/${APP}"
BIN_DIR="$HOME/.local/bin"
YES=0
DRY_RUN=0
PRUNE_ONLY=0

usage() {
    cat <<EOF
Usage: ./uninstall.sh [options]

Removes the ${APP} software installation (environment, private micromamba,
package cache and launchers). Your analysis results are NOT touched.

Options:
  --prefix DIR        installation directory (default: ${PREFIX})
  --bin-dir DIR       launcher directory (default: ${BIN_DIR})
  --prune-old-envs    only remove older, inactive environments; keep the installation
  --dry-run           show what would be removed, remove nothing
  -y, --yes           do not ask for confirmation
  -h, --help          show this help
EOF
}

die() { printf 'ERROR: %s\n' "$1" >&2; [ $# -gt 1 ] && printf '  What to do: %s\n' "$2" >&2; exit 1; }

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) [ $# -ge 2 ] || die "--prefix needs a directory"; PREFIX="$2"; shift 2 ;;
        --bin-dir) [ $# -ge 2 ] || die "--bin-dir needs a directory"; BIN_DIR="$2"; shift 2 ;;
        --prune-old-envs) PRUNE_ONLY=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -y|--yes) YES=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; die "unknown option: $1" ;;
    esac
done

while [ "${#PREFIX}" -gt 1 ] && [ "${PREFIX%/}" != "$PREFIX" ]; do PREFIX="${PREFIX%/}"; done
case "$PREFIX" in ""|"/") die "refusing to remove '${PREFIX:-<empty>}'" ;; esac
case "$PREFIX" in /*) ;; *) PREFIX="$(pwd -P)/$PREFIX" ;; esac
case "$PREFIX" in
    "$HOME"|"/usr"|"/usr/local"|"/opt"|"/home"|"/root"|"/etc")
        die "refusing to remove '$PREFIX'" ;;
esac
[ -d "$PREFIX" ] || die "no installation found at $PREFIX" "Pass --prefix if you installed elsewhere."
[ -f "$PREFIX/$ROOT_MARK" ] || die "$PREFIX does not contain the installer marker ($ROOT_MARK); refusing to delete it" \
    "This safety check prevents deleting directories the installer did not create."

confirm() {
    [ "$DRY_RUN" = 1 ] && { echo "Dry run: nothing was removed."; exit 0; }
    [ "$YES" = 1 ] && return 0
    if [ ! -t 0 ]; then die "confirmation required" "Re-run with --yes to remove without asking."; fi
    printf 'Proceed? [y/N] '
    read -r answer
    case "$answer" in y|Y|yes|YES) ;; *) echo "Nothing removed."; exit 0 ;; esac
}

if [ "$PRUNE_ONLY" = 1 ]; then
    active="$(readlink "$PREFIX/env" 2>/dev/null || true)"
    targets=()
    for d in "$PREFIX"/envs/tpp-*; do
        [ -d "$d" ] || continue
        [ "envs/$(basename "$d")" = "$active" ] && continue
        [ -f "$d.managed" ] || continue
        targets+=("$d")
    done
    if [ ${#targets[@]} -eq 0 ]; then echo "No inactive environments to remove."; exit 0; fi
    echo "Inactive environments to remove (active: ${active:-none}):"
    printf '  %s\n' "${targets[@]}"
    confirm
    for d in "${targets[@]}"; do rm -rf -- "$d" "$d.managed"; done
    echo "Removed ${#targets[@]} inactive environment(s)."
    exit 0
fi

launchers=()
for f in "$BIN_DIR/$CMD" "$PREFIX/bin/$CMD"; do
    if [ -f "$f" ] && grep -q "TPP-LAUNCHER" "$f" 2>/dev/null; then
        launchers+=("$f")
    elif [ -e "$f" ]; then
        echo "Leaving $f in place: it was not created by install.sh."
    fi
done

echo "The following will be removed:"
echo "  $PREFIX  ($(du -sh "$PREFIX" 2>/dev/null | cut -f1) - environment, micromamba, package cache, install log)"
for f in "${launchers[@]}"; do echo "  $f  (launcher)"; done
echo "Your pipeline output directories are NOT affected."
confirm
for f in "${launchers[@]}"; do rm -f -- "$f"; done
rm -rf -- "$PREFIX"
echo "Uninstalled. The repository directory $(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P) was left in place; delete it yourself if you no longer need it."
