"""Small, dependency-free helpers: atomic writes, hashing, sizes, disk space."""

from __future__ import annotations

import contextlib
import datetime as _dt
import hashlib
import json
import os
import shutil
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

_UNITS = ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]


def utc_now_iso() -> str:
    """Current UTC time as an ISO-8601 string with second precision."""
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


def human_bytes(n: int | float | None) -> str:
    """Format a byte count using binary units (``1536 -> '1.5 KiB'``)."""
    if n is None:
        return "unknown"
    value = float(n)
    for unit in _UNITS:
        if abs(value) < 1024 or unit == _UNITS[-1]:
            return f"{int(value)} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    raise AssertionError("unreachable")  # pragma: no cover


def parse_size(text: str) -> int:
    """Parse a size such as ``'500M'``, ``'20G'`` or ``'1024'`` into bytes.

    Suffixes K, M, G, T use powers of 1024. Raises ``ValueError`` on bad input.
    """
    s = str(text).strip().upper().removesuffix("B").removesuffix("I")
    if not s:
        raise ValueError("empty size")
    factors = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
    factor = 1
    if s[-1] in factors:
        factor = factors[s[-1]]
        s = s[:-1]
    value = float(s)
    if value < 0 or value != value:  # negative or NaN
        raise ValueError(f"invalid size: {text!r}")
    return int(value * factor)


def free_bytes(path: Path) -> int | None:
    """Free bytes on the filesystem holding *path* (or its nearest existing parent)."""
    p = Path(path).absolute()
    while not p.exists():
        if p.parent == p:
            return None
        p = p.parent
    try:
        return shutil.disk_usage(p).free
    except OSError:
        return None


def atomic_write_text(path: Path, text: str, encoding: str = "utf-8") -> None:
    """Write *text* to *path* atomically (temporary file + ``os.replace``)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding=encoding, newline="") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(tmp)
        raise


def atomic_write_json(path: Path, data: Any) -> None:
    """Serialise *data* as indented JSON and write it atomically."""
    atomic_write_text(path, json.dumps(data, indent=2, sort_keys=False, default=str) + "\n")


def sha256_file(path: Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    """SHA-256 hex digest of a file, read in chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def stable_hash(data: Any) -> str:
    """Deterministic SHA-256 of a JSON-serialisable structure."""
    blob = json.dumps(data, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def file_signature(path: Path) -> dict[str, Any]:
    """Cheap identity of a file: size and modification time (nanoseconds)."""
    st = Path(path).stat()
    return {"size": st.st_size, "mtime_ns": st.st_mtime_ns}


def remove_if_exists(path: Path) -> None:
    """Remove a single file if present. Never removes directories."""
    p = Path(path)
    if p.is_dir() and not p.is_symlink():
        raise IsADirectoryError(f"refusing to remove directory with remove_if_exists: {p}")
    with contextlib.suppress(FileNotFoundError):
        p.unlink()


def is_within(path: Path, root: Path) -> bool:
    """True if *path* resolves to a location inside *root*."""
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
    except ValueError:
        return False
    return True


def safe_rmtree(path: Path, root: Path) -> None:
    """Recursively delete *path* only if it lies strictly inside *root*.

    Guards against deleting the root itself, anything outside it, and
    following symlinks out of it.
    """
    p = Path(path)
    if p.is_symlink():
        p.unlink()
        return
    if not p.exists():
        return
    rp, rr = p.resolve(), Path(root).resolve()
    if rp == rr or not is_within(rp, rr):
        raise PermissionError(f"refusing to delete {rp}: not strictly inside {rr}")
    shutil.rmtree(rp)


def check_writable_dir(path: Path) -> str | None:
    """Return an error message if *path* cannot be created/written, else ``None``."""
    p = Path(path)
    try:
        p.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        return f"cannot create directory {p}: {exc.strerror or exc}"
    if not p.is_dir():
        return f"{p} exists and is not a directory"
    try:
        with tempfile.NamedTemporaryFile(dir=p, prefix=".write-test-", delete=True):
            pass
    except OSError as exc:
        return f"directory {p} is not writable: {exc.strerror or exc}"
    return None


@contextlib.contextmanager
def working_directory(path: Path) -> Iterator[None]:
    """Temporarily change the current working directory."""
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)
