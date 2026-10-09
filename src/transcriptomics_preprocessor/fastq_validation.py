"""Structural validation of FASTQ files (plain or gzip-compressed).

Checks performed for every record (``full`` mode) or for the first *N*
records (``quick`` mode):

* four-line record structure: ``@header``, sequence, ``+``, quality;
* sequence and quality strings have identical length;
* sequence characters are IUPAC nucleotide codes; quality characters are in
  the Phred+33 printable range (``!`` .. ``~``);
* gzip stream integrity (truncated or corrupt archives are reported);
* for paired-end files, R1 and R2 are read in lockstep: equal record counts
  and matching read identifiers (ignoring ``/1`` ``/2`` suffixes and
  comments after the first whitespace).

Counts produced in ``full`` mode are exact and are used as the reference read
counts in reports. In ``quick`` mode counts are *not* reported.
"""

from __future__ import annotations

import gzip
import io
import zlib
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import IO, Any

from .errors import IntegrityError, InterruptedRun
from .runner import STOP_EVENT

_SEQ_CHARS = frozenset(b"ACGTNacgtnRYKMSWBDHVUrykmswbdhvu.")
# bytes.translate(None, delete) removes allowed bytes in C; anything left is invalid.
_SEQ_DELETE = bytes(sorted(_SEQ_CHARS))
_QUAL_DELETE = bytes(range(33, 127))  # Phred+33 printable range '!'..'~'
_GZIP_MAGIC = b"\x1f\x8b"
_STOP_CHECK_EVERY = 200_000


@dataclass
class FastqStats:
    """Measured properties of one FASTQ file."""

    path: str
    records: int | None
    bases: int | None
    min_length: int | None
    max_length: int | None
    mean_length: float | None
    compressed: bool
    complete: bool  # False in quick mode (only a prefix was parsed)
    file_size: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _open_binary(path: Path) -> tuple[IO[bytes], bool]:
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == _GZIP_MAGIC:
        return io.BufferedReader(gzip.open(path, "rb"), buffer_size=1 << 20), True
    return open(path, "rb", buffering=1 << 20), False


def _read_id(header: bytes) -> bytes:
    rid = header[1:].split(None, 1)[0] if len(header) > 1 else b""
    if rid.endswith((b"/1", b"/2")):
        rid = rid[:-2]
    return rid


class _RecordReader:
    """Iterates over FASTQ records, raising :class:`IntegrityError` on malformed input."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.fh, self.compressed = _open_binary(path)
        self.record_no = 0

    def close(self) -> None:
        self.fh.close()

    def _fail(self, what: str) -> IntegrityError:
        return IntegrityError(
            f"{self.path.name}: record {self.record_no}: {what}",
            hint="The file is malformed or incomplete; re-run the conversion step "
            "(use --resume after deleting the file, or --overwrite).",
        )

    def __iter__(self) -> Iterator[tuple[bytes, bytes, bytes]]:
        readline = self.fh.readline
        try:
            while True:
                header = readline()
                if not header:
                    return
                self.record_no += 1
                seq = readline()
                plus = readline()
                qual = readline()
                if not qual:
                    raise self._fail("truncated record (file ends inside a record)")
                header = header.rstrip(b"\r\n")
                seq = seq.rstrip(b"\r\n")
                plus = plus.rstrip(b"\r\n")
                qual = qual.rstrip(b"\r\n")
                if not header.startswith(b"@") or len(header) < 2:
                    raise self._fail(
                        "header line does not start with '@' followed by an identifier"
                    )
                if not plus.startswith(b"+"):
                    raise self._fail("third line does not start with '+'")
                if len(seq) != len(qual):
                    raise self._fail(
                        f"sequence length ({len(seq)}) differs from quality length ({len(qual)})"
                    )
                if not seq:
                    # Zero-length reads are legal FASTQ but useless downstream; flag them.
                    raise self._fail("zero-length read")
                bad_seq = seq.translate(None, _SEQ_DELETE)
                if bad_seq:
                    bad = sorted({chr(c) for c in bad_seq})
                    raise self._fail(f"invalid sequence character(s): {''.join(bad)!r}")
                if qual.translate(None, _QUAL_DELETE):
                    raise self._fail("quality characters outside the Phred+33 range")
                yield header, seq, qual
        except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
            raise IntegrityError(
                f"{self.path.name}: compressed data is corrupt or truncated ({exc})",
                hint="The file was probably cut off during writing or copying; regenerate it.",
            ) from exc


class _Accumulator:
    def __init__(self) -> None:
        self.records = 0
        self.bases = 0
        self.min_len: int | None = None
        self.max_len: int | None = None

    def add(self, length: int) -> None:
        self.records += 1
        self.bases += length
        if self.min_len is None or length < self.min_len:
            self.min_len = length
        if self.max_len is None or length > self.max_len:
            self.max_len = length

    def stats(self, path: Path, compressed: bool, complete: bool) -> FastqStats:
        return FastqStats(
            path=str(path),
            records=self.records if complete else None,
            bases=self.bases if complete else None,
            min_length=self.min_len,
            max_length=self.max_len,
            mean_length=(self.bases / self.records) if self.records else None,
            compressed=compressed,
            complete=complete,
            file_size=path.stat().st_size,
        )


def _precheck(path: Path) -> None:
    if not path.exists():
        raise IntegrityError(f"expected FASTQ file is missing: {path}")
    if not path.is_file():
        raise IntegrityError(f"not a regular file: {path}")
    if path.stat().st_size == 0:
        raise IntegrityError(f"FASTQ file is empty (0 bytes): {path}")


def _check_stop() -> None:
    if STOP_EVENT.is_set():
        raise InterruptedRun("interrupted during FASTQ validation")


def validate_fastq(path: Path, max_records: int | None = None) -> FastqStats:
    """Validate a single-end FASTQ file. ``max_records`` enables quick mode."""
    path = Path(path)
    _precheck(path)
    reader = _RecordReader(path)
    acc = _Accumulator()
    complete = True
    try:
        for _h, seq, _q in reader:
            acc.add(len(seq))
            if acc.records % _STOP_CHECK_EVERY == 0:
                _check_stop()
            if max_records is not None and acc.records >= max_records:
                complete = False
                break
    finally:
        reader.close()
    if acc.records == 0:
        raise IntegrityError(f"{path.name}: no FASTQ records found (only empty lines?)")
    return acc.stats(path, reader.compressed, complete)


def validate_fastq_pair(
    r1: Path, r2: Path, max_records: int | None = None
) -> tuple[FastqStats, FastqStats]:
    """Validate paired files in lockstep (counts and mate identifiers must match)."""
    r1, r2 = Path(r1), Path(r2)
    _precheck(r1)
    _precheck(r2)
    if r1.resolve() == r2.resolve():
        raise IntegrityError(f"R1 and R2 refer to the same file: {r1}")
    rd1, rd2 = _RecordReader(r1), _RecordReader(r2)
    a1, a2 = _Accumulator(), _Accumulator()
    complete = True
    it1, it2 = iter(rd1), iter(rd2)
    try:
        while True:
            rec1 = next(it1, None)
            rec2 = next(it2, None)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                shorter = r1.name if rec1 is None else r2.name
                raise IntegrityError(
                    f"paired files have different numbers of reads: {shorter} ends first "
                    f"(after {min(a1.records, a2.records)} records)",
                    hint="R1/R2 must contain the same reads in the same order; "
                    "regenerate the files instead of using them.",
                )
            if _read_id(rec1[0]) != _read_id(rec2[0]):
                raise IntegrityError(
                    f"read {a1.records + 1}: mate identifiers differ between R1 and R2 "
                    f"({_read_id(rec1[0])!r} vs {_read_id(rec2[0])!r}); pairing is broken",
                )
            a1.add(len(rec1[1]))
            a2.add(len(rec2[1]))
            if a1.records % _STOP_CHECK_EVERY == 0:
                _check_stop()
            if max_records is not None and a1.records >= max_records:
                complete = False
                break
    finally:
        rd1.close()
        rd2.close()
    if a1.records == 0:
        raise IntegrityError(f"{r1.name}/{r2.name}: no FASTQ records found")
    return a1.stats(r1, rd1.compressed, complete), a2.stats(r2, rd2.compressed, complete)
