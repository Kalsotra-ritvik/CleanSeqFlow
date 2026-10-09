from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from transcriptomics_preprocessor.errors import IntegrityError
from transcriptomics_preprocessor.fastq_validation import validate_fastq, validate_fastq_pair
from transcriptomics_preprocessor.synthetic import write_synthetic_fastq

REC = "@r{i} 1 length=8\nACGTNACG\n+\nIIIIIIII\n"


def _write(path: Path, text: str) -> Path:
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


@pytest.mark.parametrize("name", ["x.fastq", "x.fastq.gz"])
def test_counts_are_exact(tmp_path: Path, name: str) -> None:
    p = _write(tmp_path / name, "".join(REC.format(i=i) for i in range(10)))
    st = validate_fastq(p)
    assert (st.records, st.bases, st.min_length, st.max_length) == (10, 80, 8, 8)
    assert st.compressed == name.endswith(".gz") and st.complete


def test_synthetic_paired_counts(tmp_path: Path) -> None:
    r1, r2 = tmp_path / "a_R1.fastq.gz", tmp_path / "a_R2.fastq.gz"
    write_synthetic_fastq(r1, r2, n_reads=321)
    s1, s2 = validate_fastq_pair(r1, r2)
    assert s1.records == s2.records == 321
    assert s1.bases == 321 * 100


def test_crlf_line_endings_accepted(tmp_path: Path) -> None:
    p = _write(tmp_path / "x.fastq", REC.format(i=1).replace("\n", "\r\n"))
    assert validate_fastq(p).records == 1


@pytest.mark.parametrize(
    ("text", "fragment"),
    [
        ("r1\nACGT\n+\nIIII\n", "'@'"),
        ("@r1\nACGT\n-\nIIII\n", "'+'"),
        ("@r1\nACGT\n+\nIII\n", "differs from quality length"),
        ("@r1\nACGX\n+\nIIII\n", "invalid sequence character"),
        ("@r1\nACGT\n+\nII I\n", "Phred+33"),
        ("@r1\nACGT\n+\n", "truncated record"),
        ("@r1\n\n+\n\n", "zero-length"),
        ("@\nACGT\n+\nIIII\n", "'@'"),
        ("\n\n\n", "truncated record"),  # 3 lines: an incomplete record
        ("\n\n\n\n", "'@'"),  # 4 blank lines: a record without a header
    ],
)
def test_structural_errors(tmp_path: Path, text: str, fragment: str) -> None:
    p = _write(tmp_path / "bad.fastq", text)
    with pytest.raises(IntegrityError) as e:
        validate_fastq(p)
    assert fragment in e.value.message


def test_empty_and_missing_files(tmp_path: Path) -> None:
    empty = tmp_path / "e.fastq.gz"
    empty.write_bytes(b"")
    with pytest.raises(IntegrityError, match="empty"):
        validate_fastq(empty)
    with pytest.raises(IntegrityError, match="missing"):
        validate_fastq(tmp_path / "nope.fastq")
    gz_empty = _write(tmp_path / "e2.fastq.gz", "")
    with pytest.raises(IntegrityError, match="no FASTQ records"):
        validate_fastq(gz_empty)


def test_truncated_gzip_detected(tmp_path: Path) -> None:
    p = tmp_path / "t.fastq.gz"
    write_synthetic_fastq(p, None, n_reads=2000)
    data = p.read_bytes()
    p.write_bytes(data[: len(data) // 2])
    with pytest.raises(IntegrityError, match="corrupt or truncated"):
        validate_fastq(p)


def test_corrupt_gzip_payload_detected(tmp_path: Path) -> None:
    p = tmp_path / "c.fastq.gz"
    write_synthetic_fastq(p, None, n_reads=500)
    data = bytearray(p.read_bytes())
    for i in range(200, 260):
        data[i] ^= 0xFF
    p.write_bytes(bytes(data))
    with pytest.raises(IntegrityError):
        validate_fastq(p)


def test_pair_count_mismatch(tmp_path: Path) -> None:
    r1 = _write(tmp_path / "r1.fastq", "".join(REC.format(i=i) for i in range(5)))
    r2 = _write(tmp_path / "r2.fastq", "".join(REC.format(i=i) for i in range(4)))
    with pytest.raises(IntegrityError, match="different numbers of reads"):
        validate_fastq_pair(r1, r2)


def test_pair_name_mismatch(tmp_path: Path) -> None:
    r1 = _write(tmp_path / "r1.fastq", REC.format(i=1) + REC.format(i=2))
    r2 = _write(tmp_path / "r2.fastq", REC.format(i=1) + REC.format(i=3))
    with pytest.raises(IntegrityError, match="pairing is broken"):
        validate_fastq_pair(r1, r2)


def test_pair_suffixes_ignored(tmp_path: Path) -> None:
    r1 = _write(tmp_path / "r1.fastq", "@x/1\nACGT\n+\nIIII\n")
    r2 = _write(tmp_path / "r2.fastq", "@x/2 extra\nACGT\n+\nIIII\n")
    s1, _ = validate_fastq_pair(r1, r2)
    assert s1.records == 1


def test_same_file_as_both_mates(tmp_path: Path) -> None:
    r1 = _write(tmp_path / "r1.fastq", REC.format(i=1))
    with pytest.raises(IntegrityError, match="same file"):
        validate_fastq_pair(r1, r1)


def test_quick_mode_does_not_report_counts(tmp_path: Path) -> None:
    p = _write(tmp_path / "x.fastq", "".join(REC.format(i=i) for i in range(50)))
    st = validate_fastq(p, max_records=10)
    assert st.records is None and st.bases is None and not st.complete


def test_interrupt_is_honoured(tmp_path: Path) -> None:
    from transcriptomics_preprocessor.errors import InterruptedRun
    from transcriptomics_preprocessor.runner import STOP_EVENT

    p = tmp_path / "big.fastq"
    p.write_text(REC.format(i=0) * 250_000)
    STOP_EVENT.set()
    with pytest.raises(InterruptedRun):
        validate_fastq(p)
