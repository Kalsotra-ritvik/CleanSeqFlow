"""Deterministic synthetic RNA-seq-like FASTQ generator.

Used by ``doctor --smoke-test`` and by the test suite. The reads are *not*
biologically meaningful; they are designed so that preprocessing has
something measurable to do:

* a small "transcriptome" in which a few transcripts are highly expressed
  (producing duplicated sequences, as in real RNA-seq);
* a fraction of short fragments, so reads run into the Illumina TruSeq
  adapter (adapter read-through);
* a fraction of reads with uniformly poor base qualities;
* occasional ambiguous ``N`` bases;
* quality scores that decline towards the 3' end.

Read names follow the SRA ``fasterq-dump`` default definition line
(``@<accession>.<spot> <spot> length=<len>``) with identical names for mates.
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass
from pathlib import Path
from typing import IO

#: Illumina TruSeq adapter sequences as seen at the 3' end of R1 and R2.
TRUSEQ_R1_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRUSEQ_R2_ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class SyntheticSummary:
    """What was written; useful for test assertions."""

    reads_per_file: int
    read_length: int
    paired: bool
    adapter_reads: int
    low_quality_reads: int


def _open(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "wt", compresslevel=6)
    return open(path, "w")


def _quality(rng: random.Random, length: int, poor: bool) -> str:
    chars = []
    for i in range(length):
        if poor:
            q = rng.randint(2, 14)
        else:
            # Mean quality ~37 at the 5' end, falling to ~25 at the 3' end.
            mean = 37 - 12 * (i / max(1, length - 1)) ** 2
            q = round(rng.gauss(mean, 2.5))
            q = max(2, min(41, q))
        chars.append(chr(q + 33))
    return "".join(chars)


def _maybe_add_n(rng: random.Random, seq: str, qual: str) -> tuple[str, str]:
    if rng.random() < 0.03:
        pos = rng.randrange(len(seq))
        seq = seq[:pos] + "N" + seq[pos + 1 :]
        qual = qual[:pos] + "#" + qual[pos + 1 :]
    return seq, qual


def write_synthetic_fastq(
    r1_path: Path,
    r2_path: Path | None = None,
    *,
    n_reads: int = 2000,
    read_length: int = 100,
    accession: str = "SRR0000001",
    seed: int = 1,
    adapter_fraction: float = 0.2,
    low_quality_fraction: float = 0.05,
) -> SyntheticSummary:
    """Write single-end (``r2_path is None``) or paired-end synthetic FASTQ.

    Output is gzip-compressed when the file name ends in ``.gz``.
    """
    if n_reads < 1 or read_length < 30:
        raise ValueError("n_reads must be >= 1 and read_length >= 30")
    rng = random.Random(seed)
    transcripts = [
        "".join(rng.choice("ACGT") for _ in range(rng.randint(600, 1500))) for _ in range(40)
    ]
    # A few highly expressed transcripts dominate, as in real libraries.
    weights = [50.0 if i < 3 else 1.0 for i in range(len(transcripts))]

    fh1 = _open(Path(r1_path))
    fh2 = _open(Path(r2_path)) if r2_path is not None else None
    adapter_reads = 0
    low_quality = 0
    try:
        for spot in range(1, n_reads + 1):
            tx = rng.choices(transcripts, weights=weights, k=1)[0]
            if rng.random() < adapter_fraction:
                frag_len = rng.randint(35, read_length - 15)
                adapter_reads += 1
            else:
                frag_len = rng.randint(read_length + 20, min(len(tx), read_length + 300))
            start = rng.randint(0, len(tx) - frag_len)
            fragment = tx[start : start + frag_len]
            if rng.random() < 0.5:
                fragment = reverse_complement(fragment)
            poor = rng.random() < low_quality_fraction
            low_quality += int(poor)

            r1 = (fragment + TRUSEQ_R1_ADAPTER + "A" * read_length)[:read_length]
            q1 = _quality(rng, read_length, poor)
            r1, q1 = _maybe_add_n(rng, r1, q1)
            name = f"{accession}.{spot} {spot} length={read_length}"
            fh1.write(f"@{name}\n{r1}\n+{name}\n{q1}\n")
            if fh2 is not None:
                r2 = (reverse_complement(fragment) + TRUSEQ_R2_ADAPTER + "A" * read_length)[
                    :read_length
                ]
                q2 = _quality(rng, read_length, poor)
                r2, q2 = _maybe_add_n(rng, r2, q2)
                fh2.write(f"@{name}\n{r2}\n+{name}\n{q2}\n")
    finally:
        fh1.close()
        if fh2 is not None:
            fh2.close()
    return SyntheticSummary(
        reads_per_file=n_reads,
        read_length=read_length,
        paired=r2_path is not None,
        adapter_reads=adapter_reads,
        low_quality_reads=low_quality,
    )
