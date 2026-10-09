"""Validation of user input: accessions, sample sheets (metadata CSV) and paths.

All input is treated as untrusted. Sample names become file names, so they
are restricted to a conservative character set; accessions must be SRA *run*
accessions (SRR/ERR/DRR), the only identifiers ``prefetch`` can turn
directly into a single sequencing run.
"""

from __future__ import annotations

import csv
import io
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from .errors import InputValidationError

LAYOUTS = ("paired", "single")
_LAYOUT_ALIASES = {
    "paired": "paired",
    "pe": "paired",
    "paired-end": "paired",
    "paired_end": "paired",
    "single": "single",
    "se": "single",
    "single-end": "single",
    "single_end": "single",
}

RUN_ACCESSION_RE = re.compile(r"^(SRR|ERR|DRR)[0-9]{6,12}$")
SAMPLE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,99}$")

#: Recognised identifiers that are *not* run accessions, with guidance.
_NON_RUN_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"^GSE[0-9]+$"), "a GEO Series (a whole study)"),
    (re.compile(r"^GSM[0-9]+$"), "a GEO Sample"),
    (re.compile(r"^GDS[0-9]+$"), "a GEO DataSet"),
    (re.compile(r"^(SRP|ERP|DRP)[0-9]+$"), "an SRA Study"),
    (re.compile(r"^(SRX|ERX|DRX)[0-9]+$"), "an SRA Experiment"),
    (re.compile(r"^(SRS|ERS|DRS)[0-9]+$"), "an SRA Sample"),
    (re.compile(r"^(SRA|ERA|DRA)[0-9]+$"), "an SRA Submission"),
    (re.compile(r"^PRJ(NA|EB|DB)[0-9]+$"), "a BioProject"),
    (re.compile(r"^SAM(N|EA|D)[0-9]+$"), "a BioSample"),
]
_RUN_HINT = (
    "Only SRA run accessions (SRR..., ERR..., DRR...) are accepted. "
    "Find the runs of a study with the NCBI SRA Run Selector "
    "(https://www.ncbi.nlm.nih.gov/Traces/study/) and list them in a metadata CSV."
)

REQUIRED_COLUMNS = ("sample", "accession")
OPTIONAL_KNOWN_COLUMNS = ("layout", "condition")


@dataclass(frozen=True)
class Sample:
    """One sequencing run to process.

    ``name`` is the user's sample identifier (used for output file names);
    ``accession`` is the SRA run accession the data are downloaded from.
    """

    name: str
    accession: str
    layout: str
    metadata: dict[str, str] = field(default_factory=dict)

    @property
    def paired(self) -> bool:
        return self.layout == "paired"


def normalise_layout(value: str | None) -> str | None:
    if value is None:
        return None
    return _LAYOUT_ALIASES.get(value.strip().lower())


def check_accession(value: str) -> str | None:
    """Return a problem description for an invalid run accession, else ``None``."""
    acc = value.strip()
    if not acc:
        return "empty accession"
    if RUN_ACCESSION_RE.match(acc):
        return None
    upper = acc.upper()
    if RUN_ACCESSION_RE.match(upper):
        return f"'{acc}' must be upper case ('{upper}')"
    for pattern, what in _NON_RUN_PATTERNS:
        if pattern.match(upper):
            return f"'{acc}' is {what}, not an SRA run accession"
    if "://" in acc or acc.startswith("www."):
        return f"'{acc}' looks like a URL; URLs are not supported as input"
    return f"'{acc}' is not a recognised SRA run accession (expected e.g. SRR12345678)"


def check_sample_name(name: str) -> str | None:
    if not name:
        return "empty sample name"
    if not SAMPLE_NAME_RE.match(name):
        return (
            f"sample name '{name}' is not allowed: use 1-100 characters from letters, digits, "
            "'.', '_' and '-', starting with a letter or digit (names become file names)"
        )
    if name.endswith((".fastq", ".gz", ".fq")):
        return f"sample name '{name}' must not end in a file extension"
    return None


def _check_duplicates(samples: Sequence[Sample], problems: list[str]) -> None:
    seen_names: dict[str, int] = {}
    seen_lower: dict[str, str] = {}
    seen_acc: dict[str, str] = {}
    for s in samples:
        seen_names[s.name] = seen_names.get(s.name, 0) + 1
        low = s.name.lower()
        if low in seen_lower and seen_lower[low] != s.name:
            problems.append(
                f"sample names '{seen_lower[low]}' and '{s.name}' differ only in letter case; "
                "this is unsafe on case-insensitive file systems"
            )
        seen_lower.setdefault(low, s.name)
        if s.accession in seen_acc:
            problems.append(
                f"accession {s.accession} is listed twice (samples '{seen_acc[s.accession]}' "
                f"and '{s.name}'); each run must be processed once"
            )
        else:
            seen_acc[s.accession] = s.name
    for name, count in seen_names.items():
        if count > 1:
            problems.append(
                f"sample name '{name}' appears {count} times. Sample names must be unique; "
                "merging several runs into one sample is not supported in this version"
            )


def samples_from_accessions(accessions: Iterable[str], layout: str | None) -> list[Sample]:
    """Build samples from ``--accessions``; the accession doubles as the sample name."""
    problems: list[str] = []
    lay = normalise_layout(layout)
    if layout is None:
        problems.append("--layout (paired or single) is required with --accessions")
    elif lay is None:
        problems.append(f"unsupported layout '{layout}' (use paired or single)")
    samples: list[Sample] = []
    accs = [a.strip() for a in accessions]
    if not accs:
        problems.append("no accessions given")
    for acc in accs:
        problem = check_accession(acc)
        if problem:
            problems.append(problem)
        else:
            samples.append(Sample(name=acc, accession=acc, layout=lay or "paired"))
    _check_duplicates(samples, problems)
    if problems:
        raise InputValidationError("invalid accession input", problems=problems, hint=_RUN_HINT)
    return samples


def read_metadata_csv(path: Path, default_layout: str | None = None) -> list[Sample]:
    """Parse and validate a sample sheet.

    Required columns: ``sample``, ``accession``. ``layout`` is required unless
    ``default_layout`` is given (``--layout``). Any other columns (e.g.
    ``condition``) are preserved as metadata. Comma, tab and semicolon
    delimiters are detected automatically.
    """
    p = Path(path)
    if not p.exists():
        raise InputValidationError(f"metadata file not found: {p}")
    if not p.is_file():
        raise InputValidationError(f"metadata path is not a file: {p}")
    try:
        raw = p.read_bytes()
    except OSError as exc:
        raise InputValidationError(f"cannot read metadata file {p}: {exc.strerror}") from exc
    if len(raw) > 50 * 1024 * 1024:
        raise InputValidationError(f"metadata file {p} is unexpectedly large (>50 MB)")
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputValidationError(
            f"metadata file {p} is not UTF-8 text (byte {exc.start})",
            hint="Save the sheet as 'CSV UTF-8' from your spreadsheet program.",
        ) from exc
    if not text.strip():
        raise InputValidationError(f"metadata file {p} is empty")

    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    header_line = lines[0] if lines else ""
    delimiter = max((",", "\t", ";"), key=header_line.count)
    reader = csv.reader(io.StringIO("\n".join(lines)), delimiter=delimiter, strict=True)
    try:
        rows = list(reader)
    except csv.Error as exc:
        raise InputValidationError(f"metadata file {p} is not valid CSV: {exc}") from exc

    header = [h.strip() for h in rows[0]]
    lower = [h.lower() for h in header]
    problems: list[str] = []
    for col in REQUIRED_COLUMNS:
        if col not in lower:
            problems.append(f"required column '{col}' is missing (found: {', '.join(header)})")
    if "layout" not in lower and default_layout is None:
        problems.append("column 'layout' is missing; add it or pass --layout paired|single")
    if any(not h for h in header):
        problems.append("the header contains an empty column name")
    dups = sorted({h for h in lower if lower.count(h) > 1 and h})
    if dups:
        problems.append(f"duplicate column name(s): {', '.join(dups)}")
    if problems:
        raise InputValidationError(f"metadata file {p} has an invalid header", problems=problems)

    default_lay = normalise_layout(default_layout) if default_layout else None
    if default_layout is not None and default_lay is None:
        raise InputValidationError(
            f"unsupported --layout '{default_layout}' (use paired or single)"
        )

    samples: list[Sample] = []
    if len(rows) < 2:
        raise InputValidationError(f"metadata file {p} contains a header but no samples")
    for line_no, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            problems.append(f"row {line_no}: expected {len(header)} fields, found {len(row)}")
            continue
        rec = {lower[i]: row[i].strip() for i in range(len(header))}
        name, acc = rec.get("sample", ""), rec.get("accession", "")
        row_problems = []
        if not name:
            row_problems.append("empty 'sample'")
        elif (msg := check_sample_name(name)) is not None:
            row_problems.append(msg)
        if not acc:
            row_problems.append("empty 'accession'")
        elif (msg := check_accession(acc)) is not None:
            row_problems.append(msg)
        layout_value = rec.get("layout", "")
        if layout_value:
            lay = normalise_layout(layout_value)
            if lay is None:
                row_problems.append(f"unsupported layout '{layout_value}' (use paired or single)")
        else:
            lay = default_lay
            if lay is None:
                row_problems.append("empty 'layout' and no --layout default given")
        if default_lay and layout_value and lay and lay != default_lay:
            row_problems.append(f"layout '{layout_value}' conflicts with --layout {default_lay}")
        if row_problems:
            problems.extend(f"row {line_no}: {m}" for m in row_problems)
            continue
        metadata = {
            header[i]: row[i].strip()
            for i in range(len(header))
            if lower[i] not in ("sample", "accession", "layout")
        }
        assert lay is not None
        samples.append(Sample(name=name, accession=acc, layout=lay, metadata=metadata))
    _check_duplicates(samples, problems)
    if problems:
        raise InputValidationError(
            f"metadata file {p} contains {len(problems)} problem(s)",
            problems=problems,
            hint=_RUN_HINT if any("run accession" in m for m in problems) else None,
        )
    return samples


def check_output_path(path: Path) -> list[str]:
    """Static checks on an output directory path (no filesystem changes)."""
    problems: list[str] = []
    text = str(path)
    if not text.strip():
        problems.append("output directory path is empty")
        return problems
    if "\x00" in text or "\n" in text:
        problems.append("output directory path contains control characters")
    resolved = Path(path).expanduser().absolute()
    if resolved == Path("/") or resolved == Path.home():
        problems.append(
            f"refusing to use {resolved} as output directory; choose a project subdirectory"
        )
    if resolved.exists() and not resolved.is_dir():
        problems.append(f"output path {resolved} exists and is not a directory")
    return problems
