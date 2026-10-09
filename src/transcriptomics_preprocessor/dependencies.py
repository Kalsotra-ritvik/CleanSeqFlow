"""Discovery and version detection of the external programs the pipeline needs."""

from __future__ import annotations

import re
import shutil
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .errors import DependencyError
from .runner import probe_output


@dataclass(frozen=True)
class ToolSpec:
    """How to find and identify one external program."""

    key: str
    executable: str
    version_args: tuple[str, ...]
    version_regex: str
    purpose: str
    #: Version the conda environment pins; a different version is a warning.
    expected: str
    #: Pipeline stages that need this tool.
    stages: tuple[str, ...]


TOOLS: dict[str, ToolSpec] = {
    spec.key: spec
    for spec in (
        ToolSpec(
            "prefetch",
            "prefetch",
            ("--version",),
            r"prefetch\s*:\s*([0-9][0-9.]*)",
            "download SRA runs (SRA Toolkit)",
            "3.4.1",
            ("download",),
        ),
        ToolSpec(
            "fasterq-dump",
            "fasterq-dump",
            ("--version",),
            r"fasterq-dump\s*:\s*([0-9][0-9.]*)",
            "convert SRA runs to FASTQ (SRA Toolkit)",
            "3.4.1",
            ("convert",),
        ),
        ToolSpec(
            "pigz",
            "pigz",
            ("--version",),
            r"pigz\s+([0-9][0-9.]*)",
            "parallel gzip compression of FASTQ",
            "2.8",
            ("convert",),
        ),
        ToolSpec(
            "fastqc",
            "fastqc",
            ("--version",),
            r"FastQC\s+v?([0-9][0-9.]*)",
            "per-file read quality reports",
            "0.12.1",
            ("fastqc_raw", "fastqc_clean"),
        ),
        ToolSpec(
            "fastp",
            "fastp",
            ("--version",),
            r"fastp\s+([0-9][0-9.]*)",
            "adapter trimming and quality filtering",
            "1.3.6",
            ("fastp",),
        ),
        ToolSpec(
            "multiqc",
            "multiqc",
            ("--version",),
            r"version\s+([0-9][0-9.a-z]*)",
            "aggregate QC reports across samples",
            "1.35",
            ("multiqc_raw", "multiqc_clean"),
        ),
    )
}


@dataclass(frozen=True)
class ToolStatus:
    key: str
    found: bool
    path: str | None
    version: str | None
    expected: str
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.found and self.version is not None and self.error is None

    @property
    def version_matches(self) -> bool:
        return self.version == self.expected

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["ok"] = self.ok
        d["version_matches_pin"] = self.version_matches
        return d


def detect_tool(spec: ToolSpec) -> ToolStatus:
    """Locate *spec* on ``PATH`` and determine its version by running it."""
    path = shutil.which(spec.executable)
    if path is None:
        return ToolStatus(spec.key, False, None, None, spec.expected, "not found on PATH")
    code, out = probe_output([path, *spec.version_args], timeout=120)
    match = re.search(spec.version_regex, out)
    if match is None:
        snippet = " ".join(out.split())[:200]
        return ToolStatus(
            spec.key,
            True,
            path,
            None,
            spec.expected,
            f"could not determine version (exit code {code}): {snippet}",
        )
    return ToolStatus(spec.key, True, path, match.group(1).rstrip("."), spec.expected)


def detect_all(keys: list[str] | None = None) -> dict[str, ToolStatus]:
    return {k: detect_tool(TOOLS[k]) for k in (keys or list(TOOLS))}


def require(statuses: dict[str, ToolStatus], keys: list[str]) -> None:
    """Raise :class:`DependencyError` listing every unusable tool among *keys*."""
    problems = []
    for key in keys:
        st = statuses.get(key) or detect_tool(TOOLS[key])
        if not st.ok:
            problems.append(f"{key}: {st.error or 'unusable'}")
    if problems:
        raise DependencyError(
            "required external program(s) are missing or broken: " + "; ".join(problems),
            hint="Run ./install.sh (or ./install.sh --repair) and then "
            "'transcriptomics-preprocess doctor'. If you run the Python package directly, "
            "activate the environment created by install.sh first.",
        )


def missing_message(statuses: dict[str, ToolStatus], keys: list[str]) -> list[str]:
    return [f"{k}: {statuses[k].error}" for k in keys if k in statuses and not statuses[k].ok]


def environment_info() -> dict[str, Any]:
    """Interpreter and environment details recorded for reproducibility."""
    import platform

    from . import __version__

    prefix = Path(sys.prefix)
    return {
        "transcriptomics_preprocessor": __version__,
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "environment_prefix": str(prefix),
        "conda_environment": (prefix / "conda-meta").is_dir(),
        "platform": platform.platform(),
        "machine": platform.machine(),
    }
