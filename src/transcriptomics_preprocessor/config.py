"""Configuration loading, merging and validation.

Precedence: packaged defaults < user ``--config`` YAML < command-line options.
All values are validated against :data:`SCHEMA`; unknown keys are rejected so
that a misspelt option can never be silently ignored.
"""

from __future__ import annotations

import copy
import os
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

from .errors import ConfigError
from .utils import parse_size


@dataclass(frozen=True)
class Field:
    """Validation rule for one configuration value."""

    kind: str  # int | float | bool | str | size | choice | path | dna
    minimum: float | None = None
    maximum: float | None = None
    choices: tuple[str, ...] = ()
    nullable: bool = False


_DNA = frozenset("ACGTN")

#: Validation rules. Keys are dotted configuration paths.
SCHEMA: dict[str, Field] = {
    "threads": Field("int", 1, 512),
    "jobs": Field("int", 1, 64),
    "tmpdir": Field("path", nullable=True),
    "on_error": Field("choice", choices=("continue", "stop")),
    "download.max_size": Field("size", 1),
    "download.retries": Field("int", 0, 20),
    "download.retry_wait_seconds": Field("float", 0, 3600),
    "download.max_parallel": Field("int", 1, 16),
    "download.timeout_hours": Field("float", 0.01, 720),
    "download.keep_sra": Field("bool"),
    "download.min_free_space": Field("size", 0),
    "conversion.space_factor": Field("float", 1, 100),
    "conversion.keep_unpaired": Field("bool"),
    "conversion.compression_level": Field("int", 1, 9),
    "conversion.timeout_hours": Field("float", 0.01, 720),
    "validation.mode": Field("choice", choices=("full", "quick")),
    "validation.quick_records": Field("int", 1000, 100_000_000),
    "fastqc.memory_mb": Field("int", 100, 10000),
    "fastqc.timeout_hours": Field("float", 0.01, 720),
    "fastp.adapter_trimming": Field("bool"),
    "fastp.adapter_sequence": Field("dna", nullable=True),
    "fastp.adapter_sequence_r2": Field("dna", nullable=True),
    "fastp.detect_adapter_for_pe": Field("bool"),
    "fastp.quality_filtering": Field("bool"),
    "fastp.qualified_quality_phred": Field("int", 0, 60),
    "fastp.unqualified_percent_limit": Field("int", 0, 100),
    "fastp.n_base_limit": Field("int", 0, 1000),
    "fastp.average_qual": Field("int", 0, 60),
    "fastp.length_filtering": Field("bool"),
    "fastp.length_required": Field("int", 1, 100000),
    "fastp.length_limit": Field("int", 0, 1000000),
    "fastp.cut_front": Field("bool"),
    "fastp.cut_tail": Field("bool"),
    "fastp.cut_window_size": Field("int", 1, 1000),
    "fastp.cut_mean_quality": Field("int", 1, 36),
    "fastp.trim_front1": Field("int", 0, 10000),
    "fastp.trim_tail1": Field("int", 0, 10000),
    "fastp.trim_front2": Field("int", 0, 10000),
    "fastp.trim_tail2": Field("int", 0, 10000),
    "fastp.trim_poly_g": Field("choice", choices=("auto", "on", "off")),
    "fastp.trim_poly_x": Field("bool"),
    "fastp.low_complexity_filter": Field("bool"),
    "fastp.complexity_threshold": Field("int", 0, 100),
    "fastp.overrepresentation_analysis": Field("bool"),
    "fastp.compression_level": Field("int", 1, 9),
    "fastp.timeout_hours": Field("float", 0.01, 720),
    "multiqc.timeout_hours": Field("float", 0.01, 720),
    "screening.min_retained_percent": Field("float", 0, 100),
    "screening.min_clean_reads": Field("int", 0, 10**12),
    "screening.max_clean_adapter_percent": Field("float", 0, 100),
    "screening.min_clean_q30_rate": Field("float", 0, 1),
    "output.checksums": Field("bool"),
}


def load_defaults() -> dict[str, Any]:
    """Return the packaged default configuration as a nested dictionary."""
    text = (
        resources.files("transcriptomics_preprocessor")
        .joinpath("data/default.yaml")
        .read_text(encoding="utf-8")
    )
    data = yaml.safe_load(text)
    if not isinstance(data, dict):  # pragma: no cover - packaging error
        raise ConfigError("packaged default configuration is malformed")
    return data


def load_yaml_file(path: Path) -> dict[str, Any]:
    """Load a user configuration file (YAML mapping)."""
    p = Path(path)
    if not p.is_file():
        raise ConfigError(f"configuration file not found: {p}")
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"configuration file {p} is not valid YAML: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise ConfigError(f"configuration file {p} is not UTF-8 text") from exc
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"configuration file {p} must contain a YAML mapping (key: value)")
    return data


def flatten(data: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    """``{'a': {'b': 1}} -> {'a.b': 1}``."""
    out: dict[str, Any] = {}
    for key, value in data.items():
        dotted = f"{prefix}{key}"
        if isinstance(value, dict):
            out.update(flatten(value, dotted + "."))
        else:
            out[dotted] = value
    return out


def unflatten(flat: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for dotted, value in flat.items():
        node = out
        *parents, leaf = dotted.split(".")
        for part in parents:
            node = node.setdefault(part, {})
        node[leaf] = value
    return out


def _coerce(key: str, value: Any, rule: Field) -> Any:
    """Validate (and minimally coerce) one value; raise ``ValueError`` with a reason."""
    if value is None:
        if rule.nullable:
            return None
        raise ValueError("a value is required")
    kind = rule.kind
    if kind == "bool":
        if isinstance(value, bool):
            return value
        raise ValueError(f"expected true/false, got {value!r}")
    if kind in ("int", "float"):
        if isinstance(value, bool):
            raise ValueError(f"expected a number, got {value!r}")
        if kind == "int":
            if isinstance(value, float) and value.is_integer():
                value = int(value)
            if not isinstance(value, int):
                raise ValueError(f"expected an integer, got {value!r}")
        else:
            if not isinstance(value, (int, float)):
                raise ValueError(f"expected a number, got {value!r}")
            value = float(value)
            if value != value:  # NaN
                raise ValueError("NaN is not allowed")
        if rule.minimum is not None and value < rule.minimum:
            raise ValueError(f"must be >= {rule.minimum:g}, got {value}")
        if rule.maximum is not None and value > rule.maximum:
            raise ValueError(f"must be <= {rule.maximum:g}, got {value}")
        return value
    if kind == "choice":
        if not isinstance(value, str) or value not in rule.choices:
            raise ValueError(f"must be one of {', '.join(rule.choices)}; got {value!r}")
        return value
    if kind == "size":
        if isinstance(value, bool):
            raise ValueError(f"expected a size such as 50G, got {value!r}")
        try:
            nbytes = parse_size(str(value))
        except ValueError:
            raise ValueError(f"expected a size such as 500M or 50G, got {value!r}") from None
        if rule.minimum is not None and nbytes < rule.minimum:
            raise ValueError(f"must be at least {rule.minimum:g} bytes")
        return str(value).strip()
    if kind == "path":
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"expected a path, got {value!r}")
        if "\x00" in value:
            raise ValueError("path contains a NUL character")
        return value
    if kind == "dna":
        if not isinstance(value, str):
            raise ValueError(f"expected a DNA sequence, got {value!r}")
        seq = value.strip().upper()
        if len(seq) < 6 or not set(seq) <= _DNA:
            raise ValueError("adapter must be at least 6 bases using only A, C, G, T, N")
        return seq
    if kind == "str":
        if not isinstance(value, str):
            raise ValueError(f"expected text, got {value!r}")
        return value
    raise AssertionError(f"unknown schema kind {kind} for {key}")  # pragma: no cover


def validate(flat: dict[str, Any]) -> dict[str, Any]:
    """Validate a flattened configuration; collect *all* problems before raising."""
    problems: list[str] = []
    clean: dict[str, Any] = {}
    for key, value in flat.items():
        rule = SCHEMA.get(key)
        if rule is None:
            problems.append(f"unknown configuration key '{key}'")
            continue
        try:
            clean[key] = _coerce(key, value, rule)
        except ValueError as exc:
            problems.append(f"{key}: {exc}")
    missing = sorted(set(SCHEMA) - set(flat))
    problems.extend(f"{key}: missing" for key in missing)

    if not problems:
        if clean["jobs"] > clean["threads"]:
            problems.append(
                f"jobs ({clean['jobs']}) cannot exceed threads ({clean['threads']}): "
                "every concurrent sample needs at least one thread"
            )
        if (
            clean["fastp.length_limit"]
            and clean["fastp.length_limit"] < clean["fastp.length_required"]
        ):
            problems.append("fastp.length_limit must be 0 (no limit) or >= fastp.length_required")
        if clean["fastp.adapter_sequence_r2"] and not clean["fastp.adapter_sequence"]:
            problems.append(
                "fastp.adapter_sequence_r2 requires fastp.adapter_sequence to be set too"
            )
        if not clean["fastp.adapter_trimming"] and (
            clean["fastp.adapter_sequence"] or clean["fastp.adapter_sequence_r2"]
        ):
            problems.append("adapter sequences were given but fastp.adapter_trimming is false")
    if problems:
        raise ConfigError("invalid configuration", problems=problems)
    return clean


@dataclass(frozen=True)
class PipelineConfig:
    """Validated configuration. ``values`` maps dotted keys to clean values."""

    values: dict[str, Any]

    def __getitem__(self, key: str) -> Any:
        return self.values[key]

    @property
    def threads(self) -> int:
        return int(self.values["threads"])

    @property
    def jobs(self) -> int:
        return int(self.values["jobs"])

    @property
    def threads_per_job(self) -> int:
        return max(1, self.threads // self.jobs)

    def section(self, name: str) -> dict[str, Any]:
        prefix = name + "."
        return {k[len(prefix) :]: v for k, v in self.values.items() if k.startswith(prefix)}

    def as_nested(self) -> dict[str, Any]:
        return unflatten(dict(self.values))

    def to_yaml(self) -> str:
        return str(yaml.safe_dump(self.as_nested(), sort_keys=False, default_flow_style=False))


def build_config(
    user_file: Path | None = None, overrides: dict[str, Any] | None = None
) -> PipelineConfig:
    """Merge defaults, an optional YAML file and CLI overrides, then validate."""
    flat = flatten(load_defaults())
    if user_file is not None:
        user_flat = flatten(load_yaml_file(user_file))
        unknown = sorted(k for k in user_flat if k not in SCHEMA)
        if unknown:
            raise ConfigError(
                f"unknown key(s) in {user_file}",
                problems=[f"unknown configuration key '{k}'" for k in unknown],
                hint="Compare with config/default.yaml in the repository for valid keys.",
            )
        flat.update(user_flat)
    for key, value in (overrides or {}).items():
        if value is not None:
            flat[key] = value
    return PipelineConfig(copy.deepcopy(validate(flat)))


def cpu_count() -> int:
    """CPUs available to this process (respects affinity/cgroup restrictions where exposed)."""
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):  # pragma: no cover - non-Linux
        return os.cpu_count() or 1
