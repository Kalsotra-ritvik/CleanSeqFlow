from __future__ import annotations

from pathlib import Path

import pytest

from transcriptomics_preprocessor.config import SCHEMA, build_config, flatten, load_defaults
from transcriptomics_preprocessor.errors import ConfigError


def test_defaults_cover_schema_exactly() -> None:
    assert set(flatten(load_defaults())) == set(SCHEMA)


def test_defaults_are_valid_and_match_fastp_defaults() -> None:
    cfg = build_config()
    # fastp 1.3.6 documented defaults (fastp --help): -q 15 -u 40 -n 5 -l 15, cut window 4, mean quality 20
    assert cfg["fastp.qualified_quality_phred"] == 15
    assert cfg["fastp.unqualified_percent_limit"] == 40
    assert cfg["fastp.n_base_limit"] == 5
    assert cfg["fastp.length_required"] == 15
    assert cfg["fastp.cut_window_size"] == 4 and cfg["fastp.cut_mean_quality"] == 20
    assert cfg["fastp.cut_front"] is False and cfg["fastp.cut_tail"] is False
    assert cfg.threads == 4 and cfg.jobs == 1 and cfg.threads_per_job == 4


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("threads", 0),
        ("threads", 513),
        ("threads", "four"),
        ("threads", 2.5),
        ("threads", True),
        ("jobs", 0),
        ("fastp.qualified_quality_phred", -1),
        ("fastp.qualified_quality_phred", 61),
        ("fastp.unqualified_percent_limit", 101),
        ("fastp.length_required", 0),
        ("fastp.cut_mean_quality", 0),
        ("fastp.cut_mean_quality", 37),
        ("fastp.cut_window_size", 1001),
        ("fastp.trim_poly_g", "maybe"),
        ("on_error", "ignore"),
        ("validation.mode", "deep"),
        ("download.max_size", "lots"),
        ("download.max_size", "-5G"),
        ("download.retries", 21),
        ("fastp.adapter_sequence", "AGATXX"),
        ("fastp.adapter_sequence", "ACG"),
        ("fastp.adapter_trimming", "yes"),
        ("screening.min_clean_q30_rate", 1.5),
        ("conversion.compression_level", 10),
        ("threads", None),
    ],
)
def test_out_of_range_and_wrong_types_rejected(key: str, value: object) -> None:
    with pytest.raises(ConfigError) as e:
        build_config(overrides={key: value}) if value is not None else _set_none(key)
    assert any(key in p for p in e.value.problems)


def _set_none(key: str) -> None:
    from transcriptomics_preprocessor.config import validate

    flat = flatten(load_defaults())
    flat[key] = None
    validate(flat)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("threads", 1),
        ("threads", 512),
        ("fastp.qualified_quality_phred", 0),
        ("fastp.qualified_quality_phred", 60),
        ("fastp.unqualified_percent_limit", 0),
        ("fastp.unqualified_percent_limit", 100),
        ("fastp.length_required", 1),
        ("fastp.cut_mean_quality", 1),
        ("fastp.cut_mean_quality", 36),
        ("download.retries", 0),
        ("download.max_size", "500M"),
        ("download.max_size", "1T"),
        ("fastp.adapter_sequence", "agatcggaagagc"),
        ("conversion.compression_level", 1),
        ("conversion.compression_level", 9),
    ],
)
def test_boundary_values_accepted(key: str, value: object) -> None:
    cfg = build_config(overrides={key: value})
    if key == "fastp.adapter_sequence":
        assert cfg[key] == "AGATCGGAAGAGC"
    else:
        assert cfg[key] == value


def test_incompatible_combinations() -> None:
    with pytest.raises(ConfigError, match="invalid"):
        build_config(overrides={"threads": 2, "jobs": 3})
    with pytest.raises(ConfigError) as e:
        build_config(overrides={"fastp.length_required": 50, "fastp.length_limit": 30})
    assert any("length_limit" in p for p in e.value.problems)
    with pytest.raises(ConfigError) as e:
        build_config(overrides={"fastp.adapter_sequence_r2": "AGATCGGAAGAGCGTCGTGT"})
    assert any("requires" in p for p in e.value.problems)
    with pytest.raises(ConfigError) as e:
        build_config(
            overrides={"fastp.adapter_trimming": False, "fastp.adapter_sequence": "AGATCGGAAGAGC"}
        )
    assert any("adapter_trimming is false" in p for p in e.value.problems)


def test_threads_per_job_division() -> None:
    assert build_config(overrides={"threads": 8, "jobs": 3}).threads_per_job == 2
    assert build_config(overrides={"threads": 3, "jobs": 3}).threads_per_job == 1


def test_user_file_precedence_and_unknown_keys(tmp_path: Path) -> None:
    f = tmp_path / "c.yaml"
    f.write_text("threads: 6\nfastp:\n  length_required: 25\n")
    cfg = build_config(f, {"threads": 2})
    assert cfg.threads == 2 and cfg["fastp.length_required"] == 25
    f.write_text("fastp:\n  lenght_required: 25\n")
    with pytest.raises(ConfigError) as e:
        build_config(f)
    assert any("lenght_required" in p for p in e.value.problems)


@pytest.mark.parametrize("content", ["threads: [1, 2", "- a\n- b\n"])
def test_bad_yaml(tmp_path: Path, content: str) -> None:
    f = tmp_path / "c.yaml"
    f.write_text(content)
    with pytest.raises(ConfigError):
        build_config(f)


def test_empty_yaml_is_defaults(tmp_path: Path) -> None:
    f = tmp_path / "c.yaml"
    f.write_text("")
    assert build_config(f).threads == 4
    with pytest.raises(ConfigError, match="not found"):
        build_config(tmp_path / "missing.yaml")


def test_yaml_roundtrip() -> None:
    import yaml

    cfg = build_config(overrides={"threads": 8})
    again = yaml.safe_load(cfg.to_yaml())
    assert again["threads"] == 8 and again["fastp"]["length_required"] == 15
