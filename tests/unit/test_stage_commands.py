"""Command construction for external tools, including systematic option combinations."""

from __future__ import annotations

import itertools
from pathlib import Path

import pytest

from transcriptomics_preprocessor import conversion, download, fastp, fastqc, multiqc
from transcriptomics_preprocessor.config import build_config
from transcriptomics_preprocessor.console import Console
from transcriptomics_preprocessor.context import StageContext
from transcriptomics_preprocessor.errors import IntegrityError
from transcriptomics_preprocessor.layout import OutputLayout
from transcriptomics_preprocessor.validation import Sample

PE = Sample("s1", "SRR1234567", "paired")
SE = Sample("s2", "SRR7654321", "single")


def _fastp_cmd(sample: Sample, **overrides: object) -> list[str]:
    cfg = build_config(overrides=overrides)
    ins = [Path("in_R1.fq.gz"), Path("in_R2.fq.gz")][: 2 if sample.paired else 1]
    outs = [Path("o_R1.fq.gz"), Path("o_R2.fq.gz")][: 2 if sample.paired else 1]
    return fastp.build_command(sample, ins, outs, Path("j.json"), Path("h.html"), cfg, 8)


def _opt(cmd: list[str], flag: str) -> str | None:
    return cmd[cmd.index(flag) + 1] if flag in cmd else None


def test_fastp_defaults_paired() -> None:
    cmd = _fastp_cmd(PE)
    assert _opt(cmd, "--in2") == "in_R2.fq.gz" and _opt(cmd, "--out2") == "o_R2.fq.gz"
    assert "--detect_adapter_for_pe" in cmd
    assert _opt(cmd, "--qualified_quality_phred") == "15"
    assert _opt(cmd, "--unqualified_percent_limit") == "40"
    assert _opt(cmd, "--n_base_limit") == "5"
    assert _opt(cmd, "--length_required") == "15"
    assert _opt(cmd, "--thread") == "8"
    for absent in (
        "--dedup",
        "-D",
        "--cut_front",
        "--cut_tail",
        "--trim_poly_g",
        "--disable_trim_poly_g",
        "--length_limit",
        "--average_qual",
        "--trim_front1",
        "--trim_front2",
    ):
        assert absent not in cmd  # duplicates are never removed; fastp defaults otherwise


def test_fastp_single_end_has_no_r2_options() -> None:
    cmd = _fastp_cmd(SE, **{"fastp.adapter_sequence": "AGATCGGAAGAGCACACGTCT"})
    assert "--in2" not in cmd and "--out2" not in cmd and "--detect_adapter_for_pe" not in cmd
    assert _opt(cmd, "--adapter_sequence") == "AGATCGGAAGAGCACACGTCT"


def test_fastp_thread_cap() -> None:
    cfg = build_config(overrides={"threads": 64})
    cmd = fastp.build_command(
        PE, [Path("a"), Path("b")], [Path("c"), Path("d")], Path("j"), Path("h"), cfg, 64
    )
    assert _opt(cmd, "--thread") == "16"


def test_fastp_r2_trimming_is_explicit() -> None:
    cmd = _fastp_cmd(PE, **{"fastp.trim_front1": 3})
    assert _opt(cmd, "--trim_front1") == "3" and _opt(cmd, "--trim_front2") == "0"


BOOL_SWITCHES = {
    "fastp.adapter_trimming": "--disable_adapter_trimming",
    "fastp.quality_filtering": "--disable_quality_filtering",
    "fastp.length_filtering": "--disable_length_filtering",
}


@pytest.mark.parametrize(
    "combo",
    list(
        itertools.product(
            [True, False],
            [True, False],
            [True, False],
            ["auto", "on", "off"],
            [False, True],
            [False, True],
        )
    ),
)
def test_fastp_option_combinations(combo: tuple[bool, bool, bool, str, bool, bool]) -> None:
    """Exhaustive over 3 enable switches x polyG mode x cut_front x cut_tail (96 combinations)."""
    adapt, qual, length, polyg, cut_front, cut_tail = combo
    for sample in (PE, SE):
        cmd = _fastp_cmd(
            sample,
            **{
                "fastp.adapter_trimming": adapt,
                "fastp.quality_filtering": qual,
                "fastp.length_filtering": length,
                "fastp.trim_poly_g": polyg,
                "fastp.cut_front": cut_front,
                "fastp.cut_tail": cut_tail,
            },
        )
        assert ("--disable_adapter_trimming" in cmd) == (not adapt)
        assert ("--disable_quality_filtering" in cmd) == (not qual)
        assert ("--qualified_quality_phred" in cmd) == qual
        assert ("--disable_length_filtering" in cmd) == (not length)
        assert ("--length_required" in cmd) == length
        assert ("--trim_poly_g" in cmd) == (polyg == "on")
        assert ("--disable_trim_poly_g" in cmd) == (polyg == "off")
        assert ("--cut_front" in cmd) == cut_front and ("--cut_tail" in cmd) == cut_tail
        assert ("--cut_mean_quality" in cmd) == (cut_front or cut_tail)
        assert ("--detect_adapter_for_pe" in cmd) == (adapt and sample.paired)
        # descriptions mention every enabled/disabled operation
        text = " ".join(
            fastp.describe_parameters(
                build_config(
                    overrides={
                        "fastp.adapter_trimming": adapt,
                        "fastp.quality_filtering": qual,
                        "fastp.length_filtering": length,
                    }
                ),
                sample.paired,
            )
        )
        assert ("DISABLED" in text) == (not (adapt and qual and length))
        assert "Deduplication: OFF" in text


def test_fastp_stage_parameters_ignore_irrelevant_r2_for_single_end() -> None:
    a = fastp.stage_parameters(build_config(overrides={"fastp.trim_front2": 0}), paired=False)
    b = fastp.stage_parameters(build_config(overrides={"fastp.trim_front2": 7}), paired=False)
    assert a == b
    assert "timeout_hours" not in a


def test_fastqc_command_and_stem() -> None:
    cmd = fastqc.build_command(
        [Path("x/a_R1.fastq.gz"), Path("x/a_R2.fastq.gz")], Path("o"), Path("t"), 8, 512
    )
    assert _opt(cmd, "--threads") == "2" and _opt(cmd, "--dir") == "t" and "--noextract" in cmd
    assert OutputLayout.fastqc_stem(Path("a_R1.fastq.gz")) == "a_R1"
    assert OutputLayout.fastqc_stem(Path("se.fq")) == "se"


def test_multiqc_command_has_privacy_flags() -> None:
    cmd = multiqc.build_command(Path("l.txt"), Path("o"), "r", "t")
    for flag in ("--no-version-check", "--no-ai", "--no-megaqc-upload", "--file-list", "--force"):
        assert flag in cmd


def _ctx(tmp_path: Path, **ov: object) -> StageContext:
    cfg = build_config(overrides=ov)
    return StageContext(OutputLayout(tmp_path), cfg, tmp_path / "tmp", Console(quiet=True), {}, 2)


def test_prefetch_command(tmp_path: Path) -> None:
    cmd = download.build_command(PE, _ctx(tmp_path, **{"download.max_size": "50G"}))
    assert cmd[:2] == ["prefetch", "SRR1234567"]
    assert _opt(cmd, "--max-size") == str(50 * 1024 * 1024)  # KB
    assert _opt(cmd, "--verify") == "yes" and _opt(cmd, "--resume") == "yes"


@pytest.mark.parametrize(
    ("text", "kind"),
    [
        # Real sra-tools 3.4.1 output when NCBI is unreachable (observed 2026-10-09):
        (
            "prefetch.3.4.1 int: connection not found while validating within network system module - "
            "cannot resolve remote location of 'SRR26053372'",
            "network",
        ),
        ("curl: (28) Operation timed out", "network"),
        ("err: item not found while resolving 'SRR0'", "permanent"),
        ("requires authorization (dbGaP)", "permanent"),
        ("something odd happened", "unknown"),
    ],
)
def test_prefetch_failure_classification(text: str, kind: str) -> None:
    assert download.classify_failure(text) == kind


def test_fasterq_dump_command() -> None:
    cmd = conversion.build_command(Path("a.sra"), Path("w"), Path("t"), 0)
    assert "--split-3" in cmd and "--skip-technical" in cmd and _opt(cmd, "--threads") == "1"


def test_classify_fasterq_outputs(tmp_path: Path) -> None:
    for name in ("SRR1_1.fastq", "SRR1_2.fastq", "SRR1.fastq"):
        (tmp_path / name).write_text("x")
    r1, r2, single = conversion.classify_outputs(tmp_path)
    assert r1 and r1.name == "SRR1_1.fastq" and r2 and single and single.name == "SRR1.fastq"
    (tmp_path / "SRR1_3.fastq").write_text("x")
    with pytest.raises(IntegrityError, match="additional read files"):
        conversion.classify_outputs(tmp_path)
