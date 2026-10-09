from __future__ import annotations

from pathlib import Path

import pytest

from transcriptomics_preprocessor.errors import InputValidationError
from transcriptomics_preprocessor.validation import (
    check_accession,
    check_output_path,
    check_sample_name,
    read_metadata_csv,
    samples_from_accessions,
)


@pytest.mark.parametrize("acc", ["SRR26053372", "ERR123456", "DRR000001", "SRR123456789012"])
def test_valid_run_accessions(acc: str) -> None:
    assert check_accession(acc) is None


@pytest.mark.parametrize(
    ("acc", "fragment"),
    [
        ("GSE243105", "GEO Series"),
        ("GSM123", "GEO Sample"),
        ("SRP123456", "SRA Study"),
        ("SRX123456", "SRA Experiment"),
        ("SRS123456", "SRA Sample"),
        ("PRJNA123456", "BioProject"),
        ("SAMN123456", "BioSample"),
        ("srr26053372", "upper case"),
        ("SRR123", "not a recognised"),  # too few digits
        ("https://example.org/x.fastq", "URL"),
        ("SRR26053372;rm -rf /", "not a recognised"),
        ("", "empty"),
    ],
)
def test_invalid_accessions_are_explained(acc: str, fragment: str) -> None:
    msg = check_accession(acc)
    assert msg is not None and fragment in msg


@pytest.mark.parametrize("name", ["ctrl_1", "S1", "a.b-c_d", "1abc", "x" * 100])
def test_valid_sample_names(name: str) -> None:
    assert check_sample_name(name) is None


@pytest.mark.parametrize(
    "name", ["", "_x", ".hidden", "a b", "a/b", "../x", "x" * 101, "s1.fastq", "a;b", "ü"]
)
def test_invalid_sample_names(name: str) -> None:
    assert check_sample_name(name) is not None


def test_accessions_require_layout() -> None:
    with pytest.raises(InputValidationError) as e:
        samples_from_accessions(["SRR26053372"], None)
    assert any("--layout" in p for p in e.value.problems)


def test_accessions_duplicates_rejected() -> None:
    with pytest.raises(InputValidationError) as e:
        samples_from_accessions(["SRR26053372", "SRR26053372"], "paired")
    assert any("listed twice" in p for p in e.value.problems)


def test_accessions_ok() -> None:
    s = samples_from_accessions(["SRR26053372", "ERR000111"], "single")
    assert [x.name for x in s] == ["SRR26053372", "ERR000111"]
    assert all(x.layout == "single" and not x.paired for x in s)


def _csv(tmp_path: Path, text: str, name: str = "m.csv") -> Path:
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def test_metadata_example_from_specification(tmp_path: Path) -> None:
    p = _csv(
        tmp_path,
        "sample,accession,layout,condition\n"
        "SRR26053372,SRR26053372,paired,siOGT\n"
        "SRR26053373,SRR26053373,paired,siOGT\n"
        "SRR26053375,SRR26053375,paired,siControl\n",
    )
    samples = read_metadata_csv(p)
    assert [s.accession for s in samples] == ["SRR26053372", "SRR26053373", "SRR26053375"]
    assert samples[2].metadata == {"condition": "siControl"}


def test_condition_optional_and_extra_columns_preserved(tmp_path: Path) -> None:
    p = _csv(tmp_path, "Sample,Accession,Layout,batch,donor\nA,SRR1234567,PE,1,d1\n")
    (s,) = read_metadata_csv(p)
    assert s.layout == "paired" and s.metadata == {"batch": "1", "donor": "d1"}


def test_tab_and_semicolon_delimiters_and_bom(tmp_path: Path) -> None:
    t = _csv(tmp_path, "sample\taccession\tlayout\nA\tSRR1234567\tsingle\n", "t.tsv")
    assert read_metadata_csv(t)[0].layout == "single"
    s = tmp_path / "s.csv"
    s.write_bytes("﻿sample;accession;layout\nA;SRR1234567;paired\n".encode())
    assert read_metadata_csv(s)[0].name == "A"


def test_layout_default_from_cli(tmp_path: Path) -> None:
    p = _csv(tmp_path, "sample,accession\nA,SRR1234567\n")
    assert read_metadata_csv(p, default_layout="paired")[0].paired
    with pytest.raises(InputValidationError):
        read_metadata_csv(p)  # no layout column and no default


def test_layout_conflict_with_cli(tmp_path: Path) -> None:
    p = _csv(tmp_path, "sample,accession,layout\nA,SRR1234567,single\n")
    with pytest.raises(InputValidationError) as e:
        read_metadata_csv(p, default_layout="paired")
    assert any("conflicts" in m for m in e.value.problems)


def test_all_problems_reported_together(tmp_path: Path) -> None:
    p = _csv(
        tmp_path,
        "sample,accession,layout\n"
        "A,SRR1234567,paired\n"
        "A,SRR7654321,paired\n"  # duplicate sample
        "B,GSE243105,paired\n"  # GEO id
        "C,SRR1234567,paired\n"  # duplicate accession
        "D,SRR1111111,triple\n"  # bad layout
        ",SRR2222222,paired\n"  # empty sample
        "E,,paired\n"  # empty accession
        "F,SRR3333333\n",  # wrong field count
    )
    with pytest.raises(InputValidationError) as e:
        read_metadata_csv(p)
    text = "\n".join(e.value.problems)
    for fragment in (
        "appears 2 times",
        "GEO Series",
        "listed twice",
        "unsupported layout",
        "empty 'sample'",
        "empty 'accession'",
        "expected 3 fields",
    ):
        assert fragment in text, fragment


@pytest.mark.parametrize(
    ("content", "fragment"),
    [
        ("", "empty"),
        ("sample,accession,layout\n", "no samples"),
        ("name,run\nA,SRR1234567\n", "required column 'sample'"),
        ("sample,accession,layout,layout\nA,SRR1234567,paired,paired\n", "duplicate column"),
        ('sample,accession,layout\n"A,SRR1234567,paired\n', "valid CSV"),
    ],
)
def test_malformed_csv(tmp_path: Path, content: str, fragment: str) -> None:
    p = _csv(tmp_path, content)
    with pytest.raises(InputValidationError) as e:
        read_metadata_csv(p)
    assert fragment in (e.value.message + " ".join(e.value.problems))


def test_non_utf8_and_missing_file(tmp_path: Path) -> None:
    p = tmp_path / "latin1.csv"
    p.write_bytes("sample,accession,layout\nGr\xfcn,SRR1234567,paired\n".encode("latin-1"))
    with pytest.raises(InputValidationError, match="UTF-8"):
        read_metadata_csv(p)
    with pytest.raises(InputValidationError, match="not found"):
        read_metadata_csv(tmp_path / "nope.csv")
    with pytest.raises(InputValidationError, match="not a file"):
        read_metadata_csv(tmp_path)


def test_case_insensitive_name_clash(tmp_path: Path) -> None:
    p = _csv(tmp_path, "sample,accession,layout\nCtrl,SRR1234567,paired\nctrl,SRR7654321,paired\n")
    with pytest.raises(InputValidationError) as e:
        read_metadata_csv(p)
    assert any("letter case" in m for m in e.value.problems)


def test_comment_and_blank_lines_ignored(tmp_path: Path) -> None:
    p = _csv(tmp_path, "# my study\nsample,accession,layout\n\nA,SRR1234567,paired\n\n")
    assert len(read_metadata_csv(p)) == 1


def test_output_path_checks(tmp_path: Path) -> None:
    assert check_output_path(tmp_path / "new dir with spaces") == []
    f = tmp_path / "file"
    f.write_text("x")
    assert any("not a directory" in m for m in check_output_path(f))
    assert any("refusing" in m for m in check_output_path(Path("/")))
    assert any("refusing" in m for m in check_output_path(Path.home()))
