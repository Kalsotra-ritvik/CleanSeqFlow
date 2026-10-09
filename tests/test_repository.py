"""Consistency checks between code, documentation and repository files."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

from transcriptomics_preprocessor.dependencies import TOOLS
from transcriptomics_preprocessor.explain import render_explanation

ROOT = Path(__file__).resolve().parent.parent


def test_config_default_matches_packaged_default() -> None:
    packaged = ROOT / "src" / "transcriptomics_preprocessor" / "data" / "default.yaml"
    assert (ROOT / "config" / "default.yaml").read_text() == packaged.read_text()


def test_pipeline_explained_doc_matches_cli_text() -> None:
    doc = (ROOT / "docs" / "pipeline_explained.md").read_text()
    assert render_explanation("all").strip() in doc


def test_example_files_are_valid() -> None:
    from transcriptomics_preprocessor.config import build_config
    from transcriptomics_preprocessor.validation import read_metadata_csv

    assert len(read_metadata_csv(ROOT / "examples" / "metadata.csv")) == 3
    cfg = build_config(ROOT / "examples" / "config.yaml")
    assert cfg.jobs == 2 and cfg.threads_per_job == 4


def test_pinned_versions_agree_everywhere() -> None:
    env_yml = (ROOT / "environment.yml").read_text()
    lock = (ROOT / "environment-lock-linux-64.txt").read_text()
    readme = (ROOT / "README.md").read_text()
    names = {
        "prefetch": "sra-tools",
        "fasterq-dump": "sra-tools",
        "pigz": "pigz",
        "fastqc": "fastqc",
        "fastp": "fastp",
        "multiqc": "multiqc",
    }
    for key, spec in TOOLS.items():
        pkg = names[key]
        assert re.search(
            rf"^\s*- {re.escape(pkg)}={re.escape(spec.expected)}\s*$", env_yml, re.M
        ), key
        assert re.search(rf"/{re.escape(pkg)}-{re.escape(spec.expected)}-", lock), key
        assert spec.expected in readme, key


def test_lock_file_entries_have_checksums_and_no_defaults_channel() -> None:
    urls = [
        ln
        for ln in (ROOT / "environment-lock-linux-64.txt").read_text().splitlines()
        if ln.startswith("http")
    ]
    assert len(urls) > 50
    assert all(re.search(r"#[0-9a-f]{32}$", u) for u in urls)
    assert all(
        u.startswith(
            ("https://conda.anaconda.org/conda-forge/", "https://conda.anaconda.org/bioconda/")
        )
        for u in urls
    )


def test_installer_has_no_unsafe_patterns() -> None:
    text = (ROOT / "install.sh").read_text()
    code = "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("#"))
    assert not re.search(r"\|\s*(ba)?sh\b", code)  # no curl | bash
    # sudo may appear in advice text, but never in command position
    assert not re.search(r"(^|[;&|(]|\$\()\s*sudo\b", code, re.M)
    assert ".bashrc" not in code and ".profile" not in code.replace("~/.profile", "")
    assert re.search(r'MICROMAMBA_SHA256_LINUX_64="[0-9a-f]{64}"', text)


def test_shell_scripts_parse() -> None:
    for script in ("install.sh", "uninstall.sh", "scripts/update-lock.sh"):
        r = subprocess.run(
            ["bash", "-n", str(ROOT / script)], capture_output=True, text=True, check=False
        )
        assert r.returncode == 0, r.stderr


def test_installer_rejects_dangerous_prefixes(tmp_path: Path) -> None:
    for prefix in ("/", "//", ""):
        r = subprocess.run(
            [str(ROOT / "install.sh"), "--prefix", prefix],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            check=False,
        )
        assert r.returncode != 0 and "refusing" in r.stderr
    assert not any(tmp_path.iterdir())


def test_uninstaller_refuses_unmarked_directory(tmp_path: Path) -> None:
    (tmp_path / "keep.txt").write_text("x")
    r = subprocess.run(
        [str(ROOT / "uninstall.sh"), "--prefix", str(tmp_path), "--yes"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode != 0 and "marker" in r.stderr
    assert (tmp_path / "keep.txt").exists()


def test_uninstaller_removes_marked_install_and_only_own_launchers(tmp_path: Path) -> None:
    prefix, bindir = tmp_path / "pfx", tmp_path / "bin"
    (prefix / "bin").mkdir(parents=True)
    bindir.mkdir()
    (prefix / ".tpp-install-root").write_text("x")
    (prefix / "bin" / "transcriptomics-preprocess").write_text("#!/bin/sh\n# TPP-LAUNCHER\n")
    (bindir / "transcriptomics-preprocess").write_text("#!/bin/sh\necho someone else's\n")
    base = [str(ROOT / "uninstall.sh"), "--prefix", str(prefix), "--bin-dir", str(bindir)]
    r = subprocess.run([*base, "--dry-run"], capture_output=True, text=True, check=False)
    assert r.returncode == 0 and prefix.exists()
    r = subprocess.run(base, capture_output=True, text=True, stdin=subprocess.DEVNULL, check=False)
    assert r.returncode != 0 and prefix.exists()  # no tty, no --yes: refused
    r = subprocess.run([*base, "--yes"], capture_output=True, text=True, check=False)
    assert r.returncode == 0, r.stderr
    assert not prefix.exists()
    assert (bindir / "transcriptomics-preprocess").exists()  # foreign launcher kept


def test_gitignore_excludes_data() -> None:
    gi = (ROOT / ".gitignore").read_text()
    for pattern in ("*.sra", "*.fastq.gz", "*.bam", "tmp/", "__pycache__/"):
        assert pattern in gi


def test_python_m_entry_point() -> None:
    r = subprocess.run(
        [sys.executable, "-m", "transcriptomics_preprocessor", "--version"],
        capture_output=True,
        text=True,
        check=False,
        cwd=ROOT,
    )
    assert r.returncode == 0
