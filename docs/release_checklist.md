# Release checklist

## Before tagging 0.1.0

- [ ] **Choose a license** (project owner) and add `LICENSE`; set `license` in
      `pyproject.toml`; update the README "License" section.
- [ ] Replace the placeholder author in `CITATION.cff` with real names (and ORCIDs).
- [ ] Run one real end-to-end job against NCBI on a small public run, paired and single:
      ```bash
      transcriptomics-preprocess run --accessions <small paired SRR> --layout paired --outdir real_pe --threads 4
      transcriptomics-preprocess run --accessions <small single SRR> --layout single --outdir real_se --threads 4
      ```
      Check `07_reports/summary_report.html`, `06_final_fastq/MANIFEST.csv` and
      `cd 06_final_fastq && sha256sum -c SHA256SUMS`, and record the result in
      `docs/validation_report.md`.
- [ ] Confirm that the GitHub Actions workflow passes on the default branch.
- [ ] Set the date in `CHANGELOG.md`, then tag: `git tag -a v0.1.0 -m "v0.1.0" && git push origin v0.1.0`.

## Publishing the repository (commands)

```bash
# one-time, for a new repository
git init -b main
git remote add origin https://github.com/<owner>/<repository>.git

# check exactly what will be committed (no data, environments or reports)
git status --short
git add -A
git ls-files                          # review the list
git ls-files | xargs du -k | sort -rn | head   # largest tracked files

# commit and push
git commit -m "Describe the change"
git push -u origin main
```

`.gitignore` excludes sequencing data (`*.sra`, `*.fastq*`, `*.bam`, …), pipeline output
directories, environments and caches. If `git status` shows any such file, do not commit it.
