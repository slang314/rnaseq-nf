# CLAUDE.md

This repo is one of two pipelines described in `../CLAUDE.md` (the full handoff
brief: context, status, known risks, commands). Read that first if it exists.

When this repo stands alone on GitHub, the essentials are:
- nf-core-style layout: `main.nf` -> `workflows/` -> `modules/local/`; helper scripts in `bin/` (stdlib Python only).
- Test: `nextflow run main.nf -profile test,docker` (or `conda`), `nf-test test`, `python3 -m unittest discover -s tests -v`.
- Test data is synthetic with known answers (`scripts/simulate_test_data.py`); tests assert those answers.
- Use public data only. Never commit lab or unpublished data.
