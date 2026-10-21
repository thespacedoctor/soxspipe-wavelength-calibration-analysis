# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

An ad hoc data-analysis workspace, not an installable package, for the soxspipe sub-project "skyline wavelength calibration anchoring". The work checks the wavelength calibration of the SOXS spectrograph pipeline (`soxspipe`) against reference lines (HITRAN telluric lines, stellar absorption lines, sky emission lines) and archival X-Shooter spectra, and tests whether skyline fitting can anchor the wavelength solution.

There is no build, lint, or test suite. Work consists of editing the constants at the top of a script and re-running it, exploring data in the notebook, or re-reducing data with the `soxspipe` CLI.

## Environment

- Python dependencies (astropy, numpy, pandas, matplotlib, pygam) and the `soxspipe` CLI live in conda environments on this machine (e.g. `soxspipe`, `soxspipe_0.17.1`), not in this repo. The base `python3` lacks astropy and soxspipe, so activate an env first: `conda activate soxspipe`.
- `soxspipe` is a separate, independently installed project. This repo only uses its CLI and its output FITS products.
- `shared/telluric_lines_from_hitran.py` also needs `hitran-api` (`hapi`) and network access to hitran.org. The HITRAN API is rate-limited, so it caches tables in `shared/hitran_cache/`.

## Layout

- `input/`: source material for the project. Do not write into this folder.
- `output/`: deliverables and project-adjacent files.
- `soxs_uvvis_analysis/`: the UV-VIS notebook, reference tables, VIS spectra, results, reports, and two `soxspipe` VIS working directories. See `soxs_uvvis_analysis/README.md`.
- `soxs_nir_analysis/`: the NIR script `soxs_telluric_wavelength_offsets.py`, `rerun.sh`, the HITRAN line list, NIR spectra, X-Shooter products, and `legacy/` (a superseded script). See `soxs_nir_analysis/README.md`.
- `shared/`: `project_data/` (the 18 GB `soxspipe` working directory for both arms), `export.sh`, `cp_data.sh`, `telluric_lines_from_hitran.py`, and `hitran_cache/`. See `shared/README.md`.

The READMEs hold the method, the settings, and the data sets. Read them before you change code. Do not copy their content into this file.

## Common commands

All paths in the code are relative. Run each command from inside its own folder.

```bash
cd soxs_nir_analysis && python soxs_telluric_wavelength_offsets.py
cd soxs_uvvis_analysis && jupyter lab
cd shared && python telluric_lines_from_hitran.py --wmin 8000 --wmax 8400 --min-intensity 1e-25 --out telluric_8000_8400.py
```

The scripts take no command-line arguments (except `telluric_lines_from_hitran.py`). Set inputs by editing the `# INPUTS` block of the script or cell 2 of the notebook. Run `soxspipe` and `shared/export.sh` from inside `shared/project_data/`. `soxs_nir_analysis/rerun.sh` re-reduces one NIR exposure and runs the script again.

## Pipeline-version comparison

The same standard-star exposures are kept for three pipeline versions: `soxs_extractions_0.17.1` (tagged release), `soxs_extractions_v.dev` (dev branch), and `soxs_extractions_v.dev_horne_fix`. Each arm has its own copy in `soxs_nir_analysis/data/` and `soxs_uvvis_analysis/data/`. The UV-VIS data also include three `2026_june_july_standards_*` sets (attached arcs, unattached arcs, unattached arcs with sky line fitting). Switch between sets with `SOXS_ROOT` (NIR script) or `INDIR` (notebook).

## Known stale items

- `shared/cp_data.sh` has `DEST_DIR` set to an old Desktop path, and it reads from `/Volumes/SOXS/raw_archive/soxs_commissioning_data`, which must be mounted.
- The scripts in `soxs_nir_analysis/legacy/` still use the old paths.

## Version control

Git tracks code, the notebook, documents, and small tables only. `.gitignore` excludes `*/data/`, `shared/project_data/`, `shared/hitran_cache/`, `soxs_uvvis_analysis/soxspipe_workdirs/`, `soxs_uvvis_analysis/reference/*.tsv`, `soxs_uvvis_analysis/results/plots/`, `soxs_nir_analysis/legacy/2026-01/`, and the global patterns `*.fits`, `*.png`, `*.pdf`, `*.html`, and `*.log`.

## Agent skills

### Issue tracker

Issues live in Linear, team `thespacedoctor`, project `soxspipe`. See `docs/agents/issue-tracker.md`.

### Triage labels

Triage states are Linear workflow statuses (`Triage`, `needs-info`, `ready-for-agent`, `ready-for-dave`, `wontfix`); `bug` and `enhancement` are labels. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.
