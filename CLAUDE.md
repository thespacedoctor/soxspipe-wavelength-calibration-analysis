# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

An ad hoc data-analysis workspace, not an installable package, for the soxspipe sub-project "skyline wavelength calibration anchoring". The work checks the wavelength calibration of the SOXS spectrograph pipeline (`soxspipe`) against reference lines (HITRAN telluric lines, stellar absorption lines, sky emission lines) and archival X-Shooter spectra, and tests whether skyline fitting can anchor the wavelength solution.

There is no build, lint, or test suite. Work consists of editing the constants at the top of a script and re-running it, exploring data in the notebook, or re-reducing data with the `soxspipe` CLI.

## Environment

- Python dependencies (astropy, numpy, pandas, matplotlib, pygam) and the `soxspipe` CLI live in conda environments on this machine (e.g. `soxspipe`, `soxspipe_0.17.1`), not in this repo. The base `python3` lacks astropy and soxspipe, so activate an env first: `conda activate soxspipe`.
- `soxspipe` is a separate, independently installed project. This repo only uses its CLI and its output FITS products.
- `telluric_lines_from_hitran.py` also needs `hitran-api` (`hapi`) and network access to hitran.org. The HITRAN API is rate-limited, so it caches tables in `hitran_cache/`.

## Layout

- `input/`: source material for the project. Do not write into this folder.
- `output/`: deliverables and project-adjacent files.
- `eso_soxspipe_wavelength_cal_investigation_plots/`: the analysis workspace (scripts, notebook, reference tables, extracted spectra, plots). It has its own detailed guide at `eso_soxspipe_wavelength_cal_investigation_plots/.claude/CLAUDE.md` covering each script, the analysis method, and the standard stars' spectral lines. Read it before working there.
- `2026_june_VIS_example/` and `VIS_attached_arc_example/`: two self-contained `soxspipe` VIS working directories (`raw/`, `sessions/base/`, `soxspipe.db`, with `sof`, `reduced`, `qc`, `soxspipe.yaml` and `soxspipe.log` symlinked into `sessions/base/`). Run `soxspipe` from inside the directory. Their `soxspipe.yaml` files differ only in `soxs-stare-std` → `sky-subtraction` → `subtract_sky` (`True` in the June example, `False` in the attached-arc example) and in the log path.
- `20260707T071228_soxs_absorption_lines.reg`: a DS9 region file (physical pixel coordinates) labelling stellar absorption lines (Balmer series, Ca II H) on the `20260707T071228` exposure.

## Common commands

Run from `eso_soxspipe_wavelength_cal_investigation_plots/`. The analysis scripts take no CLI arguments; set inputs by editing the `# INPUTS` block at the top of each file.

```bash
python soxs_telluric_wavelength_offsets.py
python soxs_plot_spectrum_tell_vs_xshooter_shift.py
python telluric_lines_from_hitran.py --wmin 8000 --wmax 8400 --min-intensity 1e-25 --out telluric_8000_8400.py
```

`telluric_lines_from_hitran.py` takes Angstroms and writes air wavelengths by default (`--vacuum` switches to vacuum). It writes the `--out` file plus a `.txt` copy.

Re-reduce one exposure and re-run the offsets (the pattern in `rerun.sh`):

```bash
cd project_data/
soxspipe stare sessions/base/sof/<SOF_NAME>.sof -x
cp sessions/base/reduced/<DATE>/soxs-stare-std/<NAME>_EXTRACTED_MERGED.fits ../tmp/
cd ..
python soxs_telluric_wavelength_offsets.py
```

Generate a `sof` file for a raw exposure (as in `export.sh`): `soxspipe raw sof <SOF_NAME>.sof`.

## Pipeline-version comparison

`extraction_archive/` holds the same standard-star spectra reduced in different ways, side by side, so offset results can be compared: `soxs_nir_extractions_0.17.1` (tagged release), `soxs_nir_extractions_v.dev` (dev branch), `soxs_nir_extractions_v.dev_horne_fix`, and the 2026 June–July VIS standards with attached arcs, unattached arcs, and unattached arcs with skyline fitting. Archival X-Shooter products are in `extraction_archive/xshooter_extractions/`. Switch between sets by changing `SOXS_ROOT` in the script's `# INPUTS` block.

## Known stale paths

The subfolder guide and the scripts predate the move of the extraction folders into `extraction_archive/`:

- Both scripts set `XSHOOTER_ROOT = Path("xshooter_extractions")`, but the folder is now `extraction_archive/xshooter_extractions/`.
- `soxs_telluric_wavelength_offsets.py` assigns `SOXS_ROOT` twice; the second assignment wins.
- `cp_data.sh` copies into an old Desktop path (`/Users/Dave/Dropbox/Mac/Desktop/...`) and pulls from `/Volumes/SOXS/raw_archive/soxs_commissioning_data`, which must be mounted.
- `rerun.sh` copies into `tmp/`, which no longer exists at the workspace root (it is now `extraction_archive/tmp/`).

## Version control

Git tracks code, the notebook, guides, and small tables only. It ignores FITS files, generated plots and reports (`*.png`, `*.pdf`, `*.html`), logs, the two VIS working directories, `*.reg`, and the data folders inside the analysis workspace (`2026-01/`, `extraction_archive/`, `project_data/`, `hitran_cache/`, `tmp/`, `results/plots/`, `*.tsv`). See `.gitignore`.

## Agent skills

### Issue tracker

Issues live in Linear, team `thespacedoctor`, project `soxspipe`. See `docs/agents/issue-tracker.md`.

### Triage labels

Triage states are Linear workflow statuses (`Triage`, `needs-info`, `ready-for-agent`, `ready-for-dave`, `wontfix`); `bug` and `enhancement` are labels. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: one `CONTEXT.md` and `docs/adr/` at the repo root. See `docs/agents/domain.md`.
