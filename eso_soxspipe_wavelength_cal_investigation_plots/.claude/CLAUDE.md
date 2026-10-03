# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

An ad hoc data-analysis workspace (not an installable package) for investigating the wavelength
calibration accuracy of the SOXS spectrograph pipeline (`soxspipe`). It compares SOXS-extracted
flux-standard spectra against a HITRAN telluric absorption-line list and against archival
X-Shooter spectra of the same standard stars, to quantify wavelength offsets/drift between SOXS
pipeline versions and against ESO's X-Shooter.

There is no build, lint, or test suite. Work here consists of editing the constants at the top of
a script and re-running it, or exploring data interactively in the notebook.

## Environment

- Python dependencies (astropy, numpy, pandas, matplotlib, pygam) and the `soxspipe` CLI live in
  conda environments on this machine (e.g. `soxspipe`, `soxspipe_0.17.1`), not in this repo. The
  base `python3` on PATH does not have astropy/soxspipe installed — activate the relevant conda
  env first, e.g. `conda activate soxspipe`.
- `soxspipe` itself (the reduction pipeline under test) is a separate, independently installed
  project. This repo only consumes its CLI and its output FITS products.

## Running the analysis scripts

Both analysis scripts are configured by editing the constants in their `# INPUTS` block at the top
of the file (input directories, target star match string, wavelength windows, fit/clipping
parameters, output paths) — there is no CLI argument parsing. After editing, just run:

```bash
python soxs_telluric_wavelength_offsets.py
python soxs_plot_spectrum_tell_vs_xshooter_shift.py
```

- `soxs_telluric_wavelength_offsets.py` — current/primary script. Measures SOXS emission/absorption
  line-minimum wavelengths against a HITRAN telluric line list (`telluric_5_lines.dat`), derives a
  global per-spectrum wavelength shift from a reference line, and plots per-spectrum and summary
  offset-vs-wavelength grids to `<SOXS_ROOT>/plots_soxs_telluric_offsets/`.
- `soxs_plot_spectrum_tell_vs_xshooter_shift.py` — earlier/related version of the same analysis
  that also overplots archival X-Shooter spectra for a three-way (SOXS vs HITRAN vs X-Shooter)
  comparison. Output goes to `<SOXS_ROOT>/plots_soxs_vs_xshooter/`.
- `STDCheck-Cleaned.ipynb` — interactive notebook: loads SOXS-extracted spectra plus the
  `flux_standards.fits` reference library, fits absorption-line minima with a GAM
  (`pygam.LinearGAM`), and looks at line-minimum position vs. airmass/other header properties.

Reprocessing a single raw exposure end-to-end (reduce with soxspipe, copy the extracted product
into `tmp/`, then re-run the offset analysis) follows the pattern in `rerun.sh`.

## Directory layout

- `project_data/` — a `soxspipe` working directory (raw frames, `soxspipe.yaml` recipe config,
  `sof`/session files, calibration/science reductions, `soxspipe.db`, logs). This is where
  `soxspipe` itself is actually invoked from.
- `soxs_extractions_0.17.1/` vs `soxs_extractions_v.dev/` — the same set of extracted, merged SOXS
  standard-star spectra reduced with two different pipeline versions (a tagged release vs. the dev
  branch), kept side by side so offset results can be compared across pipeline versions. Each has
  its own `plots_soxs_telluric_offsets/` output.
- `xshooter_extractions/` — archival ESO X-Shooter `SDP_TELL_SLIT_FLUX_MERGE1D_NIR.fits` products
  (one per UUID-named directory), used as the independent-instrument reference.
  `is_in_corron_path()`/`_CorrOn` filtering in the scripts excludes atmospheric-correction-on
  variants.
- `tmp/` and `2026-01/` — scratch/working copies of extracted spectra and plots used while
  iterating on a specific reduction (e.g. `tmp/hide/` holds intermediate soxspipe debug products
  such as unflattened/untuned flux and sky-subtraction FITS for one exposure).
- `telluric_5_lines.dat` — HITRAN reference line list in **vacuum nm**; scripts convert it to air
  wavelengths via `vac_to_air_um()` before comparing to SOXS/X-Shooter (both in air wavelengths).
- `flux_standards.fits` — soxspipe's library flux-standard reference table (used by the notebook).
- `SOXS_SPECTRAL_FORMAT_TAB_{NIR,VIS}.csv` — per-order wavelength coverage tables for the SOXS
  spectral format.
- `slit5_feige110_props_all_time*.csv` — precomputed per-exposure header properties (airmass,
  seeing, exptime, temperatures, rotator angle) for Feige 110 SLIT5 observations, used by the
  notebook instead of re-reading every FITS header.
- `cp_data.sh` — pulls specific raw commissioning FITS frames (by exact filename) from the
  `/Volumes/SOXS/raw_archive/soxs_commissioning_data` archive into `soxs_extractions_v.dev/`.
- `export.sh` — generates `soxspipe raw sof` files for a fixed list of raw exposures.

## Analysis method (shared by both offset scripts)

1. Load target spectra: a SOXS extracted-merged FITS table is kept only if its `OBS NAME`/`OBJECT`
   header matches `SOXS_MATCH_STRING` (e.g. `SOXS_Feige110`); X-Shooter files are matched the same
   way via `XSHOOTER_MATCH_STRING`.
2. Flux is normalized (median, by default) over `NORM_WAVE_RANGE` so spectra are comparable.
3. A single global wavelength shift per SOXS spectrum is derived by fitting a high-order polynomial
   to a strong reference feature in `SHIFT_REFERENCE_WINDOW` and comparing its fitted minimum to
   the nearest HITRAN line (`measure_reference_minimum_for_shift` / `apply_soxs_shift`).
4. For every HITRAN line inside the spectrum's coverage, a local parabola fit around the flux
   minimum gives a per-line offset (`find_local_minimum_near_line`); results are sigma-clipped
   against a rolling local median (`clipped_good_mask`) before plotting/statistics.
5. Output is a paginated grid of per-spectrum offset plots plus one summary
   offset-vs-wavelength scatter plot, saved as PNGs.

When changing analysis behavior, prefer editing the named constants in the `# INPUTS` block over
introducing CLI flags — that matches how both scripts are currently structured and how they're run
interactively/iteratively during an investigation.


## Standards Spectral Info

Here's a detailed summary of these stars with their spectral types and prominent absorption lines:

| Star | Spectral Type | Primary Absorption Lines | Wavelengths |
|------|---------------|-------------------------|-------------|
| **GD71** | DA1 | Hydrogen Balmer series | H-α: 656.3 nm; H-β: 486.1 nm; H-γ: 434.0 nm; H-δ: 410.2 nm |
| **LTT3218** | DA5 | Hydrogen Balmer series | H-α: 656.3 nm; H-β: 486.1 nm; H-γ: 434.0 nm; H-δ: 410.2 nm |
| **GD153** | DA1.5 | Hydrogen Balmer series | H-α: 656.3 nm; H-β: 486.1 nm; H-γ: 434.0 nm; H-δ: 410.2 nm |
| **EG274** | DA2 | Hydrogen Balmer series | H-α: 656.3 nm; H-β: 486.1 nm; H-γ: 434.0 nm; H-δ: 410.2 nm |
| **LTT7987** | DA | Hydrogen Balmer series | H-α: 656.3 nm; H-β: 486.1 nm; H-γ: 434.0 nm; H-δ: 410.2 nm |
| **Feige110** | DOp (sdO) | Balmer lines and weak HeII 4686Å | H-α: 656.3 nm; He II: 468.6 nm |
| **EG21** | DA | Hydrogen Balmer series | H-α: 656.3 nm; H-β: 486.1 nm; H-γ: 434.0 nm; H-δ: 410.2 nm |

## Key Notes:

- **DA white dwarfs** (GD71, LTT3218, GD153, EG274, LTT7987, EG21) are hydrogen-rich white dwarfs where the most accurate determinations are made by fitting theoretical spectra to absorption lines in the Balmer series of the optical range. The numbers (1, 5, 1.5, 2) indicate effective temperature subclasses.

- **Feige110 (sdO)**: This is a subdwarf O star—a much hotter object than the DA white dwarfs. The spectra of sdO stars are dominated by weak Balmer and strong He II lines depending on the He abundance.

All the DA stars show very broad, gravity-broadened hydrogen Balmer lines that are characteristic of their high surface gravity environments.
