# SOXS UV-VIS wavelength calibration analysis

This folder holds the notebook and data that check the wavelength calibration of the SOXS UV-VIS arm. It is for people who run or change the analysis.

## Contents

- [What the notebook measures](#what-the-notebook-measures)
- [Folder layout](#folder-layout)
- [How to run the notebook](#how-to-run-the-notebook)
- [Configuration](#configuration)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Standard stars and their spectral lines](#standard-stars-and-their-spectral-lines)
- [soxspipe working directories](#soxspipe-working-directories)
- [Known issues in the notebook](#known-issues-in-the-notebook)

## What the notebook measures

`STDCheck-Cleaned.ipynb` measures how far the wavelength scale of SOXS spectra is from a reference. It uses spectra of flux standard stars. It does the following:

1. Load every extracted spectrum (`*.fits`) in `INDIR`. Keep only spectra of the arm in `ARM` that were taken with the 1.0 arcsecond slit. Match each spectrum to a star name from its `HIERARCH ESO OBS TARG NAME` header.
2. For each star and each wavelength window in `STANDARD_STARS`, find the wavelength of the flux minimum. The notebook fits a generalized additive model (`pygam.LinearGAM`, 25 splines) to the flux in the window. It takes the minimum of the fitted curve on a grid of 2000 points (`fit_abs_line`).
3. Find the same minimum in a reference spectrum. For a stellar line, the reference is the star's model spectrum in `reference/flux_standards.fits`. For a telluric line, the reference is the telluric model spectrum in `reference/soxs_vis_r5000_telluric_spectrum.tsv`.
4. Record the shift: measured minimum minus reference minimum, in angstroms.
5. Write the shifts to a CSV file. Then calculate statistics (count, median, mean, standard deviation, minimum, maximum) for the whole detector, for each spectral order, and for each star and window.
6. Make plots for each star and window.

A window with equal minimum and maximum values (for example `[6563, 6563]`) is a telluric line. The notebook widens that window by 7 angstroms on each side. It then uses the telluric model as the reference.

## Folder layout

| Path | Content |
| --- | --- |
| `STDCheck-Cleaned.ipynb` | The analysis notebook. |
| `reference/` | Reference tables that the notebook reads. |
| `data/` | Extracted spectra (inputs). Not in version control. |
| `results/` | CSV files and `plots/` (outputs). |
| `reports/` | Exports of earlier versions of the notebook (HTML and PDF). |
| `soxspipe_workdirs/` | Two complete `soxspipe` working directories and a DS9 region file. |

## How to run the notebook

All paths in the notebook are relative. Start Jupyter from inside this folder.

```bash
conda activate soxspipe
cd soxs_uvvis_analysis
jupyter lab
```

Open `STDCheck-Cleaned.ipynb`, change the constants in cell 2 (see [Configuration](#configuration)), and run all cells. The notebook needs `pandas`, `numpy`, `matplotlib`, `astropy`, `pygam`, `tabulate` (for `to_markdown`), and `IPython`. These are not installed by this repository.

## Configuration

Cell 2 of the notebook holds all settings.

| Constant | Meaning |
| --- | --- |
| `INDIR` | Folder of extracted spectra. Three options are in the cell; two are commented out. |
| `ARM` | `'VIS'` or `'NIR'`. Selects the spectra to keep and the file `reference/SOXS_SPECTRAL_FORMAT_TAB_{ARM}.csv`. |
| `TELLURIC_SPECTRUM` | Path to the telluric model spectrum. |
| `STANDARD_STARS` | Dictionary. Each key is a star name. Each value is a list of `[min, max]` windows in angstroms. |
| `STARDARD_STAR_AKAS` | Alternative catalog names for stars (for example `LTT7987` is also `CD3017706`). The name in the code has this spelling. |
| `TELLURIC_WINDOWS` | Wide windows (angstroms) for the telluric sanity-check plots. |
| `TELLURIC_LINES` | Telluric line wavelengths (angstroms) that the telluric plots mark. |
| `TELLURIC_LINES_NAMED` | Named lines (for example `Ca II K`, `H alpha`). The plot uses the name as the label when a window contains the line. |
| `ORDERS_TO_INSPECT` | Spectral orders to keep, for example `['u', 'g', 'r', 'i']`. Use `None` to keep all. |
| `PLOT_METRICS` | Set to `True` to also write the `*_metrics.png` plots. The default is `False`. |
| `RESULTS_DIR`, `PLOTS_DIR` | Output folders. |
| `RESULTS_CSV_PATH`, `STAR_WINDOW_STATS_CSV_PATH`, `ORDER_STATS_CSV_PATH`, `DETECTOR_STATS_CSV_PATH` | Output CSV paths. |
| `SPECTRAL_FORMAT_CSV_PATH` | Per-order wavelength table for the current arm. |

The star names that `STANDARD_STARS` accepts come from the columns of `flux_standards.fits`: `GD71`, `LTT3218`, `GD153`, `EG274`, `LTT7987`, `Feige110`, and `EG21`. In the current notebook, the windows for `GD71`, `LTT3218`, and `GD153` are commented out.

## Inputs

### Reference tables (`reference/`)

| File | Use |
| --- | --- |
| `flux_standards.fits` | The `soxspipe` library of flux standard spectra. Column `WAVE` is in nm. One column per star. |
| `soxs_vis_r5000_telluric_spectrum.tsv` | Telluric model spectrum. Wavelength in angstroms, continuum-normalized flux. Not in version control. |
| `SOXS_SPECTRAL_FORMAT_TAB_VIS.csv`, `SOXS_SPECTRAL_FORMAT_TAB_NIR.csv` | Wavelength range of each spectral order. The notebook uses these to assign a window to an order (by the midpoint of the window). |
| `slit5_feige110_props_all_time.csv` | Header properties (airmass, seeing, exposure time, temperatures, rotator angle) of Feige 110 exposures with the 5.0 arcsecond slit. The notebook reads it in cell 4. |
| `slit5_feige110_props_all_time_orig.csv` | Earlier copy of the table above. No code reads it. |

### Spectra (`data/`)

Each folder holds extracted, merged spectra (`*_EXTRACTED_MERGED.fits`). The notebook reads the columns `WAVE` (nm) and `FLUX_COUNTS`, and header keywords such as airmass, seeing, and temperatures.

| Folder | Content |
| --- | --- |
| `2026_june_july_standards_attached_arcs/` | 28 VIS spectra. Reduced with attached arcs. |
| `2026_june_july_standards_unattached_arcs/` | The same 28 spectra. Reduced with unattached arcs. |
| `2026_june_july_standards_unattached_arcs_w_skyline_fitting/` | The same 28 spectra. Reduced with unattached arcs and sky line fitting. This set is the default `INDIR`. |
| `soxs_extractions_0.17.1/`, `soxs_extractions_v.dev/`, `soxs_extractions_v.dev_horne_fix/` | The VIS spectra of three versions of the pipeline: the tagged release 0.17.1, the development branch, and the development branch with a change to the Horne extraction. The NIR spectra of the same sets are in `../soxs_nir_analysis/data/`. |
| `tmp/` | One VIS spectrum (`20260628T053000...`). A scratch copy. |

The notebook skips spectra of other arms and spectra that were not taken with the 1.0 arcsecond slit. It therefore uses fewer spectra than a folder holds.

## Outputs

The notebook writes all outputs to `results/`.

| File | Content |
| --- | --- |
| `std_check_wavelength_shifts.csv` | One row for each spectrum, star, and window. Columns: `Star`, `WindowMinAngstrom`, `WindowMaxAngstrom`, `Filename`, `LibraryMinAngstrom`, `MeasuredMinAngstrom`, `ShiftAngstrom`, `Airmass`, `MJD-OBS`, `Seeing`, `T_com`, `T_ccd`, `Rotang`. |
| `std_check_stats_by_detector.csv` | Shift statistics for the whole arm. |
| `std_check_stats_by_order.csv` | Shift statistics for each spectral order. |
| `std_check_stats_by_star_window.csv` | Shift statistics for each star, order, and window. |
| `plots/{star}_{min}-{max}A_comparison.png` | The SOXS spectra on top of the star's model spectrum in one window, with a histogram of the measured minima. |
| `plots/{star}_telluric_line_{wavelength}A_comparison.png` | The same plot for a telluric line, with the telluric model as the reference. |
| `plots/{star}_{min}-{max}A_metrics.png` | Shift against airmass, MJD, seeing, temperatures, and rotator angle. Written only when `PLOT_METRICS = True`. |
| `plots/{star}_{min}-{max}A_telluric.png` | Telluric sanity check. The SOXS spectra on top of the telluric model in a wide window. |

The CSV files and plots are overwritten each time you run the notebook. With the default `PLOT_METRICS = False`, the existing `*_metrics.png` files are not updated. Statistics and plots use only the orders in `ORDERS_TO_INSPECT`.

The folder `reports/` holds HTML and PDF exports of earlier notebook runs. Their titles are "SOXS VIS Wavelength Calibration Investigation - ONE" and "- TWO". The current notebook header says "THREE".

## Standard stars and their spectral lines

This table comes from the earlier guide for this workspace. The spectral types and line wavelengths are not checked against the code or a catalog. Confirm them before you use them in a publication.

| Star | Spectral type | Main absorption lines | Wavelengths |
| --- | --- | --- | --- |
| GD71 | DA1 | Hydrogen Balmer series | H-alpha 656.3 nm; H-beta 486.1 nm; H-gamma 434.0 nm; H-delta 410.2 nm |
| LTT3218 | DA5 | Hydrogen Balmer series | Same as GD71 |
| GD153 | DA1.5 | Hydrogen Balmer series | Same as GD71 |
| EG274 | DA2 | Hydrogen Balmer series | Same as GD71 |
| LTT7987 | DA | Hydrogen Balmer series | Same as GD71 |
| Feige110 | DOp (sdO) | Balmer lines and He II | H-alpha 656.3 nm; He II 468.6 nm |
| EG21 | DA | Hydrogen Balmer series | Same as GD71 |

- DA stars are white dwarfs with hydrogen-rich atmospheres. Their Balmer lines are very broad.
- Feige110 is a subdwarf O star. It is hotter than the DA stars.

## soxspipe working directories

`soxspipe_workdirs/` holds two complete `soxspipe` working directories for VIS data: `2026_june_VIS_example/` and `VIS_attached_arc_example/`. Each one has `raw/`, `sessions/base/`, `soxspipe.db`, and `soxspipe.yaml`. The folders `sof`, `reduced`, `qc` and the files `soxspipe.yaml` and `soxspipe.log` are symlinks into `sessions/base/`. Run `soxspipe` from inside the directory.

The two `soxspipe.yaml` files differ in two places only:

- `soxs-stare-std` → `sky-subtraction` → `subtract_sky` is `True` in `2026_june_VIS_example` and `False` in `VIS_attached_arc_example`.
- The log file path (`filename`).

The file `20260707T071228_soxs_absorption_lines.reg` is a DS9 region file in physical pixel coordinates. It labels stellar absorption lines (the Balmer series and Ca II H) on the exposure `20260707T071228`.

The folder `soxspipe_workdirs/` is not in version control.

## Known issues in the notebook

- Cell 2 has a `TODO` for the stars `GD71`, `LTT3218`, and `GD153`. These stars have no windows defined.
