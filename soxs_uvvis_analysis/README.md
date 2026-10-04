# SOXS UV-VIS wavelength calibration analysis

This folder holds the notebook (created by Geza, edited by Dave) and data that check the wavelength calibration of the SOXS UV-VIS arm. 

## Contents

- [What the notebook measures](#what-the-notebook-measures)
- [Method](#method)
- [Folder layout](#folder-layout)
- [How to run the notebook](#how-to-run-the-notebook)
- [Configuration](#configuration)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Current results](#current-results)
- [Caveats](#caveats)
- [Tests](#tests)
- [Standard stars and their spectral lines](#standard-stars-and-their-spectral-lines)
- [soxspipe working directories](#soxspipe-working-directories)
- [Known issues in the notebook](#known-issues-in-the-notebook)

## What the notebook measures

`STDCheck-Cleaned.ipynb` finds the wavelength correction that each spectral order (u, g, r, and i) of the SOXS UV-VIS arm needs. It uses spectra of flux standard stars. For each order, it gives:

- The **median and the mean shift** of that order. The correction to apply is <u>*minus*</u> the shift.
- The scatter of the shifts, in three forms:
  - **Internal scatter** is the line-to-line scatter inside one exposure, after the notebook removes that exposure's own median shift. It is empty (NaN) when each exposure has only one line.
  - **Exposure scatter** is the standard deviation of the per-exposure order medians.
  - **Star scatter** is the standard deviation of the per-star medians. 
- A robust version of the internal scatter, based on the median absolute deviation (MAD).
- The number of lines clipped. The notebook removes, in one pass, values that are more than 3 sigma from the median, where sigma comes from the MAD.

The statistics code is in `wavecal_stats.py`. The notebook imports it.

## Method

### Step 1. Measure the line shifts

1. Load every extracted spectrum (`*.fits`) in `INDIR`. Keep only spectra of the arm in `ARM` that were taken with the 1.0 arcsecond slit. Match each spectrum to a star name from its `HIERARCH ESO OBS TARG NAME` header.
2. For each line window, find the wavelength of the flux minimum. The notebook fits a generalized additive model (`pygam.LinearGAM`, 25 splines) to the flux in the window. It takes the minimum of the fitted curve on a grid of 2000 points (`fit_abs_line`).
3. Find the same minimum in a reference spectrum.
4. Record the shift: measured minimum minus reference minimum, in angstroms.

The notebook uses two kinds of line.

| Line type | Windows | Reference spectrum |
| --- | --- | --- |
| Stellar | The `[min, max]` windows in `STANDARD_STARS`. | The same star's rest-frame spectrum in `reference/flux_standards.fits`. The library wavelengths are air wavelengths. |
| Telluric | Each wavelength in `TELLURIC_LINES` (6277.70, 6868.47, 7252.95, and 8228.27 Å), plus and minus `TELLURIC_HALF_WINDOW_ANGSTROM` (7 Å). | The telluric model in `reference/soxs_vis_r5000_telluric_spectrum.tsv`. The notebook measures these lines in every exposure. |

Each window gets the order whose wavelength range contains the window midpoint (`assign_spectral_order`).

### Step 2. Barycentric correction

The Earth moves, and this motion Doppler-shifts the light of a star. The orbital motion is up to about 30 km/s. The rotation of the Earth adds about 0.46 km/s. At 5000 Å, a velocity of 30 km/s moves a line by about 0.5 Å. This shift is not an error of the wavelength solution, so the notebook must remove it.

The barycentric correction is the velocity that you add to the velocity seen from the telescope to get the velocity seen from the barycenter (the center of mass of the Solar System). `barycentric_correction_kms` calculates it for each exposure with `astropy`. It uses the `RA`, `DEC`, and `MJD-OBS` header keywords and the position of the NTT at La Silla. The result is the `BaryCorrKms` column.

Telluric lines form in the atmosphere of the Earth. They do not move with the Earth's motion relative to the star. For this reason, they give a fixed wavelength reference.

### Step 3. Fit the velocity of each star

1. For each exposure and each r or i order, the notebook calculates the median of the telluric shifts in that order, after 3 sigma MAD clipping (`telluric_order_offsets`). The notebook then compares the offsets of all exposures in that order and drops an exposure whose offset is more than 3 sigma from the others. This catches bad fits in the r order, which has only one telluric line in each exposure.
2. For each r and i stellar line in an exposure that has a telluric order offset, the notebook calculates a velocity:

   ```text
   vStar = c * (shift - telluric order offset) / wavelength + baryCorr
   ```

   Here `c` is the speed of light (299792.458 km/s), `shift` and `wavelength` are in angstroms, and `baryCorr` is the barycentric correction in km/s.
3. The velocity of the star is the median of these values after 3 sigma MAD clipping (`fit_star_velocities`). The velocity includes the radial velocity and the gravitational redshift of the star.

### Step 4. Calculate the residual

For each line, the notebook calculates the shift that the velocity of the star should cause, and removes it (`velocity_corrected_residuals`):

```text
residual = shift - wavelength * (vStar - baryCorr) / c
```

For a telluric line, the expected shift is 0, so the residual is the shift. The residual is the error of the wavelength solution.

### Step 5. Calculate the correction for each order

`clipped_order_stats` calculates the statistics of the residuals for each order. `ORDER_ANCHORS` sets which line type the notebook uses for each order:

| Order | Anchor | Reason |
| --- | --- | --- |
| u | Stellar | The order has no telluric lines. |
| g | Stellar | The order has no telluric lines. |
| r | Telluric | The order has one telluric line (6278 Å). |
| i | Telluric | The order has three telluric lines (6868, 7253, and 8228 Å). |

The r and i stellar lines fit the star velocity, so their residuals are close to 0 by construction. The notebook does not use them for the order corrections.

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
| `STANDARD_STARS` | Dictionary. Each key is a star name. Each value is a list of `[min, max]` stellar-line windows in angstroms (air). `BALMER_WINDOWS` is the shared list of Balmer-line windows. |
| `STANDARD_STAR_AKAS` | Alternative catalog names for stars (for example `LTT7987` is also `CD3017706`). |
| `TELLURIC_LINES_NAMED` | Dictionary of the telluric lines that the notebook measures in every exposure. Each key is a line name (for example `O2 B 6868`). Each value is its wavelength in angstroms (air): 6277.70, 6868.47, 7252.95, and 8228.27. |
| `TELLURIC_LINES` | List of the wavelengths in `TELLURIC_LINES_NAMED`. Do not edit it directly. |
| `TELLURIC_HALF_WINDOW_ANGSTROM` | Half-width of the window around each telluric line. The default is 7.0. |
| `ORDER_ANCHORS` | Dictionary that maps each order to the line type for its correction, `'stellar'` or `'telluric'`. The default is `u` and `g` stellar, `r` and `i` telluric. |
| `TELLURIC_WINDOWS` | Wide windows (angstroms) for the telluric sanity-check plots only. They do not affect the corrections. |
| `STELLAR_LINES_NAMED` | Dictionary of named stellar lines (for example `H epsilon`, `He II 4686`) in angstroms (air). It is used only to label the plots. |
| `REFERENCE_LINES_NAMED` | `STELLAR_LINES_NAMED` and `TELLURIC_LINES_NAMED` together. The comparison plot labels the line with the name of the named line in the window that is closest to the measured minimum. |
| `ORDERS_TO_INSPECT` | Spectral orders to keep, for example `['u', 'g', 'r', 'i']`. Use `None` to keep all. |
| `PLOT_METRICS` | Set to `True` to also write the `*_metrics.png` plots. The default is `False`. |
| `RESULTS_DIR`, `PLOTS_DIR` | Output folders. |
| `RESULTS_CSV_PATH`, `ORDER_CORRECTIONS_CSV_PATH`, `ORDER_CORRECTIONS_RAW_CSV_PATH`, `STAR_VELOCITIES_CSV_PATH`, `STAR_WINDOW_STATS_CSV_PATH` | Output CSV paths. |
| `SPECTRAL_FORMAT_CSV_PATH` | Per-order wavelength table for the current arm. |

The star names that `STANDARD_STARS` accepts come from the columns of `flux_standards.fits`: `GD71`, `LTT3218`, `GD153`, `EG274`, `LTT7987`, `Feige110`, and `EG21`. In the current notebook, the windows for `GD71`, `LTT3218`, and `GD153` are commented out.

## Inputs

### Reference tables (`reference/`)

| File | Use |
| --- | --- |
| `flux_standards.fits` | The `soxspipe` library of flux standard spectra. Column `WAVE` is in nm. One column per star. |
| `soxs_vis_r5000_telluric_spectrum.tsv` | Telluric model spectrum. Wavelength in angstroms, continuum-normalized flux. Not in version control. |
| `SOXS_SPECTRAL_FORMAT_TAB_VIS.csv`, `SOXS_SPECTRAL_FORMAT_TAB_NIR.csv` | Wavelength range of each spectral order. The notebook uses these to assign a window to an order (by the midpoint of the window). |
| `slit5_feige110_props_all_time.csv` | Header properties (airmass, seeing, exposure time, temperatures, rotator angle) of Feige 110 exposures with the 5.0 arcsecond slit. The current notebook does not read it. |
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
| `std_check_wavelength_shifts.csv` | One row for each spectrum, star, and line. Columns: `Star`, `LineType`, `Order`, `WindowMinAngstrom`, `WindowMaxAngstrom`, `Filename`, `LibraryMinAngstrom`, `LineWavelengthAngstrom`, `MeasuredMinAngstrom`, `ShiftAngstrom`, `BaryCorrKms`, `Airmass`, `MJD-OBS`, `Seeing`, `T_com`, `T_ccd`, `Rotang`, `ExpectedVelocityShiftAngstrom`, `ResidualAngstrom`. |
| `std_check_order_corrections.csv` | The main result. One row for each order, with the velocity-corrected residuals. Columns: `Order`, `Anchor`, `NLines`, `NExposures`, `MedianShiftAngstrom`, `MeanShiftAngstrom`, `StdShiftAngstrom`, `SEMeanAngstrom`, `InternalScatterStd`, `InternalScatterMadStd`, `ExposureScatterStd`, `NStars`, `StarScatterStd`, `NClipped`. |
| `std_check_order_corrections_raw.csv` | The same statistics for the raw shifts (`ShiftAngstrom`), with no velocity correction and all line types. It has no `Anchor` column. It is for comparison with the earlier method. |
| `std_check_star_velocities.csv` | The fitted velocity of each star. Columns: `Star`, `VStarKms`, `VStarErrKms`, `VStarMadStdKms`, `NLines`. |
| `std_check_stats_by_star_window.csv` | Statistics for each star and line. Columns: `Star`, `LineType`, `Order`, `WindowMinAngstrom`, `WindowMaxAngstrom`, `Count`, `MedianShiftAngstrom`, `MedianResidualAngstrom`, `StdResidualAngstrom`. |
| `plots/{star}_{min}-{max}A_comparison.png` | The SOXS spectra on top of the star's model spectrum in one window, with a histogram of the measured minima. Each SOXS spectrum is shifted into the star's rest frame by its expected velocity shift (fitted star velocity minus barycentric correction). A subtitle states the correction and the star velocity. The legend gives the median and standard deviation of the velocity-corrected shift (`ResidualAngstrom`). If the star has no fitted velocity, the plot shows the raw spectra and the raw shift, and the subtitle says so. |
| `plots/{star}_telluric_line_{wavelength}A_comparison.png` | The same plot for a telluric line, with the telluric model as the reference. The spectra are not shifted, and a box gives the median and standard deviation of the raw shift. |
| `plots/star_velocity_diagnostic.png` | The velocity of each stellar line against wavelength, one panel for each star. The dashed line is the fitted star velocity. A pure velocity offset gives a flat line. An error of the wavelength solution in one order gives a step. |
| `plots/{star}_{min}-{max}A_metrics.png` | Shift against airmass, MJD, seeing, temperatures, and rotator angle. Written only when `PLOT_METRICS = True`. |
| `plots/{star}_{min}-{max}A_telluric.png` | Telluric sanity check. The SOXS spectra on top of the telluric model in a wide window. |

The CSV files and plots are overwritten each time you run the notebook. With the default `PLOT_METRICS = False`, the existing `*_metrics.png` files are not updated. Statistics and plots use only the orders in `ORDERS_TO_INSPECT`.

## Current results

These values come from the CSV files in `results/`. They are for the default `INDIR` (`2026_june_july_standards_unattached_arcs_w_skyline_fitting/`), 19 exposures of four stars. Run the notebook again to update them.

### Star velocities

The velocities are barycentric, in km/s. `Error` is the MAD standard deviation divided by the square root of the number of lines.

| Star | Velocity (km/s) | Error (km/s) | MAD standard deviation (km/s) | Lines |
| --- | --- | --- | --- | --- |
| EG21 | 44.87 | 4.68 | 8.11 | 3 |
| EG274 | 30.09 | 2.05 | 5.02 | 6 |
| Feige110 | 14.16 | 3.93 | 10.40 | 7 |
| LTT7987 | 60.20 | 7.41 | 14.82 | 4 |

The error is a lower limit, because lines in one exposure share that exposure's offset.

### Corrections by order

All values are in angstroms, from `std_check_order_corrections.csv`. The `Clipped` column is a number of lines.

| Order | Anchor | Lines | Exposures | Median | Mean | SE mean | Internal std | Internal MAD std | Exposure std | Star std | Clipped |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| u | Stellar | 38 | 19 | 0.204 | 0.180 | 0.102 | 0.259 | 0.275 | 0.442 | 0.224 | 0 |
| g | Stellar | 41 | 19 | 0.140 | 0.014 | 0.149 | 0.208 | 0.172 | 0.649 | 0.411 | 0 |
| r | Telluric | 15 | 15 | 0.104 | 0.079 | 0.060 | NaN | NaN | 0.232 | 0.259 | 4 |
| i | Telluric | 55 | 19 | -0.013 | -0.032 | 0.054 | 0.203 | 0.067 | 0.234 | 0.150 | 2 |

Before the velocity correction, the u and g medians were 0.740 Å and 0.510 Å (`std_check_order_corrections_raw.csv`). Most of that was the velocity of the stars.

## Caveats

- **The u and g corrections depend on the fitted star velocities.** The velocity comes mostly from Hα in the r order. For Feige110, it also comes from the He lines. An error `σv` in the velocity gives an error of about `λ·σv/c` in the correction. For example, 4 km/s at 4000 Å is about 0.05 Å.
- **The r order has one telluric line (6278 Å) in each exposure.** Its internal scatter cannot be measured, so it is NaN. The notebook clipped 4 of 19 r exposures as bad fits, so the r row uses 15 lines.
- **The u and g results still depend on the star.** After the velocity correction, the median u residual is 0.03 to 0.50 Å and the median g residual is −0.26 to 0.64 Å, depending on the star. The star scatter (0.22 Å in u, 0.41 Å in g) is larger than the SE of the mean. Possible causes are: a velocity from Hα that does not fit the higher Balmer lines, too few lines per star, or a real error between the r and g orders. This data cannot separate them. Treat the u and g corrections as uncertain to about the star scatter.
- **The r and i stellar residuals are about 0 by construction.** The notebook uses these lines to fit the star velocity. It does not use them for the order corrections.
- **The exposure-to-exposure scatter in u and g is large.** It is about 0.4 to 0.65 Å, and the internal scatter is about 0.2 to 0.25 Å. Most of the scatter is from one exposure to the next, not from line to line. The cause is not known.
- **The header keywords `ESO QC SKY SHIFT O2` and `ESO QC SKY SHIFT O3` are exactly 0.0 in all files.** The reason is not known. Check it.
- **The star velocities have few lines.** Each star has 3 to 7 lines. The errors in the table above show the effect.

## Tests

The statistics code in `wavecal_stats.py` has unit tests in `tests/`. Run them from this folder:

```bash
cd soxs_uvvis_analysis
python -m pytest tests -q
```

The tests need `pytest`, `numpy`, `pandas`, and `astropy`.

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
