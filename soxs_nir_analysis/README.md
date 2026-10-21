# SOXS NIR wavelength calibration analysis

This folder holds the script and data that measure the wavelength offsets of SOXS NIR spectra against HITRAN telluric lines. It is for people who run or change the analysis.

## Contents

- [What the script does](#what-the-script-does)
- [Folder layout](#folder-layout)
- [How to run the script](#how-to-run-the-script)
- [Settings](#settings)
- [Wavelength conversion](#wavelength-conversion)
- [Outputs](#outputs)
- [Data sets](#data-sets)
- [Compare pipeline versions](#compare-pipeline-versions)
- [Re-reduce one exposure with rerun.sh](#re-reduce-one-exposure-with-rerunsh)

## What the script does

`soxs_telluric_wavelength_offsets.py` compares the telluric absorption lines in SOXS NIR spectra with a list of HITRAN line positions. It also loads archival X-Shooter spectra of the same kind of target and shows them in the overview figure. The script does these steps:

1. Load the HITRAN line list and convert it from vacuum to air wavelengths (see [Wavelength conversion](#wavelength-conversion)). Choose the reference line: the line closest to the center of `SHIFT_REFERENCE_WINDOW`.
2. Load the SOXS spectra (`*_EXTRACTED_ME*.fits` in `SOXS_ROOT`). Keep a spectrum only if its `OBS NAME` (or `OBJECT`) header contains `SOXS_MATCH_STRING`. Skip files with `NOTFLAT` or `ERROR` in the name. Load the X-Shooter spectra in the same way with `XSHOOTER_MATCH_STRING`.
3. Divide the flux of each spectrum by its median (or mean, or maximum) in `NORM_WAVE_RANGE`.
4. Find one global wavelength shift for each SOXS spectrum. Fit a polynomial of degree `SHIFT_POLY_DEGREE` to the flux in `SHIFT_REFERENCE_WINDOW`. Take the wavelength of the minimum of the fit. The shift is the reference HITRAN wavelength minus that minimum.
5. For each HITRAN line in the wavelength range of the spectrum, find the flux minimum within `SEARCH_HALF_WINDOW_NM` of the line. Refine the minimum with a parabola through the three lowest points. The offset is the minimum minus the HITRAN wavelength. The script does this for the original and for the shifted wavelengths. It rejects an offset if it is larger than `MAX_OFFSET_TO_ACCEPT_NM`, if the window has fewer than `MIN_POINTS_IN_WINDOW` points, or (when `MIN_DEPTH` is set) if the line is too shallow.
6. Remove outliers. Compare each offset with a rolling median of 7 neighboring offsets. Apply sigma clipping to the differences (`clipped_good_mask`).
7. Plot the results: one page of panels for each group of spectra, and one summary plot of offset against wavelength.

The script prints the shift for each spectrum and the median and standard deviation of all shifts.

## Folder layout

| Path | Content |
| --- | --- |
| `soxs_telluric_wavelength_offsets.py` | The analysis script. |
| `rerun.sh` | Re-reduces one exposure with `soxspipe`, then runs the script. |
| `reference/telluric_5_lines.dat` | HITRAN telluric lines in vacuum nm (766 lines, 1097 nm to 2478 nm). One wavelength on each line. |
| `data/` | Spectra and plots. Not in version control. |
| `legacy/` | A superseded script and its old plots. |

## How to run the script

The script has no command-line arguments. Edit the constants in the `# INPUTS` block at the top of the file. Run it from inside this folder, because the paths are relative.

```bash
conda activate soxspipe
cd soxs_nir_analysis
python soxs_telluric_wavelength_offsets.py
```

The script needs `numpy`, `pandas`, `matplotlib`, and `astropy`. With `SHOW_FIG = True` it opens figure windows and saves the PNG files. Set `SHOW_FIG = False` to run without windows.

## Settings

All settings are at the top of the script.

| Constant | Current value | Meaning |
| --- | --- | --- |
| `SOXS_ROOT` | `data/soxs_extractions_v.dev_horne_fix` | Folder of SOXS spectra. A commented-out option is `data/tmp`. |
| `TELL_FILE` | `reference/telluric_5_lines.dat` | HITRAN line list (vacuum nm). |
| `XSHOOTER_ROOT` | `data/xshooter_extractions` | Folder of X-Shooter products. |
| `SOXS_MATCH_STRING` | `SOXS_LTT3218` | Text that the `OBS NAME` header must contain. `SOXS_Feige110` is a commented-out option. |
| `SOXS_FLUX_COL` | `FLUX_COUNTS` | Flux column. The script uses `extractedFluxOptimal` if this column is missing. |
| `XSHOOTER_MATCH_STRING` | `Hip012377_TELL_SLT` | Text that the X-Shooter header must contain. |
| `XSHOOTER_FILENAME` | `SDP_TELL_SLIT_FLUX_MERGE1D_NIR.fits` | File name to search for under `XSHOOTER_ROOT`. |
| `MAX_XSHOOTER_FILES` | `20` | Maximum number of X-Shooter spectra to load. |
| `NORMALIZE`, `NORM_MODE`, `NORM_WAVE_RANGE` | `True`, `"median"`, `(1671, 1676)` | Flux normalization (wavelengths in nm). |
| `TOPPANEL_WAVE_RANGE` | `(1700, 1850)` | Wavelength range of the top panel of the overview figure. |
| `SHIFT_REFERENCE_WINDOW` | `(1750.0, 1752.0)` | Window (nm) used to find the global shift. |
| `SHIFT_POLY_DEGREE` | `12` | Polynomial degree for the shift fit. |
| `SEARCH_HALF_WINDOW_NM` | `0.25` | Half-width (nm) of the search window around each line. |
| `MAX_OFFSET_TO_ACCEPT_NM` | `SEARCH_HALF_WINDOW_NM` | Largest accepted offset (nm). |
| `OFFSET_CLIP_NM` | `0.19` | Clipping limit for offsets (nm). |
| `MIN_POINTS_IN_WINDOW` | `5` | Fewest points for a line fit. |
| `PARABOLA_REFINE` | `True` | Refine each minimum with a parabola. |
| `MIN_DEPTH` | `None` | Smallest line depth. `None` turns the test off. |
| `N_ROWS`, `N_COLS` | `3`, `4` | Panels in each page of the grid plot. |
| `XLIM_NM`, `YLIM_OFFSET_NM` | `None`, `(-0.3, 0.3)` | Axis limits. |
| `SAVE_FIG`, `SHOW_FIG` | `True`, `True` | Save the PNG files. Show the figures. |
| `OUTDIR` | `SOXS_ROOT / "plots_soxs_telluric_offsets"` | Output folder. |
| `OUTFILE_PREFIX` | `SOXS_LTT3218_telluric_offsets_grid` | Start of each output file name. |

If you change `SOXS_MATCH_STRING`, also change `OUTFILE_PREFIX`. Otherwise the new plots overwrite the old plots.

## Wavelength conversion

The HITRAN line list is in vacuum wavelengths. The SOXS and X-Shooter spectra are in air wavelengths. The function `vac_to_air_um` converts a vacuum wavelength in micrometers to an air wavelength. It divides the wavelength by the refractive index `n`. The index uses the constants `a = 8.34213e-5`, `b1 = 2.406030e-2`, `c1 = 130.0`, `b2 = 1.5997e-4`, and `c2 = 38.9`, as written in the script.

`load_hitran_air_nm` reads the file, converts nm to micrometers, converts to air, converts back to nm, and removes duplicates. The script does this once, before it uses the lines.

## Outputs

The script writes PNG files to `<SOXS_ROOT>/plots_soxs_telluric_offsets/`:

- `<OUTFILE_PREFIX>_page01.png`, `_page02.png`, and so on. Each page is a grid of panels with the offset to HITRAN against HITRAN wavelength for each spectrum. Each panel shows the original and the shifted offsets.
- `<OUTFILE_PREFIX>_summary.png`. All offsets in one plot.

Each data folder in `data/` already holds the plots of an earlier run in this folder. A new run overwrites files that have the same name.

## Data sets

Each folder in `data/` holds NIR extracted, merged spectra (`*_NIR_*_EXTRACTED_MERGED.fits`) and the folder `plots_soxs_telluric_offsets/`.

| Folder | Content |
| --- | --- |
| `soxs_extractions_0.17.1/` | Spectra reduced with the tagged release 0.17.1 of `soxspipe`. |
| `soxs_extractions_v.dev/` | The same exposures, reduced with the development branch. |
| `soxs_extractions_v.dev_horne_fix/` | The same exposures, reduced with the development branch and a change to the Horne extraction. This set is the default `SOXS_ROOT`. |
| `xshooter_extractions/` | Archival ESO X-Shooter products. One `SDP_TELL_SLIT_FLUX_MERGE1D_NIR.fits` file in each folder named with a UUID. Files in a `_CorrOn` folder are skipped (`is_in_corron_path`). |
| `tmp/` | Working copies of single spectra, for example the output of `rerun.sh`. `tmp/hide/` holds intermediate `soxspipe` debug products for one exposure. |

The folders `soxs_extractions_0.17.1/` and `xshooter_extractions/` each hold a `list_spectra_parameters.sh` file. It lists the `DATE-OBS` header of each file with `fitsheader`. The VIS spectra of the three `soxs_extractions_*` sets are in `../soxs_uvvis_analysis/data/`.

## Compare pipeline versions

The three `soxs_extractions_*` folders hold the same standard-star exposures. Each folder has a different version of the pipeline. To compare versions:

1. Set `SOXS_ROOT` to the first folder. Run the script.
2. Set `SOXS_ROOT` to the next folder. Run the script.
3. Compare the median shift in the terminal and the plots in each `plots_soxs_telluric_offsets/` folder.

Each of the three folders holds 12 NIR spectra. The script keeps only the spectra that match `SOXS_MATCH_STRING`, so the count that it reports can differ between folders.

## Re-reduce one exposure with rerun.sh

`rerun.sh` re-reduces one exposure and runs the analysis again. The file names in the script are fixed. Edit them for a different exposure. The script does these steps:

1. Go to `../shared/project_data/`.
2. Run `soxspipe stare sessions/base/sof/<SOF_NAME>.sof -x` for the exposure `20260111T063325_NIR_3_STARE_STD_FLUX_SLIT5_0_150_0S_SOXS`.
3. Copy the `*_EXTRACTED_MERGED.fits` product from `sessions/base/reduced/2026-01-10/soxs-stare-std/` to `data/tmp/` in this folder.
4. Go back to this folder and run `soxs_telluric_wavelength_offsets.py`.

Run it from inside this folder:

```bash
conda activate soxspipe
bash rerun.sh
```

To analyze the new file, set `SOXS_ROOT = Path("data/tmp")` in the script first. To make the `sof` file for a raw exposure, see `../shared/README.md`.

## Legacy

`legacy/` holds an older script and its plots. They are superseded. Do not use them.
