# SOXS NIR wavelength calibration analysis

This folder contains the script and data used to measure wavelength offsets in SOXS NIR spectra against HITRAN telluric lines. The script was originally written by Michaël Marsset, and Dave has edited it.

## What the script does

`soxs_telluric_wavelength_offsets.py` compares the telluric absorption lines in SOXS NIR spectra with a list of HITRAN line positions. It also loads archival X-Shooter spectra of the same kind of target. It uses them to check the HITRAN line list and to show an overview figure. The script does these steps:

1. Load the HITRAN line list and convert it from vacuum to air wavelengths (see [Wavelength conversion](#wavelength-conversion)).
2. Load the SOXS spectra (`*_EXTRACTED_ME*.fits` in `SOXS_ROOT`). Keep a spectrum only if its `OBS NAME` (or `OBJECT`) header contains `SOXS_MATCH_STRING`. Skip files with `NOTFLAT` or `ERROR` in the name. Load the X-Shooter spectra the same way, using `XSHOOTER_MATCH_STRING`.
3. Divide each spectrum's flux by its median in `NORM_WAVE_RANGE`.
4. Measure each X-Shooter spectrum against all HITRAN lines (see step 5). Keep a HITRAN line only if it is accepted, not clipped, and has an offset smaller than `XSH_MAX_OFFSET_NM` in at least a fraction `MIN_XSH_FRACTION` of the X-Shooter spectra that cover it. This needs at least `MIN_XSH_SPECTRA` X-Shooter spectra. If there are fewer, the script prints a warning, measures all HITRAN lines, and does not calculate the SOXS − X-Shooter differential.
5. For each kept line in the wavelength range of a spectrum, find the flux minimum within `SEARCH_HALF_WINDOW_NM` of the line. Refine an interior minimum with a parabola through the three points around it. The offset is the minimum minus the HITRAN wavelength. The script rejects a line and gives it one reason. It tests in this order:
   - `too_few_points`: the window has fewer than `MIN_POINTS_IN_WINDOW` points.
   - `edge_minimum`: the lowest point is the first or last point of the window.
   - `too_shallow`: the depth is less than `MIN_DEPTH` (when `MIN_DEPTH` is set). The depth is the `CONTINUUM_PERCENTILE` percentile of the window flux minus the lowest flux.
   - `saturated`: the lowest flux is less than `MAX_CORE_FLUX` (when `MAX_CORE_FLUX` is set).
   - `offset_too_large`: the absolute offset is `MAX_OFFSET_TO_ACCEPT_NM` or more.
6. Clip outliers. For each accepted line, calculate the median offset of the other accepted lines within ±`CLIP_TREND_WINDOW_NM`/2 in wavelength. The line itself is not in this median. Use the global median if the line has fewer than `MIN_TREND_NEIGHBOURS` neighbours, or if the spectrum has fewer than 10 accepted lines. Sigma-clip the differences with `CLIP_SIGMA`. The script marks the clipped lines in the `clipped` column.
7. Find one global wavelength shift for each SOXS spectrum:
   1. The coarse shift is minus the sigma-clipped median of the accepted offsets.
   2. Apply the coarse shift. Measure the lines again with the narrower window `REFINE_HALF_WINDOW_NM`, and clip outliers.
   3. The global shift is the coarse shift minus the sigma-clipped median of the remaining offsets.

   If a spectrum has no accepted lines, the script records it as failed and goes on.
8. Measure the SOXS spectrum again on the shifted wavelengths, with the same window as the refine pass (`REFINE_HALF_WINDOW_NM`). The original offsets use `SEARCH_HALF_WINDOW_NM`.
9. For each line, calculate the SOXS − X-Shooter differential: the SOXS offset minus the median X-Shooter offset of the same line.
10. Write the CSV tables, the cutouts PDF, and the plots (see [Outputs](#outputs)).

Every spread in the output is a robust σ: 1.4826 times the median absolute deviation (MAD). The script prints the shift for each spectrum and the median and standard deviation of all shifts. It also prints the number of validated lines and the X-Shooter control statistics.

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

The script needs `numpy`, `pandas`, `matplotlib`, and `astropy`. With `SHOW_FIG = True` it opens figure windows. With `SAVE_FIG = True` it saves the PNG files. Set `SHOW_FIG = False` to run without windows.

## Settings

All settings are at the top of the script.

| Constant | Current value | Meaning |
| --- | --- | --- |
| `SOXS_ROOT` | `data/soxs_extractions_0.19.0` | Folder of SOXS spectra. The script holds other choices as overwritten or commented-out lines (`data/soxs_extractions_0.17.1`, `data/tmp`). |
| `TELL_FILE` | `reference/telluric_5_lines.dat` | HITRAN line list (vacuum nm). |
| `XSHOOTER_ROOT` | `data/xshooter_extractions` | Folder of X-Shooter products. |
| `SOXS_MATCH_STRING` | `SOXS_EG21` | Text that the `OBS NAME` header must contain. `SOXS_LTT3218` and `SOXS_Feige110` are overwritten lines. |
| `SOXS_FLUX_COL` | `FLUX_COUNTS` | Flux column. The script uses `extractedFluxOptimal` if this column is missing. |
| `XSHOOTER_MATCH_STRING` | `Hip012377_TELL_SLT` | Text that the X-Shooter header must contain. |
| `XSHOOTER_FILENAME` | `SDP_TELL_SLIT_FLUX_MERGE1D_NIR.fits` | File name to search for under `XSHOOTER_ROOT`. |
| `MAX_XSHOOTER_FILES` | `20` | Maximum number of X-Shooter spectra to load. |
| `NORMALIZE`, `NORM_MODE`, `NORM_WAVE_RANGE` | `True`, `"median"`, `(1671, 1676)` | Flux normalization (wavelengths in nm). |
| `TOPPANEL_WAVE_RANGE` | `(1700, 1850)` | Wavelength range of the top panel of the overview figure. |
| `OVERVIEW_ZOOM_WINDOW` | `(1750.0, 1752.0)` | Wavelength range (nm) of the zoom panel of the overview figure. |
| `SEARCH_HALF_WINDOW_NM` | `0.25` | Half-width (nm) of the search window around each line. |
| `REFINE_HALF_WINDOW_NM` | `0.15` | Half-width (nm) of the narrower window. Used in the refine pass of the global shift and to measure the shifted offsets. |
| `MAX_OFFSET_TO_ACCEPT_NM` | `SEARCH_HALF_WINDOW_NM` | An offset must be smaller than this (nm). The script also uses it in the refine pass. |
| `MIN_POINTS_IN_WINDOW` | `5` | Fewest points in a search window. |
| `PARABOLA_REFINE` | `True` | Refine each minimum with a parabola. |
| `MIN_DEPTH` | `0.05` | Smallest line depth, in normalized flux. `None` turns the test off. |
| `MAX_CORE_FLUX` | `0.1` | A line with a lowest flux below this value is `saturated`. `None` turns the test off. |
| `CONTINUUM_PERCENTILE` | `90` | Percentile of the window flux that sets the continuum for the depth. |
| `MIN_XSH_FRACTION` | `0.5` | Smallest fraction of X-Shooter spectra that must accept a line. |
| `MIN_XSH_SPECTRA` | `3` | Fewest X-Shooter spectra for line validation. |
| `XSH_MAX_OFFSET_NM` | `0.05` | Largest absolute X-Shooter offset (nm) for a line to count as accepted. |
| `CLIP_TREND_WINDOW_NM` | `20` | Full width (nm) of the wavelength window for the local median in outlier clipping. |
| `CLIP_SIGMA` | `3.0` | Sigma limit for clipping. The script also uses it for the clipped median of the global shift. |
| `MIN_LINES_FOR_TREND` | `10` | Fewest accepted lines for a local median. With fewer, clipping uses the global median. |
| `MIN_TREND_NEIGHBOURS` | `3` | Fewest neighbours for a local median. With fewer, the line uses the global median. |
| `MAD_TO_SIGMA` | `1.4826` | Factor that changes the MAD to a σ. |
| `DIAGNOSTIC_CUTOUTS` | `True` | Write the cutouts PDF. |
| `N_ROWS`, `N_COLS` | `3`, `4` | Panels in each page of the grid plot and the cutouts PDF. |
| `XLIM_NM`, `YLIM_OFFSET_NM` | `None`, `(-0.3, 0.3)` | Axis limits. |
| `SAVE_FIG`, `SHOW_FIG` | `True`, `True` | Save the PNG files. Show the figures. |
| `OUTDIR` | `SOXS_ROOT / "plots_soxs_telluric_offsets"` | Output folder. |
| `OUTFILE_PREFIX` | `f"{SOXS_MATCH_STRING}_telluric_offsets"` (now `SOXS_EG21_telluric_offsets`) | Start of each output file name. It follows `SOXS_MATCH_STRING`. |

The code removes three older constants: `SHIFT_REFERENCE_WINDOW`, `SHIFT_POLY_DEGREE`, and `OFFSET_CLIP_NM`.

## Wavelength conversion

The HITRAN line list is in vacuum wavelengths. The SOXS and X-Shooter spectra are in air wavelengths. The function `vac_to_air_um` converts a vacuum wavelength in micrometers to an air wavelength. It divides the wavelength by the refractive index `n`. The index uses the constants `a = 8.34213e-5`, `b1 = 2.406030e-2`, `c1 = 130.0`, `b2 = 1.5997e-4`, and `c2 = 38.9`, as written in the script.

`load_hitran_air_nm` reads the file, converts nm to micrometers, converts to air, converts back to nm, and removes duplicates. The script does this once, before it uses the lines.

## Outputs

The script writes to `<SOXS_ROOT>/plots_soxs_telluric_offsets/`. It writes nothing if no SOXS spectrum gives a result.

The script always writes these CSV files, also when `SAVE_FIG = False`:

- `<OUTFILE_PREFIX>_lines.csv`. One row for each measured line of each spectrum. The SOXS spectra have an `original` and a `shifted` row for each line. The X-Shooter spectra have an `original` row. Columns: `spectrum`, `instrument`, `version`, `line_nm`, `wmin_nm`, `offset_nm`, `depth`, `core_flux`, `parabola_fallback`, `reason`, `clipped`, `differential_nm`.
- `<OUTFILE_PREFIX>_summary.csv`. One row for each spectrum, with the shift, the number of lines, the median and robust σ of the offsets and of the differentials, and a count for each rejection reason. A file that fails to load or measure gets a row with `status` set to `failed` and the error text.
- `<OUTFILE_PREFIX>_reference_lines.csv`. One row for each HITRAN line, with the X-Shooter validation result (`n_measured`, `n_accepted`, `fraction`, `median_offset_nm`, `is_valid`, `validation_available`).

When `DIAGNOSTIC_CUTOUTS = True`, the script writes `<OUTFILE_PREFIX>_cutouts.pdf`. It shows the flux around each rejected or clipped line of the first SOXS spectrum (original wavelengths). It leaves out `too_few_points` lines. It writes no file if no line is flagged.

When `SAVE_FIG = True`, the script writes these PNG files:

- `<OUTFILE_PREFIX>_page01.png`, `_page02.png`, and so on. Each page is a grid of panels with the offset to HITRAN against HITRAN wavelength for each spectrum. Each panel shows the original and the shifted offsets of the accepted, unclipped lines.
- `<OUTFILE_PREFIX>_summary.png`. Two panels. The top panel shows all offsets. The bottom panel shows the shifted SOXS − X-Shooter differentials.

A new run overwrites files that have the same name. Files from earlier runs that used a different prefix stay in the folder.

## Data sets

Each folder in `data/` holds NIR extracted, merged spectra (`*_NIR_*_EXTRACTED_MERGED.fits`) and the folder `plots_soxs_telluric_offsets/`.

| Folder | Content |
| --- | --- |
| `soxs_extractions_0.17.1/` | Spectra reduced with the tagged release 0.17.1 of `soxspipe`. |
| `soxs_extractions_0.19.0/` | The same exposures, reduced with `soxspipe` 0.19.0. It also holds `_F1` to `_F5` variants of some exposures. This set is the default `SOXS_ROOT`. |
| `xshooter_extractions/` | Archival ESO X-Shooter products. One `SDP_TELL_SLIT_FLUX_MERGE1D_NIR.fits` file in each folder named with a UUID. Files in a `_CorrOn` folder are skipped (`is_in_corron_path`). |
| `tmp/` | Working copies of single spectra, for example the output of `rerun.sh`. `tmp/hide/` holds intermediate `soxspipe` debug products for one exposure. |

The folders `soxs_extractions_0.17.1/` and `xshooter_extractions/` each hold a `list_spectra_parameters.sh` file. It lists the `DATE-OBS` header of each file with `fitsheader`. The VIS spectra of the three `soxs_extractions_*` sets are in `../soxs_uvvis_analysis/data/`.

## Compare pipeline versions

The two `soxs_extractions_*` folders hold the same standard-star exposures. Each folder has a different version of the pipeline. To compare versions:

1. Set `SOXS_ROOT` to the first folder. Run the script.
2. Set `SOXS_ROOT` to the next folder. Run the script.
3. Compare the median shift in the terminal and the plots in each `plots_soxs_telluric_offsets/` folder.

`soxs_extractions_0.17.1/` holds 12 NIR spectra and `soxs_extractions_0.19.0/` holds 24. The script keeps only the spectra that match `SOXS_MATCH_STRING`, so the count that it reports can differ between folders.

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

## Run the tests

The tests use synthetic spectra. Run them from inside this folder:

```bash
conda activate soxspipe
MPLBACKEND=Agg python -m pytest tests/
```

## Known limitation

The depth test and the saturation test can fail on some lines. The continuum is the `CONTINUUM_PERCENTILE` percentile of the flux in the search window. The saturation test uses an absolute cut (`MAX_CORE_FLUX`) on the normalized flux. Both are sensitive to a sloping continuum and to blended lines. Look at the cutouts PDF to check the rejected lines.

## Legacy

`legacy/` holds an older script and its plots. They are superseded. Do not use them.
