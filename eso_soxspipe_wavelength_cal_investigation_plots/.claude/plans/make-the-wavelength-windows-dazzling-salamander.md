# Per-star wavelength windows + shift-statistics tables for STDCheck-Cleaned.ipynb

## Context

`STDCheck-Cleaned.ipynb` currently applies one global `WAVELENGTH_WINDOWS` list to every
star in `STANDARD_STARS`, and only ever shows median/σ of the measured wavelength shift as
inline text on each comparison plot (`plot_comparison`'s `shiftText`). In practice, not every
window is meaningful for every star — e.g. the existing `std_check_wavelength_shifts.csv`
output shows the `4680-4690 Å` window only ever produced results for `EG21`, `EG274`,
`Feige110`, `LTT7987` (not `GD71`/`GD153`/`LTT3218`), because that absorption feature isn't
present/measurable in all of them. There's also currently no consolidated statistics table
(only min/max/median/mean per star+window individually via eyeballing plots), and no way to
see whether wavelength shifts are systematic across a whole SOXS spectral order or across the
VIS/NIR detector as a whole. This change makes windows per-star, adds full summary statistics
(median/mean/std/min/max) at three levels of aggregation (star+window, spectral order, whole
detector), renders them as Markdown tables in the notebook, and writes everything (CSVs +
plot PNGs) into a new `results/` folder instead of scattering outputs at the repo root.

Decisions locked in with the user:
- `STANDARD_STARS` becomes a dict `{starName: [[min, max], ...]}`; every star defaults to the
  same window(s) currently in `WAVELENGTH_WINDOWS` (the user will hand-tune per-star windows
  later — this plan only restructures the config, it does not choose new science windows).
- Comparison + metrics plots, currently only `plt.show()`n, will now also be saved as PNGs
  under `results/plots/`.
- A window is assigned to a spectral order by checking whether its **midpoint** wavelength
  falls inside that order's range — VIS uses the non-overlapping `WLMIN`/`WLMAX` columns,
  NIR uses `WLMINFUL`/`WLMAXFUL` (its only columns, per `SOXS_SPECTRAL_FORMAT_TAB_NIR.csv`).

## File changed

Only `STDCheck-Cleaned.ipynb` (via the `NotebookEdit` tool). No other file needs to change.

## Implementation

### 1. `# INPUTS` cell

- Replace `STANDARD_STARS = [...]` + `WAVELENGTH_WINDOWS = [[4680, 4690]]` with a single dict:
  ```python
  STANDARD_STARS = {
      "GD71": [[4680, 4690]],
      "LTT3218": [[4680, 4690]],
      "GD153": [[4680, 4690]],
      "EG274": [[4680, 4690]],
      "LTT7987": [[4680, 4690]],
      "Feige110": [[4680, 4690]],
      "EG21": [[4680, 4690]],
  }
  ```
- Add new path constants, all rooted under a results folder:
  ```python
  RESULTS_DIR = "./results/"
  PLOTS_DIR = RESULTS_DIR + "plots/"
  RESULTS_CSV_PATH = RESULTS_DIR + "std_check_wavelength_shifts.csv"          # unchanged per-spectrum rows
  STAR_WINDOW_STATS_CSV_PATH = RESULTS_DIR + "std_check_stats_by_star_window.csv"
  ORDER_STATS_CSV_PATH = RESULTS_DIR + "std_check_stats_by_order.csv"
  DETECTOR_STATS_CSV_PATH = RESULTS_DIR + "std_check_stats_by_detector.csv"
  SPECTRAL_FORMAT_CSV_PATH = f"./SOXS_SPECTRAL_FORMAT_TAB_{ARM}.csv"
  ```
- Add `os.makedirs(RESULTS_DIR, exist_ok=True)` and `os.makedirs(PLOTS_DIR, exist_ok=True)` in
  this cell (or immediately after it).
- Update the "HOW TO USE" markdown cell at the top to describe the new dict shape for
  `STANDARD_STARS` and mention that outputs now live under `results/`.

### 2. `plot_comparison` / `plot_metrics`

- Add a `savePath` parameter to each; call `fig.savefig(savePath, dpi=150, bbox_inches="tight")`
  right before `plt.show()`.
- Keep everything else (including the existing median/σ `shiftText` annotation) as-is.

### 3. Driver loop

- Change the loop to `for starName, windowsForStar in STANDARD_STARS.items():` and the inner
  loop to `for windowMin, windowMax in windowsForStar:`.
- Build plot filenames per (star, window), e.g.
  `f"{PLOTS_DIR}{starName}_{windowMin:.0f}-{windowMax:.0f}A_comparison.png"` and
  `..._metrics.png`, and pass them into `plot_comparison`/`plot_metrics`.

### 4. New helper: `compute_shift_stats(resultsDf, groupColumns)`

New function (snake_case, docstring per project convention) that:
- Filters to finite `ShiftAngstrom` values (same finiteness check already used in
  `plot_metrics`).
- Groups by `groupColumns` (empty list = one overall group), aggregating `ShiftAngstrom` into
  `Count`, `MedianShiftAngstrom`, `MeanShiftAngstrom`, `StdShiftAngstrom`, `MinShiftAngstrom`,
  `MaxShiftAngstrom`.
- Returns the resulting DataFrame (reset index so group columns are regular columns).

### 5. New helper: `load_spectral_format(arm)` and `assign_spectral_order(windowMin, windowMax, spectralFormatDf)`

- `load_spectral_format(arm)`: reads `SOXS_SPECTRAL_FORMAT_TAB_{arm}.csv`, returns a DataFrame
  with columns `Order`, `RangeMinAngstrom`, `RangeMaxAngstrom` — using `WLMIN`/`WLMAX` × 10 for
  VIS, `WLMINFUL`/`WLMAXFUL` × 10 for NIR (nm → Å).
- `assign_spectral_order(windowMin, windowMax, spectralFormatDf)`: computes the window midpoint
  and returns the `Order` whose range contains it, or `np.nan` if no order matches.

### 6. New "stats + tables" cell(s) after the driver loop

Using the `resultsDf` already returned by `write_shift_results_csv`:

1. **Per star+window** — `compute_shift_stats(resultsDf, ["Star", "WindowMinAngstrom", "WindowMaxAngstrom"])`.
   Save to `STAR_WINDOW_STATS_CSV_PATH`; render as one Markdown table via
   `display(Markdown(df.to_markdown(index=False)))` (needs `from IPython.display import display, Markdown`;
   `tabulate` is already installed in the `soxspipe` env, confirmed).
2. **Per spectral order** — load `spectralFormatDf = load_spectral_format(ARM)`, add an `Order`
   column to `resultsDf` via `assign_spectral_order`, then
   `compute_shift_stats(resultsDf, ["Order"])`. Save to `ORDER_STATS_CSV_PATH`; render as a
   Markdown table. Print a note listing any (star, window) pairs whose `Order` came back NaN
   (window didn't fall in any tabulated order).
3. **Whole detector (current ARM)** — `compute_shift_stats(resultsDf, [])`, with an `Arm` column
   set to `ARM` added to the single output row. Save to `DETECTOR_STATS_CSV_PATH`; render as a
   Markdown table.

### 7. `write_shift_results_csv` call

No signature change needed — just pass the new `RESULTS_CSV_PATH` (now under `results/`).

## Verification

- Run the notebook end-to-end (`jupyter nbconvert --to notebook --execute STDCheck-Cleaned.ipynb`
  or interactively) against the currently configured `INDIR`/`ARM`.
- Confirm `results/` is created containing: `std_check_wavelength_shifts.csv`,
  `std_check_stats_by_star_window.csv`, `std_check_stats_by_order.csv`,
  `std_check_stats_by_detector.csv`, and `results/plots/*.png` (comparison + metrics PNGs per
  star/window that had matching spectra).
- Confirm the three Markdown tables render in the notebook output (not raw DataFrame reprs) and
  that their numbers are sane against the existing `median shift = ...` text already annotated
  on each comparison plot.
- Note: the stale root-level `std_check_wavelength_shifts.csv` (untracked, from a prior run) is
  left as-is — flag it to the user afterward in case they want it removed now that output moves
  under `results/`.
