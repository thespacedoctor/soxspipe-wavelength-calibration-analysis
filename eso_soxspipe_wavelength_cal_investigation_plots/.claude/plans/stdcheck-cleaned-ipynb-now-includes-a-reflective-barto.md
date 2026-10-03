# Simplify telluric sanity-check plots; fix telluric-line comparison colors; label lines by name

## Context

Two loose ends from the last round of styling work on `STDCheck-Cleaned.ipynb`:

1. `plot_telluric_window()`'s bottom histogram panel (measured-minimum distribution) adds little
   value now that the top panel already overlays the continuous telluric model curve for a by-eye
   check — it's being dropped to simplify the sanity-check plots down to a single panel.
2. In `plot_comparison()`, the telluric-line branch (`telluric` truthy) currently plots the
   telluric-model reference curve in green and the SOXS spectra in red (`TELLURIC_COLOR`) — the
   opposite of the now-established convention elsewhere (`plot_telluric_window`) where the telluric
   model is red and SOXS spectra are blue. This is being brought in line, including the branch's own
   histogram panel (which is kept — only the *sanity-check* plots lose their histogram, not
   `plot_comparison`'s).

Separately, `TELLURIC_LINES` currently holds a long HITRAN-derived list of bare wavelengths, used
only inside `plot_telluric_window` to draw unlabeled vertical dashed markers. That's being replaced
with a small named dict of well-known reference lines (Ca II K/H, Hβ, Mg I b1, Na D1, Hα) for use in
`plot_comparison` instead: when one of these lines falls inside the current window, the existing
black dashed "measured minimum" line's text label switches from showing the numeric wavelength to
showing the line's name (e.g. `"Hα"` instead of `"6562.81 Å"`) — no new lines are drawn, just a
better label on the line that's already there. Windows with no matching named line keep the current
numeric label as a fallback (e.g. Hδ, or the existing telluric-line-vs-telluric-model windows, which
have no entry in the new list).

## Changes

All in `STDCheck-Cleaned.ipynb`.

### 1. Replace `TELLURIC_LINES` (constants cell)

Comment out the current active list:
```python
# TELLURIC_LINES = [6867.19, 7177.39, 7191.63, 7240.56, 7290.31, 7593.70, 7633.26, ...]
```
(the two already-commented historical versions above it are untouched) and replace it with a named
dict:
```python
TELLURIC_LINES = {
    "Ca II K": 3933.66,
    "Ca II H": 3968.47,
    "H beta": 4861.33,
    "Mg I b1": 5183.64,
    "Na D1": 5895.924,
    "H alpha": 6562.8,
}
```

### 2. `plot_telluric_window()` — drop the bottom panel and the now-unused line markers

- Change `fig.add_subplot(2, 1, 1)` to a single-axes figure (`fig, ax1 = plt.subplots(figsize=(12, 5))`
  or equivalent), and reduce `figsize` height accordingly (single panel, e.g. `(12, 5)`).
- Delete the `ax2` block entirely (the `finiteMeasured`/histogram code), and move `ax1.set_xlabel("Wavelength [Å]", ...)` onto `ax1` (currently only set on the removed `ax2`, since `ax1` shared its x-axis).
- Delete the `for telluricLine in TELLURIC_LINES: ... ax1.axvline(...)` marker loop. It's now
  meaningless here — none of the new named lines (§1) fall inside any `TELLURIC_WINDOWS` range —
  and the continuous telluric-model curve already shown in this panel conveys the same information.
- `measuredMinima`/`resultRows` computation is unchanged (still needed for the returned
  `resultRows`, just no longer histogrammed).

### 3. `plot_comparison()` — swap telluric-branch colors, unify SOXS color, relabel by line name

- Reference curve: color becomes `TELLURIC_COLOR if telluric else 'green'` (currently hardcoded
  `'green'`). Keep `lw=2.5, alpha=0.6` and the existing label logic.
- Folder-spectra loop: drop the `if telluric: color=TELLURIC_COLOR ... else: color='tab:blue'`
  branch — always plot with `color='tab:blue'` (keep `lw=1.2, alpha=0.6`).
- `ax2.hist(...)`: drop its `if telluric: color=TELLURIC_COLOR` branch too — always pass
  `color='tab:blue'`.
- Label offset: change `labelX = libraryMinAngstrom + 0.02 * windowSpan` to
  `labelX = libraryMinAngstrom + 0.01 * windowSpan` (label sits closer to its line).
- Label text: before the `ax1.text(...)` call, look up whether the window contains one of the new
  named lines and use its name instead of the numeric wavelength:
  ```python
  lineLabel = f"{libraryMinAngstrom:.2f} Å"
  for lineName, lineWaveAngstrom in TELLURIC_LINES.items():
      if windowMin < lineWaveAngstrom < windowMax:
          lineLabel = lineName
          break
  ```
  and pass `lineLabel` to `ax1.text(...)` in place of the current inline f-string. Windows with no
  matching entry (e.g. Hδ, or the telluric-line-vs-telluric-model windows) keep the numeric label.

## Verification

- Re-run the notebook top-to-bottom in the `soxspipe` conda env (`conda run -n soxspipe jupyter
  nbconvert --to notebook --execute ...`) to confirm it executes cleanly (in particular that nothing
  downstream still iterates `TELLURIC_LINES` expecting a flat list).
- Inspect a regenerated `results/plots/*_telluric.png` (e.g. `LTT7987_6860-6890A_telluric.png`) to
  confirm it's now a single panel with no bottom histogram and no vertical line markers.
- Inspect a regenerated telluric-line `results/plots/*_comparison.png` (e.g.
  `LTT7987_telluric_line_8227A_comparison.png`) to confirm the telluric-model curve is now red and
  the SOXS spectra/histogram are blue.
- Inspect a regenerated stellar-line `results/plots/*_comparison.png` whose window contains one of
  the six named lines (e.g. `LTT7987_6558-6568A_comparison.png`, which spans Hα at 6562.8) to
  confirm the dashed line's label now reads `"H alpha"` instead of the numeric wavelength, and sits
  closer to the line than before. Also check one with no matching line (e.g. the Hδ or He II
  windows) still shows the numeric label.
