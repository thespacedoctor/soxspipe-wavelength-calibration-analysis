# Changes

## Unreleased

- **FEATURE**: telluric line shift plot for every slit in the UV-VIS notebook (`PRIMARY_SLIT`, `ORDER_COLORS`, `SLIT_MARKERS`, `Slit` column); the main analysis is unchanged
- **FEATURE**: velocity-corrected per-order wavelength corrections and scatter for the UV-VIS notebook (barycentric correction, fitted star velocities, telluric anchors)
- **FEATURE**: NIR telluric script writes SOXS − X-Shooter differentials, lines/summary/reference_lines CSVs, and a diagnostic cutout PDF of rejected and clipped lines
- **ENHANCEMENT**: UV-VIS telluric and stellar shift diagnostic plots share one x-range and one key: order bands in the background, marker color for slit (`SLIT_COLORS`), marker shape for star (`STAR_MARKER_PALETTE`); the stellar plot now includes every slit and `ORDER_MARKERS` and `SLIT_MARKERS` are removed
- **ENHANCEMENT**: one two-panel UV-VIS shift diagnostic figure (telluric above stellar, `shift_diagnostic.png`) with dashed lines between the orders, a fixed y-axis range (`SHIFT_YLIM_ANGSTROM`), and one by-order table with telluric rows first
- **ENHANCEMENT**: one UV-VIS comparison plot per line and one telluric sanity-check plot per window, each with every star together in its own color (`STAR_PALETTE`, loaded stars only), replacing the per-star plots
- **ENHANCEMENT**: one panel for the UV-VIS star velocity diagnostic plot, with every star together, velocity-corrected shifts in Å, and the median and σ in the legend (`ORDER_MARKERS`)
- **ENHANCEMENT**: per-order medians in the legends of the UV-VIS telluric and velocity-corrected stellar shift diagnostic plots; legends moved outside the axes
- **ENHANCEMENT**: per-order σ in the UV-VIS line-shift plot legends, a two-line all-orders entry, comparison and telluric sanity-check legends outside the axes, and wider plots
- **ENHANCEMENT**: per-order table (N, median, σ) below each UV-VIS line-shift plot; the telluric shift plot now comes before the velocity-corrected stellar one
- **ENHANCEMENT**: draw the library and telluric model spectra on top of the SOXS spectra in the UV-VIS comparison and telluric sanity-check plots
- **ENHANCEMENT**: velocity-correct the SOXS spectra in the UV-VIS stellar-window comparison plots and show the median and σ of the corrected shift in the legend
- **ENHANCEMENT**: name stellar and telluric lines in UV-VIS plot labels (`STELLAR_LINES_NAMED`, `TELLURIC_LINES_NAMED`); fixes the Ca II H label on the H epsilon window
- **ENHANCEMENT**: add `PLOT_METRICS` switch for the UV-VIS notebook metrics plots, replacing an `if False:` block
- **ENHANCEMENT**: NIR telluric script rejects false minima (edge, shallow, saturated) and measures only HITRAN lines validated against X-Shooter spectra
- **ENHANCEMENT**: NIR telluric depth and saturation use a sideband continuum that scales with a per-spectrum line width estimate; `continuum_mode` marks one-sided fits
- **REFACTOR**: split workspace into UV-VIS, NIR and shared folders, each with a README
- **DOCS**: update the NIR README for the rewritten telluric script (method, settings, outputs, tests, known limitation)
- **FIXED**: correct the telluric-window widening comment in the UV-VIS notebook (±7 Å, not ±5 Å)
- **FIXED**: H-alpha no longer treated as a telluric line; duplicate windows removed; outliers clipped in UV-VIS order statistics
- **FIXED**: NIR global shift came from one line with a degree-12 polynomial (a spike at 1751.1 nm gave +0.13 to +0.50 nm); it is now the clipped median of per-line offsets
- **FIXED**: NIR telluric script crash on an empty window, overwritten rejection reasons, grid plot label mismatch, and `OUTFILE_PREFIX` not following `SOXS_MATCH_STRING`
- **TEST**: pytest suite on synthetic spectra for the NIR telluric script
