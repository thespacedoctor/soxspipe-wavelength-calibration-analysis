# Changes

## Unreleased

- **FEATURE**: telluric line shift plot for every slit in the UV-VIS notebook (`PRIMARY_SLIT`, `ORDER_COLORS`, `SLIT_MARKERS`, `Slit` column); the main analysis is unchanged
- **FEATURE**: velocity-corrected per-order wavelength corrections and scatter for the UV-VIS notebook (barycentric correction, fitted star velocities, telluric anchors)
- **ENHANCEMENT**: one UV-VIS comparison plot per line and one telluric sanity-check plot per window, each with every star together in its own color (`STAR_PALETTE`, loaded stars only), replacing the per-star plots
- **ENHANCEMENT**: one panel for the UV-VIS star velocity diagnostic plot, with every star together, velocity-corrected shifts in Å, and the median and σ in the legend (`ORDER_MARKERS`)
- **ENHANCEMENT**: per-order medians in the legends of the UV-VIS telluric and velocity-corrected stellar shift diagnostic plots; legends moved outside the axes
- **ENHANCEMENT**: per-order σ in the UV-VIS line-shift plot legends, a two-line all-orders entry, comparison and telluric sanity-check legends outside the axes, and wider plots
- **ENHANCEMENT**: per-order table (N, median, σ) below each UV-VIS line-shift plot; the telluric shift plot now comes before the velocity-corrected stellar one
- **ENHANCEMENT**: draw the library and telluric model spectra on top of the SOXS spectra in the UV-VIS comparison and telluric sanity-check plots
- **ENHANCEMENT**: velocity-correct the SOXS spectra in the UV-VIS stellar-window comparison plots and show the median and σ of the corrected shift in the legend
- **ENHANCEMENT**: name stellar and telluric lines in UV-VIS plot labels (`STELLAR_LINES_NAMED`, `TELLURIC_LINES_NAMED`); fixes the Ca II H label on the H epsilon window
- **ENHANCEMENT**: add `PLOT_METRICS` switch for the UV-VIS notebook metrics plots, replacing an `if False:` block
- **REFACTOR**: split workspace into UV-VIS, NIR and shared folders, each with a README
- **FIXED**: correct the telluric-window widening comment in the UV-VIS notebook (±7 Å, not ±5 Å)
- **FIXED**: H-alpha no longer treated as a telluric line; duplicate windows removed; outliers clipped in UV-VIS order statistics
