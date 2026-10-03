# Fix the UV-VIS wavelength-solution analysis in STDCheck-Cleaned.ipynb

## Context

`soxs_uvvis_analysis/STDCheck-Cleaned.ipynb` measures the quality of the soxspipe VIS wavelength solution. Dave needs, for each of the u, g, r and i orders:

- a median and a mean shift that can be used to correct that order;
- the standard deviation of the scatter in the wavelength solution.

The scatter is reported two ways: line-to-line scatter inside an order after removing each exposure's own offset (internal accuracy), and exposure-to-exposure scatter of the per-order offset (stability).

The method compares each SOXS spectrum with the same star's library spectrum (`reference/flux_standards.fits`). Dave confirmed that the library spectra are corrected for radial velocity and are in the rest frame. A read-only review of the code and of the existing `results/std_check_wavelength_shifts.csv` found the following.

### Findings: analysis that is wrong

1. **The velocity frames don't match (the main problem).**
   - The SOXS spectra are topocentric. The headers have no barycentric keywords, and `ESO TEL TARG RADVEL = 0`. The library is in the rest frame.
   - So every stellar-line "shift" includes three velocities: the star's radial velocity, the white-dwarf gravitational redshift, and the Earth's barycentric motion.
   - The data show this:
     - Stellar lines are offset by +0.5 to +1.3 Å (about +50 km/s).
     - Telluric lines in the same exposures are offset by 0 ± 0.2 Å.
     - Hα and the telluric line at 6278 Å are both in the r order, but Hα is offset by +0.75 to +1.1 Å while the telluric line is not. A wavelength-solution error would move both lines in the same order equally.
   - The barycentric velocities are +6 km/s (EG21), −11 to −16 km/s (EG274), +28 km/s (Feige110) and +7 to +11 km/s (LTT7987). These alone don't explain the offset. The rest is the stars' own velocity plus gravitational redshift.
   - Today's u and g "corrections" are therefore mostly stellar velocity, not calibration error.
   - (The library is in air wavelengths: its Hα and Hβ minima are at 6562.8 Å and 4861.3 Å. So there is no air-to-vacuum mismatch.)
2. **Hα is listed as a telluric line.** `[6563, 6563]` (EG274, LTT7987, EG21) and the 6562.806 Å entry in `TELLURIC_LINES` measure Hα, a stellar line, against the telluric model.
3. **Windows are duplicated.**
   - EG21 lists `[6558, 6568]` twice, so it is counted twice (6 rows from 3 spectra).
   - EG274, LTT7987 and EG21 each measure Hα twice (once as a stellar window, once as a "telluric" window).
4. **There is no outlier rejection.** A single bad fit dominates the standard deviation. For example, LTT7987 at 7253 Å has σ = 1.48 Å while the other telluric lines have σ ≈ 0.15 Å.
5. **The statistics weight things unevenly and are inconsistent.**
   - The pooled statistics weight stars by the number of spectra and windows they have.
   - The whole-detector table (cell 18) is computed before the `ORDERS_TO_INSPECT` filter.
   - The plots use `np.std` (ddof=0) and the tables use pandas `std` (ddof=1).
6. **Stellar-line minima are noisy.** For a single star, the offset varies from exposure to exposure by up to about ±60 km/s (EG274 ranges from −18 to +112 km/s). The GAM minimum of a 10 Å window on broad white-dwarf Balmer cores is poorly constrained. This needs a robust statistic and an uncertainty per line.

### Findings: things to tidy up

- `extract_properties` (cell 6) opens each FITS file twice and never closes it.
- Unused items: the `s` and `f` imports from pygam, `slit5_feige110_props`, `allTelluricResultRows`, and the `if True:` in cell 23.
- The `windowMin == windowMax` branch in `measure_window_minimum` is dead code.
- The magic number 7 is repeated in two cells.
- `print` noise in cell 7.
- The typo `STARDARD_STAR_AKAS`.
- Telluric lines are attached to individual stars, but they don't depend on the star. They should be measured in every exposure.

## Approach

Keep Dave's same-star library comparison, then remove the velocity part before computing the order statistics.

1. **Classify the lines.** Split the configuration into `STELLAR_WINDOWS` (per star, `[min, max]`) and `TELLURIC_LINES` (one shared list: 6277.70, 6868.47, 7252.95, 8228.27 Å; Hα removed). Remove the duplicates. Measure every telluric line in every exposure, against the telluric model. Add `LineType` (`stellar` or `telluric`) and `Order` to each result row.
2. **Barycentric correction.** For each exposure, compute `BaryCorrKms` with `SkyCoord(RA, DEC).radial_velocity_correction('barycentric', obstime=MJD-OBS, location=La Silla/NTT)`. Write it to the results CSV.
3. **Per-exposure order offsets from the telluric lines (r and i only).** For each exposure and order, the offset is the median telluric shift in that order.
4. **Star velocity.**
   - For each stellar line in r or i, subtract that exposure's telluric order offset.
   - Convert the result to a topocentric velocity, `v = c · Δλ / λ`, and then to the star's barycentric velocity, `v_star = v + BaryCorrKms`.
   - Take the robust median for each star over all of its exposures. The star's velocity is constant, so this pools the data.
   - Report each star's `v_star` with its MAD error, as a diagnostic. If a star has no usable r or i stellar line, flag it and leave it out of the u/g corrections.
5. **Residuals.**
   - For every stellar line: `ResidualAngstrom = ShiftAngstrom − λ · (v_star − BaryCorrKms) / c`.
   - For every telluric line: `ResidualAngstrom = ShiftAngstrom`.
   - The residual is the wavelength-solution error.
6. **Per-order statistics** (u, g, r, i):
   - The correction comes from the residuals after 3σ MAD-based clipping: `MedianShiftAngstrom`, `MeanShiftAngstrom`, `SEMeanAngstrom`, `NLines` and `NExposures`.
   - `InternalScatterStd` and `InternalScatterMadStd`: the scatter of the residuals after subtracting each exposure's own median for that order.
   - `ExposureScatterStd`: the standard deviation of the per-exposure order medians.
   - Use the same ddof everywhere: ddof=1.
   - r and i are anchored by the telluric lines alone. Their stellar residuals are centred on zero by construction, so they are reported but not used. u and g use the velocity-corrected stellar residuals. The table says which source each order uses.
   - Also keep the raw table (the shifts before velocity correction) for comparison.
7. **Diagnostic plot.** For each exposure, plot velocity offset (`c·Δλ/λ`) against wavelength, coloured by order, with the star velocity marked. A pure velocity offset is a flat line; an order offset is a step.
8. **Tidy-ups:** everything in the tidy-up list above. The plotting cells stay the same apart from reading the new columns.

### Code structure (testable)

- Move the pure calculations into a new module, `soxs_uvvis_analysis/wavecal_stats.py`:
  - `barycentric_correction_kms(header)`
  - `telluric_order_offsets(resultsDf)`
  - `fit_star_velocities(resultsDf, orderOffsets)`
  - `velocity_corrected_residuals(resultsDf, starVelocities)`
  - `clipped_order_stats(resultsDf)`
- Use snake_case for functions and camelCase for variables, following the notebook's style.
- The notebook imports from this module, so the analysis is unit-tested.
- Add `soxs_uvvis_analysis/tests/test_wavecal_stats.py` (pytest), written test-first with the **tdd-coder** agent. Use synthetic data:
  - A known velocity plus a known order offset must be recovered.
  - The clipping must reject an injected outlier.
  - The internal and exposure scatter must match the values injected.
- Reuse the notebook's existing functions:
  - `fit_abs_line` and `measure_window_minimum` (cells 8 and 9) for the line minima.
  - `load_spectral_format` and `assign_spectral_order` (cell 15) for the orders.
  - `get_spectra_for_star` (cell 9).

## Files

- `soxs_uvvis_analysis/STDCheck-Cleaned.ipynb`: the configuration split, the telluric measurement in every exposure, the module imports, the new statistics cells, the diagnostic plot, and the tidy-ups.
- `soxs_uvvis_analysis/wavecal_stats.py`: new.
- `soxs_uvvis_analysis/tests/test_wavecal_stats.py`: new.
- `soxs_uvvis_analysis/results/`:
  - `std_check_wavelength_shifts.csv` gets new columns.
  - New files: `std_check_order_corrections.csv` (the main deliverable), `std_check_order_corrections_raw.csv` and `std_check_star_velocities.csv`.
- `soxs_uvvis_analysis/README.md`: the method, the outputs, the velocity correction, and a short explanation of the barycentric correction. Written by the **doc-writer** agent.
- `CHANGES.md`.

## Process

- Work on the branch `feature/uvvis-order-corrections`.
- Before changing the code, run the **code-debugger** agent over the notebook to confirm the findings above, starting from the existing outputs.
- Write the module with **tdd-coder**, then make the notebook edits.
- Run **code-reviewer** and **python-reviewer** in parallel, then **test-runner**.
- Commit, then ask Dave before merging into `develop`.

## Verification

- `pytest soxs_uvvis_analysis/tests -q` passes, with at least 80% coverage of `wavecal_stats.py`.
- Run the notebook headless (in `soxs_uvvis_analysis/`): `MPLBACKEND=Agg jupyter nbconvert --to notebook --execute STDCheck-Cleaned.ipynb --output-dir <scratchpad>`. It must finish without errors.
- Sanity checks on the outputs:
  - In r and i, the median telluric residual is consistent with the raw telluric shifts (about 0 ± 0.2 Å).
  - The star velocities are positive and plausible for white dwarfs (tens of km/s). Feige110 (an sdO star) may differ.
  - After the velocity correction, the stellar residuals for the same line no longer depend on the star.
  - Each u and g correction is much smaller than the raw +0.5 to +1.3 Å, or the result is reported with its reason.
  - Each row of `std_check_order_corrections.csv` gives the median, the mean, the internal and exposure scatter, N, and the anchor used.
- Compare the old and new tables in the summary to Dave.

## Caveats to report

- The velocities of the u and g orders depend on the r and i stellar lines of each star: Hα, plus He I 7065 Å for Feige110. With only one or two lines per star, the error on `v_star` carries straight into the u and g corrections. That error is reported.
- The headers contain `ESO QC SKY SHIFT O1–O4`, soxspipe's own skyline shift for each order. O2 and O3 are exactly 0.0 in every file. This is worth checking with Dave: it may mean the skyline fit did not run for those orders. I'll report it and change nothing.
