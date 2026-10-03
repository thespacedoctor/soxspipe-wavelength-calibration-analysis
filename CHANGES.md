# Changes

## Unreleased

- **FEATURE**: velocity-corrected per-order wavelength corrections and scatter for the UV-VIS notebook (barycentric correction, fitted star velocities, telluric anchors)
- **ENHANCEMENT**: add `PLOT_METRICS` switch for the UV-VIS notebook metrics plots, replacing an `if False:` block
- **REFACTOR**: split workspace into UV-VIS, NIR and shared folders, each with a README
- **FIXED**: correct the telluric-window widening comment in the UV-VIS notebook (±7 Å, not ±5 Å)
- **FIXED**: H-alpha no longer treated as a telluric line; duplicate windows removed; outliers clipped in UV-VIS order statistics
