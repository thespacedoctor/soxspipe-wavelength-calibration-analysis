# Skyline wavelength calibration anchoring

This workspace checks the wavelength calibration of the SOXS spectrograph pipeline (`soxspipe`). It compares extracted spectra of standard stars with reference lines (HITRAN telluric lines, stellar absorption lines, sky emission lines) and archival X-Shooter spectra. It also tests whether sky line fitting can anchor the wavelength solution.

## Folders

- [`soxs_uvvis_analysis/`](soxs_uvvis_analysis/README.md): notebook, reference tables, spectra, and results for the UV-VIS arm.
- [`soxs_nir_analysis/`](soxs_nir_analysis/README.md): script, reference line list, and spectra for the NIR arm.
- [`shared/`](shared/README.md): the `soxspipe` working directory for both arms, helper scripts, and the HITRAN line-list tool.
- `input/`: source material for the project. Do not write into this folder.
- `output/`: deliverables and project-adjacent files.
- `docs/`: agent and project documents.

Python packages and the `soxspipe` command are not part of this repository. Activate the conda environment `soxspipe` before you run anything.
