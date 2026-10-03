# Reorganise the workspace into UV-VIS and NIR analysis folders

## Context

The repo is hard to navigate. All analysis lives in one folder, `eso_soxspipe_wavelength_cal_investigation_plots/` (W/). In that folder, UV-VIS and NIR code, reference tables, extracted spectra, plots, and an 18 GB soxspipe working directory sit side by side. Two VIS soxspipe working directories also sit at the repo root. Several paths are stale. Dave wants the UV-VIS and NIR analyses in separate folders, each with a README that explains its scripts. The older NIR script (`soxs_plot_spectrum_tell_vs_xshooter_shift.py`) is left out of the write-up.

Decisions confirmed with Dave:
- Use the layout below, with `reference/`, `data/` and `results/` subfolders.
- Split the three mixed `soxs_nir_extractions_*` folders by arm, using `_VIS_` / `_NIR_` in the filenames.
- Keep `project_data/` whole in `shared/`, because its `soxspipe.db` and sessions cover both arms.
- Move the older NIR script to `soxs_nir_analysis/legacy/` and leave it out of the README.

## Target layout

```
./
├── README.md                  index of the three folders (rewrite the empty stub)
├── CLAUDE.md                  update for the new layout
├── soxs_uvvis_analysis/
│   ├── README.md
│   ├── STDCheck-Cleaned.ipynb
│   ├── reference/             flux_standards.fits, soxs_vis_r5000_telluric_spectrum.tsv,
│   │                          SOXS_SPECTRAL_FORMAT_TAB_{VIS,NIR}.csv, slit5_feige110_props_all_time{,_orig}.csv
│   ├── data/                  2026_june_july_standards_{attached_arcs,unattached_arcs,unattached_arcs_w_skyline_fitting}/
│   │                          soxs_extractions_{0.17.1,v.dev,v.dev_horne_fix}/   (VIS FITS only)
│   │                          tmp/   (the single VIS file 20260628T053000 from extraction_archive/tmp)
│   ├── results/               CSVs and plots/ (unchanged)
│   ├── reports/               01_/02_ *.html, *.pdf
│   └── soxspipe_workdirs/     2026_june_VIS_example/, VIS_attached_arc_example/,
│                              20260707T071228_soxs_absorption_lines.reg
├── soxs_nir_analysis/
│   ├── README.md
│   ├── soxs_telluric_wavelength_offsets.py
│   ├── rerun.sh
│   ├── reference/             telluric_5_lines.dat
│   ├── data/                  soxs_extractions_{0.17.1,v.dev,v.dev_horne_fix}/   (NIR FITS, plots_soxs_telluric_offsets/, list_spectra_parameters.sh)
│   │                          xshooter_extractions/, tmp/ (the NIR files and hide/)
│   └── legacy/                soxs_plot_spectrum_tell_vs_xshooter_shift.py, 2026-01/ (its plots)
├── shared/
│   ├── README.md
│   ├── project_data/          soxspipe working dir with VIS and NIR (moved whole; its symlinks are relative)
│   ├── telluric_lines_from_hitran.py, hitran_cache/
│   └── cp_data.sh, export.sh
└── docs/  input/  output/
```

`W/` is removed when it is empty. The extraction folders drop "nir" from their names, because each copy now holds a single arm.

## Steps

1. **Branch:** create `feature/reorganise-by-arm` from `develop`.
2. **Move files:** use `git mv` for tracked files and plain `mv` for ignored data. Everything is on the same volume, so each move is a rename and no data is copied. To split each extraction folder, move `*_VIS_*.fits` to the UV-VIS `data/` folder and everything else to the NIR `data/` folder. Check the file counts before and after (about 38 VIS and 12 NIR per folder).
3. **Fix paths in the code.** Make only the path edits listed here; change no logic.
   - `STDCheck-Cleaned.ipynb` cell 2, cell 4 and cell 15 (the paths found in the survey):
     - `./extraction_archive/2026_june_july_*` → `./data/2026_june_july_*` (the live `INDIR` and the two commented-out options).
     - `./soxs_vis_r5000_telluric_spectrum.tsv`, `./flux_standards.fits`, `./slit5_feige110_props_all_time.csv` and `./SOXS_SPECTRAL_FORMAT_TAB_{ARM}.csv` → `./reference/…`.
     - Edit the cell sources with a small Python/JSON script, so the outputs stay as they are.
   - `soxs_telluric_wavelength_offsets.py` L13–17:
     - Remove the dead first `SOXS_ROOT` assignment (L13).
     - Set `SOXS_ROOT = Path("data/soxs_extractions_v.dev_horne_fix")`, and make the commented-out alternative `data/tmp`.
     - Set `TELL_FILE = Path("reference/telluric_5_lines.dat")`.
     - Set `XSHOOTER_ROOT = Path("data/xshooter_extractions")`. This also fixes the stale path.
   - `rerun.sh`: `cd ../shared/project_data/`. Copy to `../../soxs_nir_analysis/data/tmp/`, `cd` back, then run the script.
   - Legacy script: leave it unchanged.
   - `shared/cp_data.sh`: leave the logic as it is. Its `DEST_DIR` is an old Desktop path; the README says so.
   - The three `soxspipe.yaml` files (line 884, log `filename:`): change the stale Desktop path to each folder's new absolute path. These files are ignored by git.
4. **Fold the old subfolder guide into the new docs:**
   - Move the content of `W/.claude/CLAUDE.md` into the READMEs. The method steps go to NIR, and the standard-star line table goes to UV-VIS.
   - Move `W/.claude/plans/*.md` to the root `.claude/plans/`.
   - Move `W/.claude/settings.local.json` to the root `.claude/` if no root copy exists. If one exists, merge the `permissions` entries.
   - Delete `W/.claude/`.
5. **Write the READMEs** with the **doc-writer** agent, following `rules/common/documentation.md`:
   - `soxs_uvvis_analysis/README.md`:
     - What the notebook measures: the wavelength shift of standard-star absorption and telluric lines against `flux_standards.fits` and the telluric model.
     - The config cell: `ARM`, `INDIR`, `STANDARD_STARS` windows.
     - Inputs and outputs (`results/*.csv`, `results/plots/*`).
     - What the three `2026_june_july_*` input sets are.
     - The standard-star line table.
     - The two soxspipe working directories and how they differ (`subtract_sky`).
     - How to run: `conda activate soxspipe`, then Jupyter.
   - `soxs_nir_analysis/README.md`:
     - What `soxs_telluric_wavelength_offsets.py` does: the five method steps, the `# INPUTS` constants, the HITRAN vacuum-to-air conversion, and where outputs go.
     - The data sets in `data/`, including the pipeline-version comparison.
     - How `rerun.sh` works.
     - A one-line note that `legacy/` holds a superseded script.
   - `shared/README.md`:
     - `project_data/`: a soxspipe working directory with both arms. Run `soxspipe` from inside it.
     - `export.sh` (makes sof files; run it inside `project_data/`).
     - `cp_data.sh`, including the external volume and the stale destination.
     - `telluric_lines_from_hitran.py`: its CLI and cache.
   - Root `README.md`: a short index of the three folders.
6. **Update `.gitignore`:**
   - Replace the `/eso_soxspipe_wavelength_cal_investigation_plots/…` rules with arm-folder rules: `*/data/`, `shared/project_data/`, `shared/hitran_cache/`, `soxs_uvvis_analysis/results/plots/`, `soxs_uvvis_analysis/soxspipe_workdirs/`, `soxs_uvvis_analysis/reference/*.tsv`, `soxs_nir_analysis/legacy/2026-01/`.
   - Keep the global `*.fits`, `*.png`, `*.pdf`, `*.html` and `*.log` rules.
7. **Update the root `CLAUDE.md`:** replace the Layout, Common commands, Pipeline-version comparison and Known stale paths sections to match the new tree, and point to the three READMEs. Remove the stale-path entries that this work fixes.
8. **Review:** run **code-reviewer** + **python-reviewer** (one `.py` change and the notebook) in parallel. Then commit with `**REFACTOR**: split workspace into UV-VIS, NIR and shared folders`, and add a `CHANGES.md` line. There's no remote, so nothing is pushed and no PR is opened. Then ask Dave about merging into `develop`.

## Verification

- `git status` / `git diff --stat -M`: tracked files show as renames, plus the expected path edits.
- `find . -xtype l` returns nothing, so no symlinks are broken.
- File counts in `data/` match the counts before the split. `du -sh` totals are unchanged at about 20 GB.
- `grep -rn "extraction_archive\|eso_soxspipe_wavelength\|xshooter_extractions\"" --include=*.py --include=*.sh --include=*.ipynb --include=*.md .` matches only the deliberate mentions (none in live code).
- In `soxs_nir_analysis/`, run `conda activate soxspipe && python soxs_telluric_wavelength_offsets.py`. It must load the spectra, the X-Shooter files and the HITRAN lines, and write PNGs to `data/soxs_extractions_v.dev_horne_fix/plots_soxs_telluric_offsets/`.
- In `soxs_uvvis_analysis/`, run `jupyter nbconvert --to notebook --execute STDCheck-Cleaned.ipynb --output <scratchpad>/check.ipynb`. It must finish without errors. `git diff --stat results/` should show only the expected rewrites of the CSVs.
- `bash -n soxs_nir_analysis/rerun.sh` passes. I won't run the full re-reduction.

## Note

The folder lives in Dropbox. Renaming about 20 GB triggers a Dropbox re-index but no re-upload. Let it finish before you open the files on another machine.
