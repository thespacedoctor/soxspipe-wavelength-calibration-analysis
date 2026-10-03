# Shared tools and data

This folder holds the `soxspipe` working directory and the helper scripts that both the UV-VIS and the NIR analyses use. It is for people who re-reduce data or build line lists.

## Contents

- [project_data](#project_data)
- [export.sh](#exportsh)
- [cp_data.sh](#cp_datash)
- [telluric_lines_from_hitran.py](#telluric_lines_from_hitranpy)

## project_data

`project_data/` is a `soxspipe` working directory of about 18 GB. It holds both arms. The VIS data cover 2025-08-31 to 2025-12-10. The NIR data cover 2026-01-09 to 2026-01-15. It has the folders `raw/`, `sessions/`, `misc/`, `pre_fixed_reduced/`, and `tmp/`, and the file `soxspipe.db`. The entries `sof`, `reduced`, `qc`, `soxspipe.yaml`, and `soxspipe.log` are relative symlinks into `sessions/base/`.

Run every `soxspipe` command from inside this folder. Do not move the folder or its contents out of the folder, because the symlinks are relative.

```bash
conda activate soxspipe
cd shared/project_data
soxspipe stare sessions/base/sof/<SOF_NAME>.sof -x
```

Replace `<SOF_NAME>` with the name of a `sof` file. The NIR script `../soxs_nir_analysis/rerun.sh` uses this pattern. The folder is not in version control.

## export.sh

`export.sh` makes one `sof` file for each raw exposure in a fixed list. Each line has the form `soxspipe raw sof <SOF_NAME>.sof`. The list holds 40 VIS and 12 NIR standard-star exposures. Run the script from inside `project_data/`:

```bash
cd shared/project_data
bash ../export.sh
```

The first line of the list appears twice in the file. This is harmless.

## cp_data.sh

`cp_data.sh` copies extracted spectra from the raw archive into one folder. It finds 12 NIR and 39 VIS `*_EXTRACTED_MERGED.fits` files by exact file name with `fd`, and copies each with `cp`.

Before you run it:

- The disk `/Volumes/SOXS` must be mounted. The script reads from `/Volumes/SOXS/raw_archive/soxs_commissioning_data`.
- Install `fd`.
- Change `DEST_DIR`. It still points to an old path (`/Users/Dave/Dropbox/Mac/Desktop/eso_soxspipe_wavelength_cal_investigation_plots/soxs_extractions_v.dev`) that does not exist in this workspace. Set it to a folder in `soxs_uvvis_analysis/data/` or `soxs_nir_analysis/data/`.

## telluric_lines_from_hitran.py

`telluric_lines_from_hitran.py` builds a list of telluric absorption lines (H2O and O2) from the HITRAN database. Wavelengths are in angstroms.

It needs `numpy` and `hitran-api` (`pip install hitran-api numpy`), and network access to hitran.org. HITRAN limits how many requests you can make, so the script keeps the tables it downloads in `hitran_cache/` and uses them on the next run. The cache folder is relative to the folder you run the script from, so run it from inside `shared/`.

```bash
cd shared
python telluric_lines_from_hitran.py --wmin 8000 --wmax 8400 --min-intensity 1e-25 --out telluric_8000_8400.py
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--wmin` | `8000.0` | Smallest wavelength (angstroms). |
| `--wmax` | `8400.0` | Largest wavelength (angstroms). |
| `--min-intensity` | `1e-25` | Smallest HITRAN line intensity to keep. A lower value gives more lines. |
| `--vacuum` | off | Write vacuum wavelengths. By default the script writes air wavelengths. |
| `--out` | `telluric_8000_8400.py` | Output file. The script also writes a `.txt` file with the same name that lists wavelength, intensity, and species. |

The output `.py` file defines a list called `telluric_lines`. The folder `hitran_cache/` is not in version control.

The NIR analysis does not use this script. It reads `../soxs_nir_analysis/reference/telluric_5_lines.dat`. How that file was made is not recorded in this repository.
