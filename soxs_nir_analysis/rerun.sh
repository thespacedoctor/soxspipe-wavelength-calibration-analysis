cd ../shared/project_data/
soxspipe stare sessions/base/sof/20260111T063325_NIR_3_STARE_STD_FLUX_SLIT5_0_150_0S_SOXS.sof -x
cp sessions/base/reduced/2026-01-10/soxs-stare-std/20260111T063325_NIR_3_STARE_STD_FLUX_SLIT5_0_150_0S_SOXS_EXTRACTED_MERGED.fits ../../soxs_nir_analysis/data/tmp/
cd ../../soxs_nir_analysis/
python soxs_telluric_wavelength_offsets.py
