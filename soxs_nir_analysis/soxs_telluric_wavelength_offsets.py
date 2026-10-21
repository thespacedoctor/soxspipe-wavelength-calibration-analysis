#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from astropy.io import fits
import pandas as pd

# =========================
# INPUTS
# =========================
SOXS_ROOT = Path("data/soxs_extractions_v.dev_horne_fix")
# SOXS_ROOT = Path("data/tmp")
TELL_FILE = Path("reference/telluric_5_lines.dat")
XSHOOTER_ROOT = Path("data/xshooter_extractions")

SOXS_MATCH_STRING = "SOXS_LTT3218"
# SOXS_MATCH_STRING = "SOXS_Feige110"
SOXS_FLUX_COL = "FLUX_COUNTS"
XSHOOTER_MATCH_STRING = "Hip012377_TELL_SLT"
XSHOOTER_FILENAME = "SDP_TELL_SLIT_FLUX_MERGE1D_NIR.fits"
MAX_XSHOOTER_FILES = 20

NORMALIZE = True
NORM_MODE = "median"
NORM_WAVE_RANGE = (1671, 1676)

TOPPANEL_WAVE_RANGE = (1700, 1850)
SHIFT_REFERENCE_WINDOW = (1750.0, 1752.0)
SHIFT_POLY_DEGREE = 12

SEARCH_HALF_WINDOW_NM = 0.25
MAX_OFFSET_TO_ACCEPT_NM = SEARCH_HALF_WINDOW_NM
OFFSET_CLIP_NM = 0.19

MIN_POINTS_IN_WINDOW = 5
PARABOLA_REFINE = True
MIN_DEPTH = None

N_ROWS = 3
N_COLS = 4
PANELS_PER_FIG = N_ROWS * N_COLS

XLIM_NM = None
YLIM_OFFSET_NM = (-0.3, 0.3)

SAVE_FIG = True
SHOW_FIG = True

OUTDIR = SOXS_ROOT / "plots_soxs_telluric_offsets"
OUTFILE_PREFIX = "SOXS_LTT3218_telluric_offsets_grid"

COLOR_SOXS = "tab:cyan"
COLOR_SOXS_SHIFTED = "tab:blue"


def format_axes_no_scientific(ax):
    # DISABLE SCIENTIFIC AXIS LABEL FORMATTING.
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)


# =========================
# BASIC UTILITIES
# =========================
def vac_to_air_um(wave_um):
    # CONVERT VACUUM WAVELENGTHS (UM) TO AIR WAVELENGTHS (UM).
    wave_um = np.asarray(wave_um, dtype=float)

    a = 8.34213e-5
    b1 = 2.406030e-2
    c1 = 130.0
    b2 = 1.5997e-4
    c2 = 38.9

    n = 1.0 + a + wave_um**2 * b1 / (wave_um**2 * c1 - 1.0) + wave_um**2 * b2 / (wave_um**2 * c2 - 1.0)

    return wave_um / n


def load_hitran_air_nm(filepath: Path):
    # LOAD HITRAN LINES IN VACUUM NM AND CONVERT TO AIR NM.
    if not filepath.exists():
        raise FileNotFoundError(f"Missing HITRAN/telluric file: {filepath}")

    lines_nm_vac = np.loadtxt(filepath, dtype=float)
    lines_nm_vac = np.atleast_1d(lines_nm_vac)
    lines_nm_vac = lines_nm_vac[np.isfinite(lines_nm_vac)]

    lines_um_air = vac_to_air_um(lines_nm_vac / 1000.0)
    lines_nm_air = lines_um_air * 1000.0

    return np.unique(lines_nm_air)


def get_header_value(header, keys, default="UNKNOWN"):
    # RETURN FIRST NON-EMPTY HEADER VALUE FROM CANDIDATE KEYS.
    for key in keys:
        if key in header:
            val = str(header[key]).strip()
            if val:
                return val
    return default


def get_date_obs(header):
    # RETURN BEST-AVAILABLE OBSERVATION DATE KEY.
    return get_header_value(
        header,
        ["DATE-OBS", "DATE_OBS", "HIERARCH ESO OBS START", "MJD-OBS"],
        default="UNKNOWN",
    )


def clean_spectrum(wave, flux):
    # FILTER NON-FINITE VALUES, SORT BY WAVELENGTH, AND DROP DUPLICATES.
    wave = np.asarray(wave, dtype=float).ravel()
    flux = np.asarray(flux, dtype=float).ravel()

    good = np.isfinite(wave) & np.isfinite(flux)
    wave = wave[good]
    flux = flux[good]

    order = np.argsort(wave)
    wave = wave[order]
    flux = flux[order]

    wave, idx = np.unique(wave, return_index=True)
    flux = flux[idx]

    return wave, flux


def normalize_flux(wave_nm, flux, filepath):
    # NORMALIZE FLUX WITHIN A CONFIGURED WAVELENGTH WINDOW.
    if not NORMALIZE:
        return flux

    if NORM_WAVE_RANGE is not None:
        wmin, wmax = NORM_WAVE_RANGE
        mask = (wave_nm >= wmin) & (wave_nm <= wmax) & np.isfinite(flux)
    else:
        mask = np.isfinite(flux)

    if np.sum(mask) < 5:
        raise ValueError(f"Not enough points in normalization range {NORM_WAVE_RANGE} for {filepath}")

    if NORM_MODE == "median":
        norm = np.nanmedian(flux[mask])
    elif NORM_MODE == "mean":
        norm = np.nanmean(flux[mask])
    elif NORM_MODE == "max":
        norm = np.nanmax(flux[mask])
    else:
        raise ValueError(f"Unknown NORM_MODE: {NORM_MODE}")

    if not np.isfinite(norm) or norm == 0:
        raise ValueError(f"Invalid normalization value {norm} for {filepath}")

    return flux / norm


def chunks(seq, size):
    # YIELD FIXED-SIZE CHUNKS WITH 1-BASED PAGE INDEX.
    for i in range(0, len(seq), size):
        yield i // size + 1, seq[i : i + size]


def clipped_good_mask(results):
    # KEEP ONLY GOOD OFFSETS AND REJECT LOCAL OUTLIERS.
    results = pd.DataFrame(results)

    mask = results["ok"] == 1
    results = results.loc[mask]

    # COMPARE EACH OFFSET TO A LOCAL ROLLING MEDIAN TREND.
    results["offset_nm_median"] = np.nan
    results["offset_nm_median"] = results["offset_nm"].rolling(window=7, center=True, closed="both").median()

    results["offset_nm_median_diff"] = results["offset_nm_median"] - results["offset_nm"]
    results["offset_nm_median_diff"] = results["offset_nm_median_diff"].abs()
    mask = results["offset_nm_median_diff"].isnull() == False
    results = results.loc[mask]

    # APPLY SIGMA CLIPPING TO THE LOCAL-MEDIAN DEVIATIONS.
    from astropy.stats import sigma_clipped_stats

    mean, _, std = sigma_clipped_stats(
        results["offset_nm_median_diff"], sigma=2.0, stdfunc="std", cenfunc="mean", maxiters=3
    )

    mask = results["offset_nm_median_diff"] <= mean + 2 * std
    results = results.loc[mask]

    return results


# =========================
# READ SOXS
# =========================
def read_soxs_spectrum(filepath: Path):
    # READ ONE SOXS SPECTRUM AND STANDARDIZE ITS FIELDS.
    with fits.open(filepath) as hdul:
        header = hdul[0].header
        data = hdul[1].data
        colnames = data.columns.names

        # FILTER TO THE TARGET OBSERVATION NAME.
        obs_name = get_header_value(
            header,
            ["HIERARCH ESO OBS NAME", "OBS.NAME", "OBSNAME", "OBJECT"],
        )

        if SOXS_MATCH_STRING.upper() not in obs_name.upper():
            return None

        # CHOOSE THE BEST AVAILABLE WAVELENGTH COLUMN.
        if "WAVE" in colnames:
            wave_col = "WAVE"
        elif "wavelengthMean" in colnames:
            wave_col = "wavelengthMean"
        else:
            raise KeyError(f"No WAVE/wavelengthMean column in {filepath}")

        # CHOOSE THE PREFERRED FLUX COLUMN WITH FALLBACK.
        flux_col = SOXS_FLUX_COL
        if flux_col not in colnames:
            flux_col = "extractedFluxOptimal"

        if flux_col not in colnames:
            raise KeyError(f"No {flux_col} column in {filepath}. Columns: {colnames}")

        wave_nm = np.asarray(data[wave_col], dtype=float)
        flux = np.asarray(data[flux_col], dtype=float)
        date_obs = get_date_obs(header)

    wave_nm, flux = clean_spectrum(wave_nm, flux)
    flux = normalize_flux(wave_nm, flux, filepath)

    # LABEL THE OBSERVING MODE FROM THE FILENAME.
    if "STARE" in filepath.name.upper():
        mode = "stare"
    elif "NOD" in filepath.name.upper():
        mode = "nodding"
    else:
        mode = "unknown"

    return {"instrument": "SOXS", "wave": wave_nm, "flux": flux, "date": date_obs, "file": filepath, "mode": mode}


def find_soxs_spectra():
    # SCAN THE SOXS DIRECTORY AND LOAD MATCHING FILES.
    files = sorted(SOXS_ROOT.glob("*_EXTRACTED_ME*.fits"))
    spectra = []

    # SKIP KNOWN BAD FILES AND KEEP ONLY READABLE SPECTRA.
    for f in files:
        if "NOTFLAT" in f.name or "ERROR" in f.name:
            continue

        try:
            spec = read_soxs_spectrum(f)
        except Exception as e:
            print(f"[SOXS ERROR] {f}")
            print(f"             {e}")
            continue

        if spec is not None:
            spectra.append(spec)

    return spectra


# =========================
# MINIMUM MEASUREMENT
# =========================
def fit_local_parabola_minimum(x, y):
    # FIT A LOCAL PARABOLA AND RETURN ITS VERTEX AS THE REFINED MINIMUM.
    good = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x, dtype=float)[good]
    y = np.asarray(y, dtype=float)[good]

    if x.size < 3:
        return np.nan, np.nan, None, None

    try:
        a, b, c0 = np.polyfit(x, y, 2)

        if not np.isfinite(a) or a <= 0:
            return np.nan, np.nan, None, None

        x0 = -b / (2.0 * a)
        y0 = a * x0**2 + b * x0 + c0

        xfit = np.linspace(np.nanmin(x), np.nanmax(x), 50)
        yfit = a * xfit**2 + b * xfit + c0
        return x0, y0, xfit, yfit

    except Exception:
        return np.nan, np.nan, None, None


def find_local_minimum_near_line(
    wave_nm,
    flux,
    line_nm,
    half_window_nm=0.25,
    min_points=5,
    parabola_refine=True,
    min_depth=None,
    max_offset_nm=0.25,
    plot=False,
    ax=None,
    shift=0.0,
):
    # EXTRACT THE LOCAL WINDOW AROUND ONE REFERENCE LINE.
    m = (
        np.isfinite(wave_nm)
        & np.isfinite(flux)
        & (wave_nm >= line_nm - half_window_nm)
        & (wave_nm <= line_nm + half_window_nm)
    )

    w = wave_nm[m]
    f = flux[m]

    error = None

    # CREATE A LOCAL AXIS ONLY WHEN NEEDED.
    if plot and ax is None:
        _, ax = plt.subplots(figsize=(6, 3.5))
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel("Flux")
        format_axes_no_scientific(ax)

    # SHOW THE LOCAL SAMPLES AND THE TARGET LINE.
    if plot and ax is not None:
        ax.plot(w, f, "o-", ms=3, lw=1, color="tab:blue", label="local flux")
        ax.axvline(line_nm, color="tab:orange", ls="--", lw=1, label="reference line")
        format_axes_no_scientific(ax)

    # REQUIRE ENOUGH SAMPLES TO TRUST THE LOCAL FIT.
    if w.size < min_points:
        if plot and ax is not None:
            ax.set_title("too_few_points")
        error = "too_few_points"

    # START FROM THE LOWEST DISCRETE SAMPLE.
    i0 = np.nanargmin(f)
    wmin = w[i0]
    fmin = f[i0]

    # OPTIONALLY REJECT SHALLOW FEATURES.
    if min_depth is not None and error is None:
        depth = np.nanmedian(f) - fmin
        if not np.isfinite(depth) or depth < min_depth:
            if plot and ax is not None:
                ax.plot(wmin, fmin, "x", ms=8, mew=2, color="tab:red", label="minimum")
                ax.set_title("too_shallow")
                print("too_shallow")
            error = "too_shallow"

    # REFINE THE MINIMUM WITH A 3-POINT PARABOLA.
    if parabola_refine and 1 <= i0 <= w.size - 2 and error is None:
        xloc = w[i0 - 1 : i0 + 2]
        yloc = f[i0 - 1 : i0 + 2]
        xfit, yfit, xfit_coords, yfit_coords = fit_local_parabola_minimum(xloc, yloc)
        if plot:
            ax.plot(xloc, yloc, "o-", ms=5, color="tab:red", label="parabola points", zorder=2)
            ax.plot(xfit_coords, yfit_coords, "o", ms=1, color="tab:cyan", label="parabola fitted")

        if np.isfinite(xfit) and xloc.min() <= xfit <= xloc.max():
            wmin = xfit
            fmin = yfit
    else:
        if plot:
            ax.set_title("parabola_refine_failed")
        error = "parabola_refine_failed"

    # COMPUTE AND VALIDATE THE FINAL OFFSET.
    offset_nm = wmin - line_nm

    if error is None and (not np.isfinite(offset_nm) or abs(offset_nm) >= max_offset_nm):
        if plot and ax is not None:
            ax.plot(wmin, fmin, "x", ms=8, mew=2, color="tab:red", label="minimum")
            ax.set_title("offset_too_large")
            print("offset_too_large")
        error = "offset_too_large"

    # MARK THE ACCEPTED MINIMUM ON THE PLOT.
    if plot and ax is not None:
        ax.plot(wmin, fmin, "x", ms=8, mew=2, color="tab:green", label="accepted minimum", zorder=3)

    # LABEL SUCCESSFUL MEASUREMENTS WITH THE APPLIED SHIFT.
    if error is None and plot and ax is not None:
        if shift:
            shift_text = f" (shift applied: {shift:+.4f} nm)"
        else:
            shift_text = ""
        ax.set_title(f"ok: offset={offset_nm:+.4f} nm{shift_text}")

    if plot and ax is not None:
        ax.legend(loc="best", fontsize=8)
        plt.show()
        plt.close("all")

    if error is not None:
        return False, np.nan, np.nan, error

    return True, wmin, offset_nm, "ok"


def measure_reference_minimum_for_shift(spec, hitran_ref_nm, plot=False, ax=None):
    # MEASURE THE REFERENCE MINIMUM USED FOR GLOBAL SHIFTING.
    if "wave" in spec:
        wave = spec["wave"]
    elif "wavelengthMean" in spec:
        wave = spec["wavelengthMean"]

    flux = spec["flux"]

    # LIMIT THE FIT TO THE REFERENCE WINDOW.
    w0, w1 = SHIFT_REFERENCE_WINDOW
    mask = (wave >= w0) & (wave <= w1) & np.isfinite(flux)

    if np.sum(mask) < SHIFT_POLY_DEGREE + 2:
        raise ValueError(f"Not enough points in {SHIFT_REFERENCE_WINDOW} nm for {spec['file']}")

    # FIT A SMOOTH POLYNOMIAL MODEL TO THE LOCAL FEATURE.
    x = wave[mask]
    y = flux[mask]

    x_ref = np.nanmean(x)
    coeff = np.polyfit(x - x_ref, y, SHIFT_POLY_DEGREE)
    poly = np.poly1d(coeff)

    # EVALUATE THE MODEL ON A FINE GRID TO FIND THE MINIMUM.
    fine_wave = np.linspace(w0, w1, 3000)
    fine_flux = poly(fine_wave - x_ref)

    min_wave = fine_wave[np.nanargmin(fine_flux)]
    shift_nm = hitran_ref_nm - min_wave

    if plot:
        import matplotlib.pyplot as plt

        created_fig = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
            created_fig = True

        # PLOT ORIGINAL DATA POINTS.
        ax.scatter(x, y, s=10, color="k", alpha=0.6, label="data")

        # PLOT POLYNOMIAL FIT.
        ax.plot(fine_wave, fine_flux, color="C0", lw=1.5, label=f"poly deg={SHIFT_POLY_DEGREE}")

        # MARK FITTED MINIMUM.
        ax.plot(min_wave, poly(min_wave - x_ref), marker="x", ms=8, mew=2, color="tab:green", label="measured min")

        # MARK HITRAN REFERENCE LINE.
        ax.axvline(hitran_ref_nm, color="tab:red", linestyle="--", label="HITRAN ref")

        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel("Flux")
        try:
            format_axes_no_scientific(ax)
        except Exception:
            pass
        ax.legend(loc="best", fontsize=8)

        if created_fig:
            plt.show()
            plt.close("all")

    return min_wave, shift_nm


def nearest_hitran_reference_line(tell_nm_air):
    # PICK THE HITRAN LINE CLOSEST TO THE SHIFT-WINDOW CENTER.
    w0, w1 = SHIFT_REFERENCE_WINDOW
    center = 0.5 * (w0 + w1)

    lines = tell_nm_air[(tell_nm_air >= w0) & (tell_nm_air <= w1)]

    if lines.size == 0:
        raise RuntimeError(f"No HITRAN air line found in {SHIFT_REFERENCE_WINDOW}")

    return lines[np.argmin(np.abs(lines - center))]


def apply_soxs_shift(soxs, hitran_ref_nm, plot=False):
    # DERIVE AND APPLY A GLOBAL SHIFT TO EACH SOXS SPECTRUM.
    shifted = []

    # KEEP ONLY SPECTRA WITH A MEASURABLE REFERENCE FEATURE.
    for spec in soxs:
        try:
            min_wave, shift_nm = measure_reference_minimum_for_shift(spec, hitran_ref_nm, plot=plot)

            s = dict(spec)
            s["wave_shifted"] = spec["wave"] + shift_nm
            s["shift_nm"] = shift_nm
            s["shift_reference_min_nm"] = min_wave
            s["shift_reference_hitran_nm"] = hitran_ref_nm

            shifted.append(s)

        except Exception as e:
            print(f"[SOXS SHIFT ERROR] {spec['file']}")
            print(f"                  {e}")

    return shifted


# =========================
# MEASURE ALL HITRAN OFFSETS
# =========================
def measure_offsets_for_spectrum(wave_nm, flux, lines_nm, plot=False, shift=0.0):
    # MEASURE OFFSETS FOR ALL RELEVANT HITRAN LINES.
    out = {
        "line_nm": [],
        "wmin_nm": [],
        "offset_nm": [],
        "ok": [],
        "reason": [],
    }

    # RESTRICT THE LINE LIST TO THE SPECTRAL COVERAGE.
    if wave_nm.size == 0:
        lines_here = np.array([])
    else:
        lines_here = lines_nm[(lines_nm >= np.nanmin(wave_nm)) & (lines_nm <= np.nanmax(wave_nm))]

    # MEASURE EACH LINE INDEPENDENTLY.
    for line_nm in lines_here:
        ok, wmin, offset, reason = find_local_minimum_near_line(
            wave_nm,
            flux,
            line_nm,
            half_window_nm=SEARCH_HALF_WINDOW_NM,
            min_points=MIN_POINTS_IN_WINDOW,
            parabola_refine=PARABOLA_REFINE,
            min_depth=MIN_DEPTH,
            max_offset_nm=MAX_OFFSET_TO_ACCEPT_NM,
            plot=plot,
            shift=shift,
        )

        out["line_nm"].append(line_nm)
        out["wmin_nm"].append(wmin)
        out["offset_nm"].append(offset)
        out["ok"].append(ok)
        out["reason"].append(reason)

    out["line_nm"] = np.asarray(out["line_nm"], dtype=float)
    out["wmin_nm"] = np.asarray(out["wmin_nm"], dtype=float)
    out["offset_nm"] = np.asarray(out["offset_nm"], dtype=float)
    out["ok"] = np.asarray(out["ok"], dtype=bool)
    out["reason"] = np.asarray(out["reason"], dtype=object)

    return out


def build_analysis_results(soxs_shifted, tell_nm_air):
    # COLLECT ORIGINAL AND SHIFTED OFFSET MEASUREMENTS.
    results = []

    # STORE BOTH VERSIONS OF THE MEASUREMENTS PER SPECTRUM.
    for s in soxs_shifted:
        r_orig = measure_offsets_for_spectrum(s["wave"], s["flux"], tell_nm_air, plot=False, shift=0.0)

        r_shift = measure_offsets_for_spectrum(
            s["wave_shifted"], s["flux"], tell_nm_air, plot=False, shift=s["shift_nm"]
        )

        results.append(
            {
                "instrument": "SOXS",
                "name": s["file"].name,
                "date": s["date"],
                "file": s["file"],
                "shift_nm": s["shift_nm"],
                "reference_min_nm": s["shift_reference_min_nm"],
                "reference_hitran_nm": s["shift_reference_hitran_nm"],
                "results_original": r_orig,
                "results_shifted": r_shift,
            }
        )

    return results


# =========================
# READ XSHOOTER
# =========================
def read_xshooter_spectrum(filepath: Path):
    # READ ONE XSHOOTER SPECTRUM AND STANDARDIZE ITS FIELDS.
    with fits.open(filepath) as hdul:
        obs_name = "UNKNOWN"
        date_obs = "UNKNOWN"

        # SEARCH ALL HDUS FOR THE BEST METADATA.
        for hdu in hdul:
            header = hdu.header

            candidate_obs_name = get_header_value(
                header,
                [
                    "HIERARCH ESO OBS NAME",
                    "OBS.NAME",
                    "OBSNAME",
                    "OBJECT",
                    "HIERARCH ESO OBS TARG NAME",
                ],
            )

            if candidate_obs_name != "UNKNOWN":
                obs_name = candidate_obs_name

            candidate_date = get_date_obs(header)
            if candidate_date != "UNKNOWN":
                date_obs = candidate_date

            if obs_name != "UNKNOWN" and date_obs != "UNKNOWN":
                break

        if XSHOOTER_MATCH_STRING.upper() not in obs_name.upper():
            return None

        # IDENTIFY THE TABLE HDU AND ITS KEY COLUMNS.
        table_idx = find_flux_table_hdu(hdul)
        data = hdul[table_idx].data
        colnames = data.columns.names

        wave_col = find_column(colnames, ["WAVE", "WAVELENGTH", "LAMBDA", "WL"])
        flux_col = find_column(colnames, ["FLUX", "FLUX_REDUCED", "SCI_FLUX", "F"])

        if wave_col is None:
            raise KeyError(f"No wavelength column in {filepath}. Columns: {colnames}")

        if flux_col is None:
            raise KeyError(f"No flux column in {filepath}. Columns: {colnames}")

        wave = np.asarray(data[wave_col], dtype=float)
        flux = np.asarray(data[flux_col], dtype=float)

    # CONVERT MICRONS TO NM WHEN NEEDED.
    if np.nanmedian(wave) < 10:
        wave_nm = wave * 1000.0
    else:
        wave_nm = wave

    wave_nm, flux = clean_spectrum(wave_nm, flux)
    flux = normalize_flux(wave_nm, flux, filepath)

    return {
        "instrument": "XSHOOTER",
        "wave": wave_nm,
        "flux": flux,
        "date": date_obs,
        "file": filepath,
    }


def is_in_corron_path(path: Path):
    # EXCLUDE FILES UNDER _CORRON DIRECTORIES.
    return any(part.endswith("_CorrOn") for part in path.parts)


def find_column(colnames, candidates):
    # FIND A COLUMN BY EXACT OR PARTIAL NAME MATCH.
    upper_map = {c.upper(): c for c in colnames}

    for candidate in candidates:
        if candidate.upper() in upper_map:
            return upper_map[candidate.upper()]

    for c in colnames:
        cu = c.upper()
        for candidate in candidates:
            if candidate.upper() in cu:
                return c

    return None


def find_flux_table_hdu(hdul):
    # RETURN THE FIRST BINTABLE HDU INDEX.
    for idx, hdu in enumerate(hdul):
        if isinstance(hdu, fits.BinTableHDU):
            return idx
    raise ValueError("No BinTableHDU found")


def find_xshooter_spectra():
    # SCAN THE XSHOOTER TREE AND LOAD MATCHING FILES.
    files = sorted(XSHOOTER_ROOT.glob(f"**/{XSHOOTER_FILENAME}"))
    spectra = []

    # ENFORCE THE FILE LIMIT AND SKIP EXCLUDED PATHS.
    for f in files:
        if len(spectra) >= MAX_XSHOOTER_FILES:
            break

        if is_in_corron_path(f):
            continue

        try:
            spec = read_xshooter_spectrum(f)
        except Exception as e:
            print(f"[XSHOOTER ERROR] {f}")
            print(f"                 {e}")
            continue

        if spec is not None:
            spectra.append(spec)

    return spectra


# =========================
# PLOTS
# =========================
def apply_plot_style():
    # APPLY GLOBAL MATPLOTLIB STYLE FOR CONSISTENT FIGURES.
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
        }
    )


def make_axes_grid():
    # CREATE A FIXED-SIZE GRID OF SUBPLOTS.
    fig, axes = plt.subplots(
        N_ROWS,
        N_COLS,
        figsize=(4.5 * N_COLS, 3.0 * N_ROWS),
        squeeze=False,
    )
    return fig, axes.ravel()


def panel_title(res):
    # FORMAT PANEL TITLE WITH DATE, SHIFT, AND FILE NAME.
    return f"SOXS: {res['date']} | shift={res['shift_nm']:+.4f} nm\n" f"{res['name']}"


def add_offset_stats_box(ax, label, offsets):
    # ADD A SMALL STATS SUMMARY BOX TO A PANEL.
    offsets = np.asarray(offsets, dtype=float)
    offsets = offsets[np.isfinite(offsets)]
    offsets = offsets[np.abs(offsets) <= OFFSET_CLIP_NM]

    if offsets.size == 0:
        txt = f"{label}: N=0"
    else:
        txt = f"{label}: N={offsets.size}, " f"med={np.nanmedian(offsets):+.4f}, " f"σ={np.nanstd(offsets):.4f} nm"

    ax.text(
        0.02,
        0.98,
        txt,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75),
    )


def plot_offset_grid(all_results):
    # PLOT PER-SPECTRUM OFFSETS IN PAGINATED GRIDS.
    apply_plot_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # DRAW ONE PAGE OF PANELS AT A TIME.
    for page, page_results in chunks(all_results, PANELS_PER_FIG):
        fig, axes = make_axes_grid()

        # POPULATE EACH PANEL WITH ORIGINAL AND SHIFTED OFFSETS.
        for i, (ax, res) in enumerate(zip(axes, page_results)):
            row = i // N_COLS
            col = i % N_COLS

            is_left = col == 0
            is_bottom = row == N_ROWS - 1

            r_orig = res["results_original"]
            r_shift = res["results_shifted"]

            r_orig = clipped_good_mask(r_orig)
            r_shift = clipped_good_mask(r_shift)

            if len(r_orig.index):
                ax.plot(
                    r_orig["line_nm"],
                    r_orig["offset_nm"],
                    marker="x",
                    linestyle="-",
                    lw=0.8,
                    ms=3,
                    color=COLOR_SOXS,
                    alpha=0.75,
                    label="SOXS original - HITRAN",
                )

            if len(r_shift.index) and False:
                ax.plot(
                    r_shift["line_nm"],
                    r_shift["offset_nm"],
                    marker="o",
                    linestyle="-",
                    lw=0.8,
                    ms=3,
                    color=COLOR_SOXS_SHIFTED,
                    alpha=0.85,
                    label="SOXS shifted - HITRAN",
                )

            add_offset_stats_box(
                ax,
                "shifted",
                r_shift["offset_nm"],
            )

            ax.axhline(0.0, color="black", ls="--", lw=0.7)

            if XLIM_NM is not None:
                ax.set_xlim(*XLIM_NM)

            if YLIM_OFFSET_NM is not None:
                ax.set_ylim(*YLIM_OFFSET_NM)

            ax.set_title(panel_title(res).replace(".fits", ""), fontsize=6)

            if is_left:
                ax.set_ylabel("Offset to HITRAN [nm]")
            else:
                ax.set_ylabel("")
                ax.set_yticklabels([])

            if is_bottom:
                ax.set_xlabel("HITRAN wavelength [nm]")
            else:
                ax.set_xlabel("")
                ax.set_xticklabels([])

            ax.grid(alpha=0.2)

        # HIDE ANY UNUSED PANELS ON THE LAST PAGE.
        for ax in axes[len(page_results) :]:
            ax.axis("off")

        axes[0].legend(fontsize=7, loc="best")

        fig.suptitle(
            "SOXS telluric-line wavelength offsets w.r.t. HITRAN",
            fontsize=13,
        )

        fig.tight_layout()

        if SAVE_FIG:
            out = OUTDIR / f"{OUTFILE_PREFIX}_page{page:02d}.png"
            fig.savefig(out, dpi=200)
            print(f"[SAVED] {out}")

        if SHOW_FIG:
            plt.show()
        else:
            plt.close(fig)


def plot_summary_offset_vs_wavelength(all_results):
    # PLOT ALL OFFSET MEASUREMENTS IN ONE SUMMARY PANEL.
    apply_plot_style()

    fig, ax = plt.subplots(figsize=(13, 6))

    # OVERPLOT THE ORIGINAL AND SHIFTED MEASUREMENTS.
    for res in all_results:
        r_orig = res["results_original"]
        r_shift = res["results_shifted"]

        r_orig = clipped_good_mask(r_orig)
        r_shift = clipped_good_mask(r_shift)

        ax.plot(
            r_orig["line_nm"],
            r_orig["offset_nm"],
            marker="x",
            linestyle="none",
            ms=3,
            color=COLOR_SOXS,
            alpha=0.35,
        )

        ax.plot(
            r_shift["line_nm"],
            r_shift["offset_nm"],
            marker="o",
            linestyle="none",
            ms=3,
            color=COLOR_SOXS_SHIFTED,
            alpha=0.45,
        )

    ax.axhline(0.0, color="black", ls="--", lw=0.8)

    # BUILD A CLEAN LEGEND WITH PROXY ARTISTS.
    handles = [
        plt.Line2D(
            [],
            [],
            color=COLOR_SOXS,
            marker="x",
            ls="none",
            label="SOXS original",
        ),
        plt.Line2D(
            [],
            [],
            color=COLOR_SOXS_SHIFTED,
            marker="o",
            ls="none",
            label="SOXS shifted",
        ),
    ]

    ax.legend(handles=handles, frameon=False)
    ax.set_xlabel("HITRAN air wavelength [nm]")
    ax.set_ylabel("Measured minimum − HITRAN [nm]")
    ax.set_title("SOXS telluric-line offsets")
    ax.grid(alpha=0.25)

    if XLIM_NM is not None:
        ax.set_xlim(*XLIM_NM)

    if YLIM_OFFSET_NM is not None:
        ax.set_ylim(*YLIM_OFFSET_NM)

    fig.tight_layout()

    if SAVE_FIG:
        OUTDIR.mkdir(parents=True, exist_ok=True)
        out = OUTDIR / f"{OUTFILE_PREFIX}_summary.png"
        fig.savefig(out, dpi=200)
        print(f"[SAVED] {out}")

    if SHOW_FIG:
        plt.show()
    else:
        plt.close(fig)


def show_spectra_and_reference_line(soxs_spectra, xsh_spectra, tell_nm_air, plot=True):
    """PLOT AN EXAMPLE SOXS SPECTRUM AND MARK HITRAN REFERENCE LINES.

    Parameters
    - soxs_spectra: LIST OF SOXS SPECTRA DICTS
    - xsh_spectra: LIST OF XSHOOTER SPECTRA DICTS
    - tell_nm_air: HITRAN REFERENCE LINES
    - plot: IF TRUE SHOW FIGURE, ELSE CLOSE IT
    """

    apply_plot_style()
    # CREATE THE OVERVIEW AND ZOOM PANELS.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=False)

    # DEFINE THE ZOOM WINDOW.
    x_zoom_min = 1750.0
    x_zoom_max = 1752.0

    # DRAW ALL HITRAN REFERENCE LINES.
    for i, tn in enumerate(np.atleast_1d(tell_nm_air)):
        if i == 0:
            ax1.axvline(tn, color="tab:green", linestyle="-", lw=1.5, label="HITRAN line", alpha=0.3)
            ax2.axvline(tn, color="tab:green", linestyle="-", lw=1.5, alpha=0.3)
        else:
            ax1.axvline(tn, color="tab:green", linestyle="-", lw=1.5, alpha=0.3)
            ax2.axvline(tn, color="tab:green", linestyle="-", lw=1.5, alpha=0.3)

    # OVERPLOT THE XSHOOTER SPECTRA.
    for idx, spec in enumerate(xsh_spectra):
        wave = np.asarray(spec.get("wave", []), dtype=float)
        flux = np.asarray(spec.get("flux", []), dtype=float)

        label = "XSHOOTER spectrum" if idx == 0 else None
        ax1.plot(wave, flux, color="tab:orange", lw=0.8, label=label, alpha=0.5)
        ax2.plot(wave, flux, color="tab:orange", lw=0.8, alpha=0.5)

    stare = False
    nod = False
    # OVERPLOT THE SOXS SPECTRA BY OBSERVING MODE.
    for spec in soxs_spectra:
        wave = np.asarray(spec.get("wave", []), dtype=float)
        flux = np.asarray(spec.get("flux", []), dtype=float)
        mode = spec.get("mode")

        if mode == "stare":
            color = "tab:blue"
            if not stare:
                label = "SOXS stare"
                stare = True
            else:
                label = None
        else:
            color = "tab:red"
            if not nod:
                label = "SOXS nodding"
                nod = True
            else:
                label = None

        ax1.plot(wave, flux, color=color, lw=0.8, label=label)
        ax2.plot(wave, flux, color=color, lw=0.8)

    # HIGHLIGHT THE ZOOM WINDOW ON THE OVERVIEW PANEL.
    ax1.axvspan(x_zoom_min, x_zoom_max, color="lightgrey", alpha=0.5, zorder=0)

    # MATCH THE ZOOM PANEL BACKGROUND TO THE HIGHLIGHT.
    try:
        ax2.patch.set_facecolor("lightgrey")
        ax2.patch.set_alpha(0.5)
    except Exception:
        pass

    # SET LABELS AND VIEW LIMITS.
    ax2.set_xlabel("Wavelength [nm]")
    ax1.set_ylabel("Flux / Median")
    ax2.set_ylabel("Flux / Median")

    # CHOOSE A DATE LABEL FROM THE AVAILABLE SPECTRA.
    date_str = ""
    if len(soxs_spectra) > 0:
        date_str = soxs_spectra[0].get("date", "")
    elif len(xsh_spectra) > 0:
        date_str = xsh_spectra[0].get("date", "")
    ax1.set_title(f"Example SOXS spectrum: {date_str}")
    ax1.set_xlim(TOPPANEL_WAVE_RANGE[0], TOPPANEL_WAVE_RANGE[1])
    ax1.set_ylim(0, 2.0)
    ax2.set_xlim(x_zoom_min, x_zoom_max)
    ax2.set_ylim(0.2, 1.5)

    try:
        format_axes_no_scientific(ax1)
        format_axes_no_scientific(ax2)
    except Exception:
        pass

    ax1.legend(loc="best")
    fig.tight_layout()

    if plot:
        plt.show()
    else:
        plt.close(fig)


# =========================
# MAIN
# =========================
def main():
    # LOAD THE HITRAN LINE LIST AND CHOOSE THE SHIFT REFERENCE.
    tell_nm_air = load_hitran_air_nm(TELL_FILE)
    print(f"Loaded {tell_nm_air.size} HITRAN air lines.")

    hitran_ref_nm = nearest_hitran_reference_line(tell_nm_air)
    print(f"Reference HITRAN line for SOXS shift: {hitran_ref_nm:.6f} nm")

    # LOAD THE INPUT SOXS AND XSHOOTER SPECTRA.
    print("Searching SOXS spectra...")
    soxs = find_soxs_spectra()

    print("Searching XSHOOTER spectra...")
    xshooter = find_xshooter_spectra()

    print(f"Found {len(soxs)} SOXS spectra.")
    print(f"Found {len(xshooter)} XSHOOTER spectra.")

    # SHOW THE SPECTRA WITH THE REFERENCE-LINE CONTEXT.
    show_spectra_and_reference_line(soxs, xshooter, tell_nm_air, plot=True)

    # MEASURE AND APPLY THE GLOBAL SOXS SHIFTS.
    print("Applying SOXS shifts...")
    soxs_shifted = apply_soxs_shift(soxs, hitran_ref_nm, plot=False)

    # REPORT THE MEASURED SHIFT PER SPECTRUM.
    print("\n=== SOXS SHIFTS ===")
    for s in soxs_shifted:
        print(
            f"{s['date']}  "
            f"min={s['shift_reference_min_nm']:.6f} nm  "
            f"HITRAN={s['shift_reference_hitran_nm']:.6f} nm  "
            f"shift={s['shift_nm']:+.6f} nm  "
            f"{s['file'].name}"
        )

    allShifts = []
    for s in soxs_shifted:
        allShifts.append(s["shift_nm"])
        print(f"{s['shift_nm']:+.3f}")
    medianShift = np.median(allShifts)
    stdShift = np.std(allShifts)
    print(f"\nMedian shift: {medianShift:.3f} nm, Standard deviation: {stdShift:.3f} nm")

    # MEASURE ALL LINE OFFSETS AND GENERATE THE FIGURES.
    print("\nMeasuring offsets for all HITRAN lines...")
    all_results = build_analysis_results(soxs_shifted, tell_nm_air)

    if not all_results:
        print("No results to plot.")
        return

    plot_offset_grid(all_results)
    plot_summary_offset_vs_wavelength(all_results)

    print("Done.")


if __name__ == "__main__":
    main()
