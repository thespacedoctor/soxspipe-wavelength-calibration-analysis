#!/usr/bin/env python3
import warnings
from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
from scipy.optimize import OptimizeWarning, curve_fit

# =========================
# INPUTS
# =========================
SOXS_ROOT = Path("data/soxs_extractions_0.17.1")
SOXS_ROOT = Path("data/soxs_extractions_0.19.0")
# SOXS_ROOT = Path("data/tmp")
TELL_FILE = Path("reference/telluric_5_lines.dat")
XSHOOTER_ROOT = Path("data/xshooter_extractions")

SOXS_MATCH_STRING = "SOXS_LTT3218"
SOXS_MATCH_STRING = "SOXS_Feige110"
SOXS_MATCH_STRING = "SOXS_EG21"
SOXS_FLUX_COL = "FLUX_COUNTS"
XSHOOTER_MATCH_STRING = "Hip012377_TELL_SLT"
XSHOOTER_FILENAME = "SDP_TELL_SLIT_FLUX_MERGE1D_NIR.fits"
MAX_XSHOOTER_FILES = 20

NORMALIZE = True
NORM_MODE = "median"
NORM_WAVE_RANGE = (1671, 1676)

TOPPANEL_WAVE_RANGE = (1700, 1850)
OVERVIEW_ZOOM_WINDOW = (1750.0, 1752.0)

SEARCH_HALF_WINDOW_NM = 0.25
REFINE_HALF_WINDOW_NM = 0.15
MAX_OFFSET_TO_ACCEPT_NM = SEARCH_HALF_WINDOW_NM

MIN_POINTS_IN_WINDOW = 5
PARABOLA_REFINE = True
MIN_DEPTH = 0.05
MAX_CORE_RATIO = 0.1

# LOCAL CONTINUUM: A STRAIGHT LINE FITTED TO SIDEBANDS OUTSIDE THE LINE CORE, WHICH SCALES WITH THE LINE WIDTH.
DEFAULT_LINE_SIGMA_NM = 0.1
CORE_SIGMA_FACTOR = 2.5
NEIGHBOUR_MASK_SIGMA_FACTOR = 1.0
SIDEBAND_WIDTH_NM = 0.25
MIN_SIDEBAND_POINTS = 3
MIN_ONE_SIDED_POINTS = 3
MIN_POINTS_PER_SIDE_AFTER_REJECTION = 2
CONTINUUM_REJECT_SIGMA = 2.0
CONTINUUM_MIN_SIGMA_FRACTION = 1e-3
CONTINUUM_TOLERANCE = 1e-9

# LINE WIDTH: MEDIAN GAUSSIAN SIGMA OF ISOLATED, MODERATELY DEEP LINES FROM A COARSE FIRST PASS.
MIN_WIDTH_LINES = 10
WIDTH_DEPTH_RANGE = (0.2, 0.8)
WIDTH_FIT_HALF_WINDOW_NM = 0.6
WIDTH_SIGMA_BOUNDS_NM = (0.02, 0.6)
MIN_WIDTH_FIT_POINTS = 8
MIN_GAUSSIAN_AMPLITUDE = 1e-6
SIGMA_BOUND_LOW_FACTOR = 1.01
SIGMA_BOUND_HIGH_FACTOR = 0.99
OTHER_LINE_TOLERANCE_NM = 1e-6

# X-SHOOTER REFERENCE VALIDATION OF THE HITRAN LINE LIST.
MIN_XSH_FRACTION = 0.5
MIN_XSH_SPECTRA = 3
XSH_MAX_OFFSET_NM = 0.05

# OUTLIER CLIPPING OF OFFSETS AGAINST A LOCAL MEDIAN TREND.
CLIP_TREND_WINDOW_NM = 20
CLIP_SIGMA = 3.0
MIN_LINES_FOR_TREND = 10
MIN_TREND_NEIGHBOURS = 3
MAD_TO_SIGMA = 1.4826

DIAGNOSTIC_CUTOUTS = True

N_ROWS = 3
N_COLS = 4
PANELS_PER_FIG = N_ROWS * N_COLS

XLIM_NM = None
YLIM_OFFSET_NM = (-0.3, 0.3)

SAVE_FIG = True
SHOW_FIG = True

OUTDIR = SOXS_ROOT / "plots_soxs_telluric_offsets"
OUTFILE_PREFIX = f"{SOXS_MATCH_STRING}_telluric_offsets"

RESULT_COLUMNS = [
    "wmin_nm",
    "offset_nm",
    "depth",
    "core_flux",
    "continuum",
    "core_ratio",
    "continuum_mode",
    "parabola_fallback",
    "ok",
    "reason",
]
LINES_TABLE_COLUMNS = [
    "line_nm",
    "wmin_nm",
    "offset_nm",
    "depth",
    "core_flux",
    "continuum",
    "core_ratio",
    "continuum_mode",
    "parabola_fallback",
    "reason",
    "clipped",
]
CONTINUUM_MODES = ["two_sided", "one_sided"]
LOADER_ERRORS = (OSError, KeyError, ValueError, TypeError)
REASONS = ["ok", "too_few_points", "edge_minimum", "no_continuum", "too_shallow", "saturated", "offset_too_large"]

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


def try_read(reader: Callable[[Path], dict | None], filepath: Path, instrument: str) -> tuple[dict | None, dict | None]:
    # READ ONE FILE, RECORDING EXPECTED FITS/DATA ERRORS INSTEAD OF RAISING THEM.
    try:
        return reader(filepath), None
    except LOADER_ERRORS as e:
        error = f"{type(e).__name__}: {e}"
        print(f"[{instrument} ERROR] {filepath}")
        print(f"             {error}")
        return None, {"instrument": instrument, "file": filepath, "error": error}


def find_soxs_spectra() -> tuple[list, list]:
    # SCAN THE SOXS DIRECTORY AND LOAD MATCHING FILES; RETURN (SPECTRA, FAILED FILES).
    files = sorted(SOXS_ROOT.glob("*_EXTRACTED_ME*.fits"))
    spectra, failures = [], []

    # SKIP KNOWN BAD FILES AND KEEP ONLY READABLE SPECTRA.
    for f in files:
        if "NOTFLAT" in f.name or "ERROR" in f.name:
            continue

        spec, failure = try_read(read_soxs_spectrum, f, "SOXS")

        if failure is not None:
            failures.append(failure)
        elif spec is not None:
            spectra.append(spec)

    return spectra, failures


# =========================
# MINIMUM MEASUREMENT
# =========================
def fit_local_parabola_minimum(x, y) -> tuple:
    # FIT A LOCAL PARABOLA AND RETURN ITS VERTEX AS THE REFINED MINIMUM.
    good = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x, dtype=float)[good]
    y = np.asarray(y, dtype=float)[good]

    if x.size < 3:
        return np.nan, np.nan, None, None

    # TREAT A RANK-DEFICIENT FIT (RankWarning) AS A FAILED FIT.
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", np.exceptions.RankWarning)
            a, b, c0 = np.polyfit(x, y, 2)
    except (np.linalg.LinAlgError, ValueError, np.exceptions.RankWarning):
        return np.nan, np.nan, None, None

    if not np.isfinite(a) or a <= 0:
        return np.nan, np.nan, None, None

    x0 = -b / (2.0 * a)
    y0 = a * x0**2 + b * x0 + c0

    xfit = np.linspace(np.nanmin(x), np.nanmax(x), 50)
    yfit = a * xfit**2 + b * xfit + c0
    return x0, y0, xfit, yfit


def refine_minimum_position(w, f, i0: int) -> tuple[float, bool]:
    # REFINE THE DISCRETE MINIMUM WITH A 3-POINT PARABOLA; RETURN (POSITION, FELL_BACK_TO_SAMPLE).
    xloc = w[i0 - 1 : i0 + 2]
    yloc = f[i0 - 1 : i0 + 2]
    xfit, _, _, _ = fit_local_parabola_minimum(xloc, yloc)

    if np.isfinite(xfit) and xloc.min() <= xfit <= xloc.max():
        return float(xfit), False
    return float(w[i0]), True


def classify_minimum(
    i0: int,
    n_points: int,
    continuum: float,
    core_ratio: float,
    offset_nm: float,
    min_depth: float | None,
    max_core_ratio: float | None,
    max_offset_nm: float,
) -> str:
    # RETURN THE FIRST FAILED CHECK, OR "ok" WHEN THE MINIMUM IS A TRUSTWORTHY LINE.
    # THE ORDER MATTERS: A LATER CHECK ONLY RUNS WHEN THE EARLIER ONES PASS, AND THE DEPTH AND SATURATION
    # CHECKS NEED A CONTINUUM, SO NO_CONTINUUM MUST COME BEFORE THEM. EDGE_MINIMUM NEEDS NO CONTINUUM.
    if i0 == 0 or i0 == n_points - 1:
        return "edge_minimum"

    if not np.isfinite(continuum):
        return "no_continuum"

    if min_depth is not None and 1.0 - core_ratio < min_depth:
        return "too_shallow"

    if max_core_ratio is not None and core_ratio < max_core_ratio:
        return "saturated"

    if not np.isfinite(offset_nm) or abs(offset_nm) >= max_offset_nm:
        return "offset_too_large"

    return "ok"


def core_ratio(core_flux: float, continuum: float) -> float:
    # LOWEST FLUX DIVIDED BY THE CONTINUUM (NAN WITHOUT A CONTINUUM).
    return core_flux / continuum if np.isfinite(continuum) else np.nan


def minimum_result(
    reason: str,
    wmin_nm: float = np.nan,
    offset_nm: float = np.nan,
    core_flux: float = np.nan,
    continuum: float = np.nan,
    parabola_fallback: bool = False,
    continuum_mode: str = "none",
) -> dict:
    # ASSEMBLE ONE PER-LINE MEASUREMENT RECORD; DEPTH AND CORE RATIO FOLLOW FROM THE CORE FLUX AND CONTINUUM.
    coreRatio = core_ratio(core_flux, continuum)

    return {
        "ok": reason == "ok",
        "wmin_nm": wmin_nm,
        "offset_nm": offset_nm,
        "depth": 1.0 - coreRatio,
        "core_flux": core_flux,
        "continuum": continuum,
        "core_ratio": coreRatio,
        "continuum_mode": continuum_mode,
        "parabola_fallback": parabola_fallback,
        "reason": reason,
    }


def extract_window(wave_nm, flux, line_nm: float, half_window_nm: float) -> tuple:
    # RETURN THE FINITE SAMPLES WITHIN +/- HALF_WINDOW_NM OF THE LINE.
    inWindow = (
        np.isfinite(wave_nm)
        & np.isfinite(flux)
        & (wave_nm >= line_nm - half_window_nm)
        & (wave_nm <= line_nm + half_window_nm)
    )
    return wave_nm[inWindow], flux[inWindow]


def has_both_sides(x: np.ndarray, min_points: int, min_per_side: int = 1) -> bool:
    # TRUE WHEN THERE ARE ENOUGH POINTS WITH ENOUGH ON EACH SIDE OF THE CENTRE.
    return x.size >= min_points and int((x < 0).sum()) >= min_per_side and int((x > 0).sum()) >= min_per_side


def mask_near_lines(wave_nm: np.ndarray, lines_nm: np.ndarray | None, half_width_nm: float) -> np.ndarray:
    # TRUE FOR PIXELS FARTHER THAN HALF_WIDTH_NM FROM EVERY LINE IN LINES_NM.
    keep = np.ones(wave_nm.size, dtype=bool)

    if lines_nm is None or wave_nm.size == 0:
        return keep

    lines_nm = np.asarray(lines_nm, dtype=float)
    nearby = lines_nm[(lines_nm >= wave_nm.min() - half_width_nm) & (lines_nm <= wave_nm.max() + half_width_nm)]

    if nearby.size:
        keep = np.abs(wave_nm[:, None] - nearby[None, :]).min(axis=1) > half_width_nm

    return keep


def select_sideband_points(
    wave_nm,
    flux,
    center_nm: float,
    core_half_width_nm: float,
    sideband_width_nm: float = SIDEBAND_WIDTH_NM,
    masked_lines_nm: np.ndarray | None = None,
    mask_half_width_nm: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    # SIDEBAND PIXELS (OFFSET FROM CENTER_NM, FLUX): FROM THE CORE EDGE OUTWARD BY SIDEBAND_WIDTH_NM ON EACH SIDE.
    wave_nm = np.asarray(wave_nm, dtype=float)
    y = np.asarray(flux, dtype=float)
    x = wave_nm - center_nm
    maskHalf = core_half_width_nm if mask_half_width_nm is None else mask_half_width_nm

    # DROP NON-FINITE PIXELS, THE CORE, THE OUTER EDGE, AND PIXELS CLOSE TO OTHER HITRAN LINES.
    usable = np.isfinite(x) & np.isfinite(y) & (np.abs(x) > core_half_width_nm)
    usable &= np.abs(x) <= core_half_width_nm + sideband_width_nm
    usable &= mask_near_lines(wave_nm, masked_lines_nm, maskHalf)

    return x[usable], y[usable]


def fit_line_with_rejection(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    # STRAIGHT-LINE FIT WITH ONE REJECTION PASS; RETURN (SLOPE, INTERCEPT AT X = 0).
    slope, intercept = np.polyfit(x, y, 1)

    # DROP POINTS FAR BELOW THE FIT (BLENDED ABSORPTION) AND REFIT IF ENOUGH REMAIN ON BOTH SIDES.
    resid = y - (intercept + slope * x)
    madSigma = MAD_TO_SIGMA * np.median(np.abs(resid - np.median(resid)))
    sigma = max(madSigma, CONTINUUM_MIN_SIGMA_FRACTION * abs(intercept))
    keep = resid >= -(CONTINUUM_REJECT_SIGMA * sigma + CONTINUUM_TOLERANCE * np.median(np.abs(y)))
    if not keep.all() and has_both_sides(x[keep], MIN_SIDEBAND_POINTS, MIN_POINTS_PER_SIDE_AFTER_REJECTION):
        slope, intercept = np.polyfit(x[keep], y[keep], 1)

    return float(slope), float(intercept)


def one_sided_continuum(y: np.ndarray) -> float:
    # FLAT CONTINUUM FROM ONE SIDEBAND: THE MEDIAN AFTER ONE PASS THAT DROPS POINTS FAR BELOW IT.
    median = float(np.median(y))
    sigma = max(MAD_TO_SIGMA * np.median(np.abs(y - median)), CONTINUUM_MIN_SIGMA_FRACTION * abs(median))
    keep = y >= median - (CONTINUUM_REJECT_SIGMA * sigma + CONTINUUM_TOLERANCE * abs(median))

    if not keep.all() and keep.sum() >= MIN_ONE_SIDED_POINTS:
        return float(np.median(y[keep]))

    return median


def sideband_continuum(x: np.ndarray, y: np.ndarray) -> tuple[float, str]:
    # CONTINUUM AT X = 0 AND ITS MODE: ONE-SIDED WHEN ONE SIDE HAS 0-1 POINTS AND THE OTHER ENOUGH, ELSE TWO-SIDED.
    left, right = y[x < 0], y[x > 0]
    short, long = sorted((left, right), key=len)

    if short.size <= 1 and long.size >= MIN_ONE_SIDED_POINTS:
        return one_sided_continuum(long), "one_sided"

    if has_both_sides(x, MIN_SIDEBAND_POINTS):
        return fit_line_with_rejection(x, y)[1], "two_sided"

    return np.nan, "none"


def continuum_at_minimum(
    wave_nm,
    flux,
    wmin_nm: float,
    core_half_width_nm: float,
    sideband_width_nm: float = SIDEBAND_WIDTH_NM,
    masked_lines_nm: np.ndarray | None = None,
    mask_half_width_nm: float | None = None,
) -> tuple[float, str]:
    # SIDEBAND CONTINUUM AND ITS MODE AT THE MINIMUM: THE CORE AND THE SIDEBANDS ARE BOTH CENTRED ON THE MINIMUM.
    reach = core_half_width_nm + sideband_width_nm
    w, f = extract_window(wave_nm, flux, wmin_nm, reach)
    x, y = select_sideband_points(
        w, f, wmin_nm, core_half_width_nm, sideband_width_nm, masked_lines_nm, mask_half_width_nm
    )
    continuum, mode = sideband_continuum(x, y)

    # A NAN CONTINUUM FAILS THE > 0 TEST, SO IT ALSO MAPS TO NONE.
    return (continuum, mode) if continuum > 0 else (np.nan, "none")


def find_local_minimum_near_line(
    wave_nm,
    flux,
    line_nm: float,
    half_window_nm: float = SEARCH_HALF_WINDOW_NM,
    min_points: int = MIN_POINTS_IN_WINDOW,
    parabola_refine: bool = PARABOLA_REFINE,
    min_depth: float | None = MIN_DEPTH,
    max_core_ratio: float | None = MAX_CORE_RATIO,
    max_offset_nm: float = MAX_OFFSET_TO_ACCEPT_NM,
    line_sigma_nm: float = DEFAULT_LINE_SIGMA_NM,
    sideband_width_nm: float = SIDEBAND_WIDTH_NM,
    other_lines_nm: np.ndarray | None = None,
) -> dict:
    # REQUIRE ENOUGH SAMPLES BEFORE LOOKING FOR A MINIMUM.
    w, f = extract_window(wave_nm, flux, line_nm, half_window_nm)
    if w.size < min_points:
        return minimum_result("too_few_points")

    # START FROM THE LOWEST DISCRETE SAMPLE.
    i0 = int(np.nanargmin(f))
    coreFlux = float(f[i0])

    # REFINE ONLY INTERIOR MINIMA.
    wmin, isFallback = float(w[i0]), False
    if parabola_refine and 1 <= i0 <= w.size - 2:
        wmin, isFallback = refine_minimum_position(w, f, i0)

    # MEASURE THE CORE AGAINST THE LOCAL CONTINUUM; THE CORE IS A FIXED NUMBER OF LINE SIGMAS WIDE.
    coreHalfWidth = CORE_SIGMA_FACTOR * line_sigma_nm
    maskHalfWidth = NEIGHBOUR_MASK_SIGMA_FACTOR * line_sigma_nm
    continuum, mode = continuum_at_minimum(
        wave_nm, flux, wmin, coreHalfWidth, sideband_width_nm, other_lines_nm, maskHalfWidth
    )
    offset = wmin - line_nm
    coreRatio = core_ratio(coreFlux, continuum)
    reason = classify_minimum(i0, w.size, continuum, coreRatio, offset, min_depth, max_core_ratio, max_offset_nm)

    return minimum_result(reason, wmin, offset, coreFlux, continuum, isFallback, mode)


# =========================
# LINE WIDTH
# =========================
def gaussian_with_slope(u, c0, c1, amplitude, centre, sigma):
    # LINEAR CONTINUUM MINUS A GAUSSIAN ABSORPTION LINE.
    return c0 + c1 * u - amplitude * np.exp(-0.5 * ((u - centre) / sigma) ** 2)


def fit_gaussian_line(wave_nm, flux, wmin_nm: float, half_window_nm: float = WIDTH_FIT_HALF_WINDOW_NM) -> dict | None:
    # FIT A GAUSSIAN PLUS LINEAR CONTINUUM AROUND THE MINIMUM; RETURN SIGMA AND FRACTIONAL DEPTH, OR NONE.
    w, f = extract_window(wave_nm, flux, wmin_nm, half_window_nm)

    if w.size < MIN_WIDTH_FIT_POINTS:
        return None

    u = w - wmin_nm
    level = float(np.percentile(f, 90))
    lowSigma, highSigma = WIDTH_SIGMA_BOUNDS_NM
    p0 = [level, 0.0, max(level - float(f.min()), MIN_GAUSSIAN_AMPLITUDE), 0.0, DEFAULT_LINE_SIGMA_NM]
    bounds = ([-np.inf, -np.inf, 0.0, -0.1, lowSigma], [np.inf, np.inf, np.inf, 0.1, highSigma])

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, _ = curve_fit(gaussian_with_slope, u, f, p0=p0, bounds=bounds, maxfev=2000)
    except (RuntimeError, ValueError):
        return None

    c0, c1, amplitude, centre, sigma = popt
    continuum = c0 + c1 * centre

    # A SIGMA AT A BOUND OR A NON-POSITIVE CONTINUUM MEANS THE FIT DID NOT FIND A LINE.
    atBound = sigma <= lowSigma * SIGMA_BOUND_LOW_FACTOR or sigma >= highSigma * SIGMA_BOUND_HIGH_FACTOR
    if continuum <= 0 or atBound:
        return None

    return {"sigma_nm": float(sigma), "depth": float(amplitude / continuum)}


def estimate_line_sigma(
    wave_nm,
    flux,
    lines_nm,
    all_lines_nm: np.ndarray | None = None,
    min_lines: int = MIN_WIDTH_LINES,
    depth_range: tuple[float, float] = WIDTH_DEPTH_RANGE,
) -> dict:
    # MEDIAN GAUSSIAN SIGMA OF ISOLATED, MODERATELY DEEP, TWO-SIDED LINES ACCEPTED IN A COARSE PASS WITH THE DEFAULT
    # WIDTH. THE COARSE PASS HAS NO OFFSET ESTIMATE YET, SO ITS NEIGHBOUR MASK USES THE HITRAN POSITIONS UNSHIFTED.
    # ALSO RETURNS THE COARSE MEDIAN OFFSET, WHICH LATER PASSES USE TO ALIGN THE NEIGHBOUR MASK.
    allLines = np.asarray(lines_nm if all_lines_nm is None else all_lines_nm, dtype=float)
    coarse = measure_offsets_for_spectrum(wave_nm, flux, lines_nm, mask_lines_nm=allLines)
    accepted = coarse.loc[coarse["ok"]]
    medianOffset = clipped_median(accepted["offset_nm"])
    sigmas, nFailed = [], 0

    # A ONE-SIDED CONTINUUM IS FLAT, SO SLOPED LINES WOULD BIAS THE WIDTH: LEAVE THOSE LINES OUT.
    for row in accepted.loc[accepted["continuum_mode"] == "two_sided"].itertuples():
        others = allLines[np.abs(allLines - row.line_nm) > OTHER_LINE_TOLERANCE_NM]

        # SKIP LINES WITH ANOTHER HITRAN LINE INSIDE THE FIT WINDOW.
        if np.any(np.abs(others - row.line_nm) <= WIDTH_FIT_HALF_WINDOW_NM):
            continue

        fit = fit_gaussian_line(wave_nm, flux, row.wmin_nm)
        if fit is None:
            nFailed += 1
        elif depth_range[0] <= fit["depth"] <= depth_range[1]:
            sigmas.append(fit["sigma_nm"])

    sigmaMad = MAD_TO_SIGMA * float(np.median(np.abs(sigmas - np.median(sigmas)))) if sigmas else np.nan
    sigma = float(np.median(sigmas)) if len(sigmas) >= min_lines else DEFAULT_LINE_SIGMA_NM

    return {
        "sigma_nm": sigma,
        "n_fits": len(sigmas),
        "n_failed_fits": nFailed,
        "sigma_mad_nm": sigmaMad,
        "median_offset_nm": medianOffset if np.isfinite(medianOffset) else 0.0,
    }


# =========================
# STATS AND CLIPPING
# =========================
def robust_stats(values: np.ndarray | list) -> dict:
    # RETURN N, MEDIAN AND 1.4826 * MAD OF THE FINITE VALUES.
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return {"n": 0, "median": np.nan, "sigma": np.nan}

    median = float(np.median(values))
    sigma = float(MAD_TO_SIGMA * np.median(np.abs(values - median)))

    return {"n": int(values.size), "median": median, "sigma": sigma}


def clipped_median(values: np.ndarray | list, sigma: float = CLIP_SIGMA) -> float:
    # RETURN THE SIGMA-CLIPPED MEDIAN OF THE FINITE VALUES (NAN WHEN EMPTY).
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.nan

    _, median, _ = sigma_clipped_stats(values, sigma=sigma, maxiters=5, cenfunc="median", stdfunc="mad_std")

    return float(median)


def local_median_trend(lines_nm: np.ndarray, offsets: np.ndarray, window_nm: float) -> np.ndarray:
    # MEDIAN OFFSET OF THE OTHER LINES WITHIN +/- HALF THE WINDOW (WAVELENGTH, NOT INDEX).
    lines_nm = np.asarray(lines_nm, dtype=float)
    offsets = np.asarray(offsets, dtype=float)
    near = (np.abs(lines_nm[:, None] - lines_nm[None, :]) <= 0.5 * window_nm) & ~np.eye(lines_nm.size, dtype=bool)
    globalMedian = np.median(offsets)

    # TOO FEW LOCAL NEIGHBOURS FALL BACK TO THE GLOBAL MEDIAN.
    return np.array([np.median(offsets[row]) if row.sum() >= MIN_TREND_NEIGHBOURS else globalMedian for row in near])


def flag_outliers(
    df: pd.DataFrame, trend_window_nm: float = CLIP_TREND_WINDOW_NM, sigma: float = CLIP_SIGMA
) -> np.ndarray:
    # FLAG ACCEPTED LINES THAT DEVIATE FROM THE LOCAL TREND; REJECTED ROWS ARE NEVER FLAGGED.
    outliers = np.zeros(len(df), dtype=bool)
    okIdx = np.flatnonzero(df["ok"].to_numpy(dtype=bool))

    if okIdx.size == 0:
        return outliers

    lines = df["line_nm"].to_numpy(dtype=float)[okIdx]
    offsets = df["offset_nm"].to_numpy(dtype=float)[okIdx]

    # FALL BACK TO A GLOBAL MEDIAN WHEN TOO FEW LINES DEFINE A TREND.
    if okIdx.size >= MIN_LINES_FOR_TREND:
        trend = local_median_trend(lines, offsets, trend_window_nm)
    else:
        trend = np.full(offsets.size, np.median(offsets))

    clipped = sigma_clip(offsets - trend, sigma=sigma, maxiters=5, cenfunc="median", stdfunc="mad_std")
    outliers[okIdx] = np.ma.getmaskarray(clipped)

    return outliers


def add_clip_flags(df: pd.DataFrame) -> pd.DataFrame:
    # RETURN A COPY OF THE RESULTS WITH A BOOLEAN "clipped" COLUMN.
    return df.assign(clipped=flag_outliers(df))


def good_mask(df: pd.DataFrame) -> pd.Series:
    # ACCEPTED AND NOT CLIPPED.
    return df["ok"].astype(bool) & ~df["clipped"].astype(bool)


def good_offsets(df: pd.DataFrame) -> np.ndarray:
    # OFFSETS OF ACCEPTED, UNCLIPPED LINES.
    return df.loc[good_mask(df), "offset_nm"].to_numpy(dtype=float)


def reject_counts(df: pd.DataFrame) -> dict:
    # COUNT LINES BY OUTCOME REASON, INCLUDING ZEROS FOR REASONS THAT NEVER OCCURRED.
    counts = df["reason"].value_counts()
    return {reason: int(counts.get(reason, 0)) for reason in REASONS}


# =========================
# GLOBAL SHIFT
# =========================
def derive_global_shift(
    wave_nm,
    flux,
    lines_nm,
    search_half_window_nm: float = SEARCH_HALF_WINDOW_NM,
    refine_half_window_nm: float = REFINE_HALF_WINDOW_NM,
    line_sigma_nm: float = DEFAULT_LINE_SIGMA_NM,
    mask_lines_nm: np.ndarray | None = None,
    mask_offset_nm: float = 0.0,
) -> dict:
    # COARSE PASS: SHIFT = -(SIGMA-CLIPPED MEDIAN OF ACCEPTED PER-LINE OFFSETS); MASK_OFFSET_NM ALIGNS THE MASK.
    coarse = measure_offsets_for_spectrum(
        wave_nm, flux, lines_nm, search_half_window_nm, line_sigma_nm, mask_lines_nm, mask_offset_nm
    )
    coarseShift = -clipped_median(coarse.loc[coarse["ok"], "offset_nm"])

    if not np.isfinite(coarseShift):
        raise ValueError("No accepted line offsets to derive a global shift")

    # REFINE PASS: REMEASURE ON THE SHIFTED AXIS (DATA LINES NOW AT HITRAN POSITIONS) WITH A NARROWER WINDOW.
    refined = add_clip_flags(
        measure_offsets_for_spectrum(
            wave_nm + coarseShift, flux, lines_nm, refine_half_window_nm, line_sigma_nm, mask_lines_nm
        )
    )
    good = good_offsets(refined)

    if good.size == 0:
        raise ValueError("No accepted line offsets after the coarse shift")

    return {
        "shift_nm": float(coarseShift - clipped_median(good)),
        "n_lines": int(good.size),
        "robust_sigma_nm": robust_stats(good)["sigma"],
    }


# =========================
# MEASURE ALL HITRAN OFFSETS
# =========================
def measure_offsets_for_spectrum(
    wave_nm,
    flux,
    lines_nm,
    half_window_nm: float = SEARCH_HALF_WINDOW_NM,
    line_sigma_nm: float = DEFAULT_LINE_SIGMA_NM,
    mask_lines_nm: np.ndarray | None = None,
    mask_offset_nm: float = 0.0,
) -> pd.DataFrame:
    # MEASURE OFFSETS FOR ALL HITRAN LINES INSIDE THE SPECTRAL COVERAGE; MASK_LINES_NM IS THE FULL LINE LIST.
    # MASK_OFFSET_NM IS THE KNOWN OFFSET OF THE DATA FROM HITRAN: THE MASKED POSITIONS ARE HITRAN + MASK_OFFSET_NM.
    wave_nm = np.asarray(wave_nm, dtype=float)
    flux = np.asarray(flux, dtype=float)
    lines_nm = np.asarray(lines_nm, dtype=float)
    mask_lines_nm = None if mask_lines_nm is None else np.asarray(mask_lines_nm, dtype=float)

    if wave_nm.size == 0:
        lines_here = np.array([])
    else:
        lines_here = lines_nm[(lines_nm >= np.nanmin(wave_nm)) & (lines_nm <= np.nanmax(wave_nm))]

    # MEASURE EACH LINE INDEPENDENTLY.
    rows = []
    for line_nm in lines_here:
        others = None
        if mask_lines_nm is not None:
            others = mask_lines_nm[np.abs(mask_lines_nm - line_nm) > OTHER_LINE_TOLERANCE_NM] + mask_offset_nm

        result = find_local_minimum_near_line(
            wave_nm, flux, line_nm, half_window_nm=half_window_nm, line_sigma_nm=line_sigma_nm, other_lines_nm=others
        )
        rows.append({"line_nm": line_nm, **result})

    return pd.DataFrame(rows, columns=["line_nm"] + RESULT_COLUMNS).astype({"line_nm": float, "ok": bool})


# =========================
# REFERENCE VALIDATION
# =========================
def validate_reference_lines(
    results: list,
    lines_nm: np.ndarray,
    min_fraction: float = MIN_XSH_FRACTION,
    max_offset_nm: float = XSH_MAX_OFFSET_NM,
    min_spectra: int = MIN_XSH_SPECTRA,
) -> pd.DataFrame:
    # KEEP LINES THAT X-SHOOTER MEASURES CLEANLY (ACCEPTED, UNCLIPPED) NEAR ZERO IN ENOUGH SPECTRA.
    lines = np.asarray(lines_nm, dtype=float)
    table = pd.DataFrame({"line_nm": lines})

    # TOO FEW REFERENCE SPECTRA: VALIDATION IS UNAVAILABLE, NO LINE IS MARKED VALID.
    if len(results) < min_spectra:
        return table.assign(
            n_measured=0,
            n_accepted=0,
            fraction=np.nan,
            median_offset_nm=np.nan,
            is_valid=False,
            validation_available=False,
        )

    offsets = np.full((len(results), lines.size), np.nan)
    measured = np.zeros((len(results), lines.size), dtype=bool)
    for i, df in enumerate(results):
        byLine = df.set_index("line_nm")
        offsets[i] = byLine["offset_nm"].where(good_mask(byLine)).reindex(lines).to_numpy(dtype=float)
        measured[i] = pd.Index(lines).isin(byLine.index)

    accepted = np.isfinite(offsets) & (np.abs(offsets) < max_offset_nm)
    nMeasured = measured.sum(axis=0)
    nAccepted = accepted.sum(axis=0)
    fraction = np.divide(nAccepted, nMeasured, out=np.zeros(lines.size), where=nMeasured > 0)
    medians = pd.DataFrame(np.where(accepted, offsets, np.nan)).median(axis=0, skipna=True).to_numpy()

    return table.assign(
        n_measured=nMeasured,
        n_accepted=nAccepted,
        fraction=fraction,
        median_offset_nm=medians,
        is_valid=(nAccepted > 0) & (fraction >= min_fraction),
        validation_available=True,
    )


def select_measurement_lines(validation: pd.DataFrame) -> np.ndarray:
    # USE THE VALIDATED LINES, OR EVERY LINE WHEN VALIDATION IS UNAVAILABLE.
    if validation["validation_available"].all():
        return validation.loc[validation["is_valid"], "line_nm"].to_numpy(dtype=float)
    return validation["line_nm"].to_numpy(dtype=float)


def differential_offsets(df: pd.DataFrame, reference: pd.DataFrame) -> pd.Series:
    # OFFSET MINUS THE MEDIAN X-SHOOTER OFFSET OF THE SAME LINE (NAN WHEN NO REFERENCE).
    referenceByLine = reference.set_index("line_nm")["median_offset_nm"]
    return df["offset_nm"] - df["line_nm"].map(referenceByLine)


def measure_spectra(spectra: list, lines_nm: np.ndarray) -> tuple[list, list]:
    # ESTIMATE THE LINE WIDTH, THEN MEASURE AND CLIP-FLAG EVERY SPECTRUM ON ITS UNSHIFTED WAVELENGTH AXIS.
    results, sigmas = [], []

    for s in spectra:
        info = estimate_line_sigma(s["wave"], s["flux"], lines_nm)
        measured = measure_offsets_for_spectrum(
            s["wave"],
            s["flux"],
            lines_nm,
            line_sigma_nm=info["sigma_nm"],
            mask_lines_nm=lines_nm,
            mask_offset_nm=info["median_offset_nm"],
        )
        results.append(add_clip_flags(measured))
        sigmas.append(info)

    return results, sigmas


def build_analysis_results(
    soxs: list, lines_nm: np.ndarray, validation: pd.DataFrame, all_lines_nm: np.ndarray | None = None
) -> tuple[list, list]:
    # COLLECT ORIGINAL AND SHIFTED OFFSETS PER SOXS SPECTRUM; RETURN (RESULTS, FAILED SPECTRA).
    # ALL_LINES_NM IS THE FULL HITRAN LIST, USED TO MASK NEIGHBOURING LINES OUT OF THE SIDEBANDS.
    results, failures = [], []
    allLines = np.asarray(lines_nm if all_lines_nm is None else all_lines_nm, dtype=float)

    for s in soxs:
        sigmaInfo = estimate_line_sigma(s["wave"], s["flux"], lines_nm, allLines)
        sigma = sigmaInfo["sigma_nm"]

        try:
            shiftInfo = derive_global_shift(
                s["wave"],
                s["flux"],
                lines_nm,
                line_sigma_nm=sigma,
                mask_lines_nm=allLines,
                mask_offset_nm=sigmaInfo["median_offset_nm"],
            )
        except ValueError as e:
            print(f"[SOXS SHIFT ERROR] {s['file']}")
            print(f"                  {e}")
            failures.append({"instrument": "SOXS", "file": s["file"], "error": str(e)})
            continue

        # MEASURE THE SHIFTED SPECTRUM WITH THE SAME WINDOW AS THE SHIFT ESTIMATOR (REFINE PASS).
        # THE UNSHIFTED DATA LINES SIT AT HITRAN MINUS THE SHIFT, SO THE NEIGHBOUR MASK MOVES BY -SHIFT.
        original = add_clip_flags(
            measure_offsets_for_spectrum(
                s["wave"],
                s["flux"],
                lines_nm,
                line_sigma_nm=sigma,
                mask_lines_nm=allLines,
                mask_offset_nm=-shiftInfo["shift_nm"],
            )
        )
        shifted = add_clip_flags(
            measure_offsets_for_spectrum(
                s["wave"] + shiftInfo["shift_nm"], s["flux"], lines_nm, REFINE_HALF_WINDOW_NM, sigma, allLines
            )
        )

        results.append(
            {
                "instrument": "SOXS",
                "name": s["file"].name,
                "date": s["date"],
                "file": s["file"],
                "wave": s["wave"],
                "flux": s["flux"],
                "shift_nm": shiftInfo["shift_nm"],
                "shift_info": shiftInfo,
                "sigma_info": sigmaInfo,
                "results_original": original,
                "results_shifted": shifted,
                "differential_original": differential_offsets(original, validation),
                "differential_shifted": differential_offsets(shifted, validation),
            }
        )

    return results, failures


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


def find_xshooter_spectra() -> tuple[list, list]:
    # SCAN THE XSHOOTER TREE AND LOAD MATCHING FILES; RETURN (SPECTRA, FAILED FILES).
    files = sorted(XSHOOTER_ROOT.glob(f"**/{XSHOOTER_FILENAME}"))
    spectra, failures = [], []

    # ENFORCE THE FILE LIMIT AND SKIP EXCLUDED PATHS.
    for f in files:
        if len(spectra) >= MAX_XSHOOTER_FILES:
            break

        if is_in_corron_path(f):
            continue

        spec, failure = try_read(read_xshooter_spectrum, f, "XSHOOTER")

        if failure is not None:
            failures.append(failure)
        elif spec is not None:
            spectra.append(spec)

    return spectra, failures


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


def make_axes_grid(n_rows: int = N_ROWS, n_cols: int = N_COLS):
    # CREATE A FIXED-SIZE GRID OF SUBPLOTS.
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4.5 * n_cols, 3.0 * n_rows),
        squeeze=False,
    )
    return fig, axes.ravel()


def panel_title(res):
    # FORMAT PANEL TITLE WITH DATE, SHIFT, AND FILE NAME.
    return f"SOXS: {res['date']} | shift={res['shift_nm']:+.4f} nm\n" f"{res['name']}"


def add_offset_stats_box(ax, label: str, offsets, y_pos: float = 0.98, color: str = "black") -> None:
    # ADD A SMALL ROBUST-STATS SUMMARY BOX (MEDIAN AND 1.4826*MAD) TO A PANEL.
    stats = robust_stats(offsets)

    if stats["n"] == 0:
        txt = f"{label}: N=0"
    else:
        txt = f"{label}: N={stats['n']}, med={stats['median']:+.4f}, σ={stats['sigma']:.4f} nm"

    ax.text(
        0.02,
        y_pos,
        txt,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        color=color,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.75},
    )


def plot_offsets_series(
    ax, df: pd.DataFrame, marker: str, color: str, label: str | None, alpha: float = 0.8, linestyle: str = "-"
) -> None:
    # PLOT THE ACCEPTED, UNCLIPPED OFFSETS OF ONE MEASUREMENT SET.
    good = df.loc[good_mask(df)]
    if len(good.index):
        ax.plot(
            good["line_nm"],
            good["offset_nm"],
            marker=marker,
            linestyle=linestyle,
            lw=0.8,
            ms=3,
            color=color,
            alpha=alpha,
            label=label,
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

            plot_offsets_series(ax, r_orig, "x", COLOR_SOXS, "SOXS original - HITRAN", alpha=0.75)
            plot_offsets_series(ax, r_shift, "o", COLOR_SOXS_SHIFTED, "SOXS shifted - HITRAN", alpha=0.85)

            add_offset_stats_box(ax, "original", good_offsets(r_orig), y_pos=0.98, color=COLOR_SOXS)
            add_offset_stats_box(ax, "shifted", good_offsets(r_shift), y_pos=0.88, color=COLOR_SOXS_SHIFTED)

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
    # PLOT ALL OFFSETS (TOP) AND SOXS-MINUS-XSHOOTER DIFFERENTIALS (BOTTOM) IN ONE SUMMARY FIGURE.
    apply_plot_style()

    fig, (ax, axDiff) = plt.subplots(2, 1, figsize=(13, 9), sharex=True)

    # OVERPLOT THE ORIGINAL AND SHIFTED MEASUREMENTS.
    for res in all_results:
        plot_offsets_series(ax, res["results_original"], "x", COLOR_SOXS, None, alpha=0.35, linestyle="none")
        plot_offsets_series(ax, res["results_shifted"], "o", COLOR_SOXS_SHIFTED, None, alpha=0.45, linestyle="none")

        shifted = res["results_shifted"]
        good = good_mask(shifted)
        axDiff.plot(
            shifted.loc[good, "line_nm"],
            res["differential_shifted"][good],
            marker="o",
            linestyle="none",
            ms=3,
            color=COLOR_SOXS_SHIFTED,
            alpha=0.45,
        )

    ax.axhline(0.0, color="black", ls="--", lw=0.8)
    axDiff.axhline(0.0, color="black", ls="--", lw=0.8)

    # BUILD A CLEAN LEGEND WITH PROXY ARTISTS.
    handles = [
        plt.Line2D([], [], color=COLOR_SOXS, marker="x", ls="none", label="SOXS original"),
        plt.Line2D([], [], color=COLOR_SOXS_SHIFTED, marker="o", ls="none", label="SOXS shifted"),
    ]

    ax.legend(handles=handles, frameon=False)
    ax.set_ylabel("Measured minimum − HITRAN [nm]")
    ax.set_title("SOXS telluric-line offsets")
    axDiff.set_xlabel("HITRAN air wavelength [nm]")
    axDiff.set_ylabel("Shifted SOXS − X-shooter median [nm]")
    axDiff.set_title("Differential offsets (validated lines)")

    for a in (ax, axDiff):
        a.grid(alpha=0.25)
        if XLIM_NM is not None:
            a.set_xlim(*XLIM_NM)
        if YLIM_OFFSET_NM is not None:
            a.set_ylim(*YLIM_OFFSET_NM)

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
    x_zoom_min, x_zoom_max = OVERVIEW_ZOOM_WINDOW

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
    ax2.patch.set_facecolor("lightgrey")
    ax2.patch.set_alpha(0.5)

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

    format_axes_no_scientific(ax1)
    format_axes_no_scientific(ax2)

    ax1.legend(loc="best")
    fig.tight_layout()

    if plot:
        plt.show()
    else:
        plt.close(fig)


# =========================
# DIAGNOSTIC CUTOUTS
# =========================
def plot_line_cutout(
    ax, wave_nm: np.ndarray, flux: np.ndarray, row: dict, half_window_nm: float = SEARCH_HALF_WINDOW_NM
) -> None:
    # DRAW THE LOCAL FLUX AROUND ONE LINE WITH THE MEASURED MINIMUM AND THE REJECTION REASON.
    near = np.abs(wave_nm - row["line_nm"]) <= half_window_nm
    ax.plot(wave_nm[near], flux[near], "o-", ms=3, lw=1, color="tab:blue")
    ax.axvline(row["line_nm"], color="tab:orange", ls="--", lw=1)

    if np.isfinite(row["wmin_nm"]):
        color = "tab:green" if row["ok"] else "tab:red"
        ax.plot(row["wmin_nm"], np.interp(row["wmin_nm"], wave_nm[near], flux[near]), "x", ms=8, mew=2, color=color)

    status = row["reason"] if not row["ok"] else "clipped"
    detail = f"offset={row['offset_nm']:+.3f} depth={row['depth']:.2f} ratio={row['core_ratio']:.2f}"
    ax.set_title(f"{row['line_nm']:.3f} nm: {status}\n{detail}", fontsize=7)
    format_axes_no_scientific(ax)


def write_diagnostic_cutouts(
    wave_nm,
    flux,
    df: pd.DataFrame,
    path: Path,
    half_window_nm: float = SEARCH_HALF_WINDOW_NM,
    panels_per_page: int = PANELS_PER_FIG,
) -> int:
    # WRITE A MULTI-PAGE PDF OF CUTOUTS FOR REJECTED AND CLIPPED LINES; RETURN THE PANEL COUNT.
    wave_nm = np.asarray(wave_nm, dtype=float)
    flux = np.asarray(flux, dtype=float)

    # LINES WITH NO SAMPLES IN THE WINDOW HAVE NOTHING TO SHOW.
    flagged = df.loc[(~df["ok"].astype(bool) | df["clipped"].astype(bool)) & (df["reason"] != "too_few_points")]

    if flagged.empty:
        return 0

    apply_plot_style()
    nCols = min(N_COLS, panels_per_page)
    nRows = int(np.ceil(panels_per_page / nCols))
    path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(path) as pdf:
        for _, pageRows in chunks(list(flagged.to_dict("records")), panels_per_page):
            fig, axes = make_axes_grid(nRows, nCols)
            for ax, row in zip(axes, pageRows):
                plot_line_cutout(ax, wave_nm, flux, row, half_window_nm)
            for ax in axes[len(pageRows) :]:
                ax.axis("off")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    print(f"[SAVED] {path}")
    return len(flagged)


# =========================
# TABLES
# =========================
def build_lines_table(
    spectrum: str, instrument: str, version: str, df: pd.DataFrame, differential: pd.Series | None = None
) -> pd.DataFrame:
    # ONE ROW PER MEASURED LINE FOR THE LINES CSV.
    table = df[LINES_TABLE_COLUMNS].copy()
    table.insert(0, "version", version)
    table.insert(0, "instrument", instrument)
    table.insert(0, "spectrum", spectrum)
    table["differential_nm"] = np.nan if differential is None else differential.to_numpy()

    return table


def stats_columns(prefix: str, values: np.ndarray | list) -> dict:
    # MEDIAN AND ROBUST SIGMA COLUMNS FOR ONE SERIES.
    stats = robust_stats(values)
    return {f"median_{prefix}_nm": stats["median"], f"sigma_{prefix}_nm": stats["sigma"]}


def build_summary_row(
    spectrum: str,
    instrument: str,
    original: pd.DataFrame,
    shifted: pd.DataFrame | None,
    shift_info: dict | None,
    differential_original: pd.Series | None = None,
    differential_shifted: pd.Series | None = None,
    sigma_info: dict | None = None,
) -> dict:
    # ONE ROW PER SPECTRUM FOR THE SUMMARY CSV (SHIFT, LINE WIDTH, STATS AND REJECT COUNTS).
    row = {
        "spectrum": spectrum,
        "instrument": instrument,
        "status": "ok",
        "error": "",
        "line_sigma_nm": np.nan if sigma_info is None else sigma_info["sigma_nm"],
        "n_width_fits": 0 if sigma_info is None else sigma_info["n_fits"],
        "n_width_failed": 0 if sigma_info is None else sigma_info["n_failed_fits"],
        "line_sigma_mad_nm": np.nan if sigma_info is None else sigma_info["sigma_mad_nm"],
        "shift_nm": np.nan if shift_info is None else shift_info["shift_nm"],
        "n_lines": 0 if shift_info is None else shift_info["n_lines"],
        "n_measured": len(original),
        "n_ok": int(original["ok"].sum()),
        "n_good_original": len(good_offsets(original)),
        "n_good_shifted": np.nan if shifted is None else len(good_offsets(shifted)),
    }
    row.update(stats_columns("original", good_offsets(original)))
    row.update(stats_columns("shifted", [] if shifted is None else good_offsets(shifted)))

    diffSets = (("diff_original", original, differential_original), ("diff_shifted", shifted, differential_shifted))
    for name, df, diff in diffSets:
        values = [] if diff is None else diff[good_mask(df)].to_numpy()
        row.update(stats_columns(name, values))

    row.update({f"reason_{reason}": count for reason, count in reject_counts(original).items()})
    row.update({f"mode_{mode}": int((original["continuum_mode"] == mode).sum()) for mode in CONTINUUM_MODES})

    return row


def build_failed_row(instrument: str, failure: dict) -> dict:
    # ONE SUMMARY ROW FOR A FILE THAT COULD NOT BE READ OR MEASURED.
    file = Path(failure["file"])
    return {
        "spectrum": f"{file.parent.name}/{file.name}",
        "instrument": instrument,
        "status": "failed",
        "error": failure["error"],
    }


def write_tables(lines_table: pd.DataFrame, summary_table: pd.DataFrame, validation: pd.DataFrame) -> None:
    # WRITE THE LINES, SUMMARY AND REFERENCE-LINE VALIDATION CSVS.
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for suffix, table in (("lines", lines_table), ("summary", summary_table), ("reference_lines", validation)):
        out = OUTDIR / f"{OUTFILE_PREFIX}_{suffix}.csv"
        table.to_csv(out, index=False)
        print(f"[SAVED] {out}")


# =========================
# MAIN
# =========================
def xshooter_name(spec: dict) -> str:
    # X-SHOOTER FILE NAMES REPEAT ACROSS DIRECTORIES, SO PREFIX THE PARENT FOLDER.
    return f"{spec['file'].parent.name}/{spec['file'].name}"


def assemble_tables(
    all_results: list,
    xshooter: list,
    xsh_results: list,
    failures: list | None = None,
    xsh_sigmas: list | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # BUILD THE LINES AND SUMMARY TABLES FOR SOXS (ORIGINAL AND SHIFTED), X-SHOOTER (CONTROL) AND FAILED FILES.
    lines, summary = [], []

    for res in all_results:
        original, shifted = res["results_original"], res["results_shifted"]
        lines.append(build_lines_table(res["name"], "SOXS", "original", original, res["differential_original"]))
        lines.append(build_lines_table(res["name"], "SOXS", "shifted", shifted, res["differential_shifted"]))
        summary.append(
            build_summary_row(
                res["name"],
                "SOXS",
                original,
                shifted,
                res["shift_info"],
                res["differential_original"],
                res["differential_shifted"],
                res["sigma_info"],
            )
        )

    sigmas = xsh_sigmas or [None] * len(xshooter)
    for spec, df, sigmaInfo in zip(xshooter, xsh_results, sigmas):
        lines.append(build_lines_table(xshooter_name(spec), "XSHOOTER", "original", df))
        summary.append(build_summary_row(xshooter_name(spec), "XSHOOTER", df, None, None, sigma_info=sigmaInfo))

    summary.extend(build_failed_row(f["instrument"], f) for f in failures or [])

    return pd.concat(lines, ignore_index=True), pd.DataFrame(summary)


def print_report(
    all_results: list, xsh_results: list, validation: pd.DataFrame, failures: list, xsh_sigmas: list | None = None
) -> None:
    # PRINT THE X-SHOOTER CONTROL STATS, VALIDATION COUNTS, FAILED FILES AND PER-SPECTRUM SOXS RESULTS.
    if validation["validation_available"].all():
        nValid = int(validation["is_valid"].sum())
        print(f"\n{nValid} of {len(validation)} HITRAN lines validated against X-shooter.")
    else:
        print(
            f"\nWARNING: only {len(xsh_results)} usable X-shooter spectra (need {MIN_XSH_SPECTRA}). "
            "Line validation unavailable: measuring all HITRAN lines and skipping the differential."
        )

    if xsh_results:
        control = robust_stats(np.concatenate([good_offsets(df) for df in xsh_results]))
        print(
            f"X-shooter control (all lines): N={control['n']}, "
            f"median={control['median']:+.4f} nm, σ={control['sigma']:.4f} nm"
        )

    if xsh_sigmas:
        widths = [info["sigma_nm"] for info in xsh_sigmas]
        print(f"X-shooter line sigma: median={np.median(widths):.3f} nm (range {min(widths):.3f} to {max(widths):.3f})")

    if failures:
        print(f"\n=== {len(failures)} FAILED FILES ===")
        for failure in failures:
            print(f"{failure['instrument']}: {failure['file']}  {failure['error']}")

    print("\n=== SOXS GLOBAL SHIFTS ===")
    for res in all_results:
        orig = robust_stats(good_offsets(res["results_original"]))
        shifted = robust_stats(good_offsets(res["results_shifted"]))
        print(
            f"{res['date']}  sigma={res['sigma_info']['sigma_nm']:.3f} ({res['sigma_info']['n_fits']} fits)  "
            f"shift={res['shift_nm']:+.4f} nm  N={res['shift_info']['n_lines']}  "
            f"original: med={orig['median']:+.4f} σ={orig['sigma']:.4f}  "
            f"shifted: med={shifted['median']:+.4f} σ={shifted['sigma']:.4f}  {res['name']}"
        )

    shifts = [res["shift_nm"] for res in all_results]
    if shifts:
        print(f"\nMedian shift: {np.median(shifts):+.4f} nm, Standard deviation: {np.std(shifts):.4f} nm")


def main() -> None:
    # LOAD THE HITRAN LINE LIST AND THE INPUT SPECTRA.
    tell_nm_air = load_hitran_air_nm(TELL_FILE)
    print(f"Loaded {tell_nm_air.size} HITRAN air lines.")

    print("Searching SOXS spectra...")
    soxs, soxsFailures = find_soxs_spectra()

    print("Searching XSHOOTER spectra...")
    xshooter, xshFailures = find_xshooter_spectra()

    print(f"Found {len(soxs)} SOXS spectra.")
    print(f"Found {len(xshooter)} XSHOOTER spectra.")

    if SHOW_FIG:
        show_spectra_and_reference_line(soxs, xshooter, tell_nm_air, plot=True)

    # VALIDATE THE LINE LIST AGAINST THE X-SHOOTER REFERENCE SPECTRA.
    xsh_results, xsh_sigmas = measure_spectra(xshooter, tell_nm_air)
    validation = validate_reference_lines(xsh_results, tell_nm_air)
    measurement_lines = select_measurement_lines(validation)

    # MEASURE ALL LINE OFFSETS AND GLOBAL SHIFTS FOR EACH SOXS SPECTRUM.
    print("\nMeasuring offsets for HITRAN lines...")
    all_results, analysisFailures = build_analysis_results(soxs, measurement_lines, validation, tell_nm_air)
    failures = soxsFailures + xshFailures + analysisFailures
    print_report(all_results, xsh_results, validation, failures, xsh_sigmas)

    if not all_results:
        print("No results to plot.")
        return

    lines_table, summary_table = assemble_tables(all_results, xshooter, xsh_results, failures, xsh_sigmas)
    write_tables(lines_table, summary_table, validation)

    if DIAGNOSTIC_CUTOUTS:
        first = all_results[0]
        cutoutPath = OUTDIR / f"{OUTFILE_PREFIX}_cutouts.pdf"
        write_diagnostic_cutouts(first["wave"], first["flux"], first["results_original"], cutoutPath)

    plot_offset_grid(all_results)
    plot_summary_offset_vs_wavelength(all_results)

    print("Done.")


if __name__ == "__main__":
    main()
