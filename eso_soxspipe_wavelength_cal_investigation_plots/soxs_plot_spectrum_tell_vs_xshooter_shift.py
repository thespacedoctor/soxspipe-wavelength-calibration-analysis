#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# =========================
# INPUTS
# =========================
SOXS_ROOT = Path("2026-01")
XSHOOTER_ROOT = Path("xshooter_extractions")
TELL_FILE = Path("telluric_5_lines.dat")

SOXS_MATCH_STRING = "SOXS_LTT3218"
XSHOOTER_MATCH_STRING = "Hip012377_TELL_SLT"
XSHOOTER_FILENAME = "SDP_TELL_SLIT_FLUX_MERGE1D_NIR.fits"

MAX_XSHOOTER_FILES = 20
SOXS_FLUX_COL = "FLUX_COUNTS"

NORMALIZE = True
NORM_MODE = "median"
NORM_WAVE_RANGE = (1671, 1676)

SHIFT_REFERENCE_WINDOW = (1750.0, 1752.0)
SHIFT_POLY_DEGREE = 12

SEARCH_HALF_WINDOW_NM = 0.25
MAX_OFFSET_TO_ACCEPT_NM = SEARCH_HALF_WINDOW_NM
OFFSET_CLIP_NM = 0.2

MIN_POINTS_IN_WINDOW = 5
PARABOLA_REFINE = True
MIN_DEPTH = None

N_ROWS = 4
N_COLS = 4
PANELS_PER_FIG = N_ROWS * N_COLS

XLIM_NM = None
YLIM_OFFSET_NM = (-0.3, 0.3)

SAVE_FIG = True
SHOW_FIG = True

OUTDIR = SOXS_ROOT / "plots_soxs_vs_xshooter"
OUTFILE_PREFIX = "LTT3218_telluric_offsets_grid"

COLOR_SOXS = "tab:cyan"
COLOR_SOXS_SHIFTED = "tab:blue"
COLOR_XSHOOTER = "tab:orange"


# =========================
# BASIC UTILITIES
# =========================
def vac_to_air_um(wave_um):
    wave_um = np.asarray(wave_um, dtype=float)

    a = 8.34213e-5
    b1 = 2.406030e-2
    c1 = 130.0
    b2 = 1.5997e-4
    c2 = 38.9

    n = 1.0 + a + wave_um**2 * b1 / (wave_um**2 * c1 - 1.0) + wave_um**2 * b2 / (wave_um**2 * c2 - 1.0)

    return wave_um / n


def load_hitran_air_nm(filepath: Path):
    if not filepath.exists():
        raise FileNotFoundError(f"Missing HITRAN/telluric file: {filepath}")

    lines_nm_vac = np.loadtxt(filepath, dtype=float)
    lines_nm_vac = np.atleast_1d(lines_nm_vac)
    lines_nm_vac = lines_nm_vac[np.isfinite(lines_nm_vac)]

    lines_um_air = vac_to_air_um(lines_nm_vac / 1000.0)
    lines_nm_air = lines_um_air * 1000.0

    return np.unique(lines_nm_air)


def get_header_value(header, keys, default="UNKNOWN"):
    for key in keys:
        if key in header:
            val = str(header[key]).strip()
            if val:
                return val
    return default


def get_date_obs(header):
    return get_header_value(
        header,
        ["DATE-OBS", "DATE_OBS", "HIERARCH ESO OBS START", "MJD-OBS"],
        default="UNKNOWN",
    )


def clean_spectrum(wave, flux):
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


def is_in_corron_path(path: Path):
    return any(part.endswith("_CorrOn") for part in path.parts)


def find_column(colnames, candidates):
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
    for idx, hdu in enumerate(hdul):
        if isinstance(hdu, fits.BinTableHDU):
            return idx
    raise ValueError("No BinTableHDU found")


def chunks(seq, size):
    for i in range(0, len(seq), size):
        yield i // size + 1, seq[i : i + size]


def clipped_good_mask(results):
    return results["ok"] & np.isfinite(results["offset_nm"]) & (np.abs(results["offset_nm"]) <= OFFSET_CLIP_NM)


# =========================
# READ SOXS
# =========================
def read_soxs_spectrum(filepath: Path):
    with fits.open(filepath) as hdul:
        header = hdul[0].header
        data = hdul[1].data
        colnames = data.columns.names

        obs_name = get_header_value(
            header,
            ["HIERARCH ESO OBS NAME", "OBS.NAME", "OBSNAME", "OBJECT"],
        )

        if SOXS_MATCH_STRING.upper() not in obs_name.upper():
            return None

        if "WAVE" not in colnames:
            raise KeyError(f"No WAVE column in {filepath}")

        if SOXS_FLUX_COL not in colnames:
            raise KeyError(f"No {SOXS_FLUX_COL} column in {filepath}. Columns: {colnames}")

        wave_nm = np.asarray(data["WAVE"], dtype=float)
        flux = np.asarray(data[SOXS_FLUX_COL], dtype=float)
        date_obs = get_date_obs(header)

    wave_nm, flux = clean_spectrum(wave_nm, flux)
    flux = normalize_flux(wave_nm, flux, filepath)

    return {
        "instrument": "SOXS",
        "wave": wave_nm,
        "flux": flux,
        "date": date_obs,
        "file": filepath,
    }


def find_soxs_spectra():
    files = sorted(SOXS_ROOT.glob("2026-*/soxs-*std/*_EXTRACTED_MERGED.fits"))
    spectra = []

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
# READ XSHOOTER
# =========================
def read_xshooter_spectrum(filepath: Path):
    with fits.open(filepath) as hdul:
        obs_name = "UNKNOWN"
        date_obs = "UNKNOWN"

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


def find_xshooter_spectra():
    files = sorted(XSHOOTER_ROOT.glob(f"**/{XSHOOTER_FILENAME}"))
    spectra = []

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
# MINIMUM MEASUREMENT
# =========================
def fit_local_parabola_minimum(x, y):
    good = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x, dtype=float)[good]
    y = np.asarray(y, dtype=float)[good]

    if x.size < 3:
        return np.nan, np.nan

    try:
        a, b, c0 = np.polyfit(x, y, 2)

        if not np.isfinite(a) or a <= 0:
            return np.nan, np.nan

        x0 = -b / (2.0 * a)
        y0 = a * x0**2 + b * x0 + c0

        return x0, y0

    except Exception:
        return np.nan, np.nan


def find_local_minimum_near_line(
    wave_nm,
    flux,
    line_nm,
    half_window_nm=0.25,
    min_points=5,
    parabola_refine=True,
    min_depth=None,
    max_offset_nm=0.25,
):
    m = (
        np.isfinite(wave_nm)
        & np.isfinite(flux)
        & (wave_nm >= line_nm - half_window_nm)
        & (wave_nm <= line_nm + half_window_nm)
    )

    w = wave_nm[m]
    f = flux[m]

    if w.size < min_points:
        return False, np.nan, np.nan, "too_few_points"

    i0 = np.nanargmin(f)
    wmin = w[i0]
    fmin = f[i0]

    if min_depth is not None:
        depth = np.nanmedian(f) - fmin
        if not np.isfinite(depth) or depth < min_depth:
            return False, np.nan, np.nan, "too_shallow"

    if parabola_refine and 1 <= i0 <= w.size - 2:
        xloc = w[i0 - 1 : i0 + 2]
        yloc = f[i0 - 1 : i0 + 2]
        xfit, yfit = fit_local_parabola_minimum(xloc, yloc)

        if np.isfinite(xfit) and xloc.min() <= xfit <= xloc.max():
            wmin = xfit
            fmin = yfit

    offset_nm = wmin - line_nm

    if not np.isfinite(offset_nm) or abs(offset_nm) >= max_offset_nm:
        return False, np.nan, np.nan, "offset_too_large"

    return True, wmin, offset_nm, "ok"


def measure_reference_minimum_for_shift(spec, hitran_ref_nm):
    wave = spec["wave"]
    flux = spec["flux"]

    w0, w1 = SHIFT_REFERENCE_WINDOW
    mask = (wave >= w0) & (wave <= w1) & np.isfinite(flux)

    if np.sum(mask) < SHIFT_POLY_DEGREE + 2:
        raise ValueError(f"Not enough points in {SHIFT_REFERENCE_WINDOW} nm for {spec['file']}")

    x = wave[mask]
    y = flux[mask]

    x_ref = np.nanmean(x)
    coeff = np.polyfit(x - x_ref, y, SHIFT_POLY_DEGREE)
    poly = np.poly1d(coeff)

    fine_wave = np.linspace(w0, w1, 3000)
    fine_flux = poly(fine_wave - x_ref)

    min_wave = fine_wave[np.nanargmin(fine_flux)]
    shift_nm = hitran_ref_nm - min_wave

    return min_wave, shift_nm


def nearest_hitran_reference_line(tell_nm_air):
    w0, w1 = SHIFT_REFERENCE_WINDOW
    center = 0.5 * (w0 + w1)

    lines = tell_nm_air[(tell_nm_air >= w0) & (tell_nm_air <= w1)]

    if lines.size == 0:
        raise RuntimeError(f"No HITRAN air line found in {SHIFT_REFERENCE_WINDOW}")

    return lines[np.argmin(np.abs(lines - center))]


def apply_soxs_shift(soxs, hitran_ref_nm):
    shifted = []

    for spec in soxs:
        try:
            min_wave, shift_nm = measure_reference_minimum_for_shift(spec, hitran_ref_nm)

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
def measure_offsets_for_spectrum(wave_nm, flux, lines_nm):
    out = {
        "line_nm": [],
        "wmin_nm": [],
        "offset_nm": [],
        "ok": [],
        "reason": [],
    }

    if wave_nm.size == 0:
        lines_here = np.array([])
    else:
        lines_here = lines_nm[(lines_nm >= np.nanmin(wave_nm)) & (lines_nm <= np.nanmax(wave_nm))]

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


def build_analysis_results(soxs_shifted, xshooter, tell_nm_air):
    results = []

    for s in soxs_shifted:
        r_orig = measure_offsets_for_spectrum(s["wave"], s["flux"], tell_nm_air)
        r_shift = measure_offsets_for_spectrum(s["wave_shifted"], s["flux"], tell_nm_air)

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

    for s in xshooter:
        r = measure_offsets_for_spectrum(s["wave"], s["flux"], tell_nm_air)

        results.append(
            {
                "instrument": "XSHOOTER",
                "name": s["file"].name,
                "date": s["date"],
                "file": s["file"],
                "shift_nm": 0.0,
                "reference_min_nm": np.nan,
                "reference_hitran_nm": np.nan,
                "results_xshooter": r,
            }
        )

    return results


# =========================
# PLOTS
# =========================
def apply_plot_style():
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
    fig, axes = plt.subplots(
        N_ROWS,
        N_COLS,
        figsize=(4.5 * N_COLS, 3.0 * N_ROWS),
        squeeze=False,
    )
    return fig, axes.ravel()


def panel_title(res):
    if res["instrument"] == "SOXS":
        return f"SOXS: {res['date']} | shift={res['shift_nm']:+.4f} nm\n" f"{res['name']}"

    return f"XSHOOTER: {res['date']}\n{res['name']}"


def add_offset_stats_box(ax, label, offsets):
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
    apply_plot_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for page, page_results in chunks(all_results, PANELS_PER_FIG):
        fig, axes = make_axes_grid()

        for i, (ax, res) in enumerate(zip(axes, page_results)):
            row = i // N_COLS
            col = i % N_COLS

            is_left = col == 0
            is_bottom = row == N_ROWS - 1

            if res["instrument"] == "SOXS":
                r_orig = res["results_original"]
                r_shift = res["results_shifted"]

                good_orig = clipped_good_mask(r_orig)
                good_shift = clipped_good_mask(r_shift)

                # if np.any(good_orig):
                #     ax.plot(
                #         r_orig["line_nm"][good_orig],
                #         r_orig["offset_nm"][good_orig],
                #         marker="x",
                #         linestyle="-",
                #         lw=0.8,
                #         ms=3,
                #         color=COLOR_SOXS,
                #         alpha=0.75,
                #         label="SOXS original - HITRAN",
                #     )

                if np.any(good_shift):
                    ax.plot(
                        r_shift["line_nm"][good_shift],
                        r_shift["offset_nm"][good_shift],
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
                    r_shift["offset_nm"][good_shift],
                )

            else:
                r = res["results_xshooter"]
                good = clipped_good_mask(r)

                if np.any(good):
                    ax.plot(
                        r["line_nm"][good],
                        r["offset_nm"][good],
                        marker="o",
                        linestyle="-",
                        lw=0.8,
                        ms=3,
                        color=COLOR_XSHOOTER,
                        alpha=0.85,
                        label="XSHOOTER - HITRAN",
                    )

                add_offset_stats_box(
                    ax,
                    "XSHOOTER",
                    r["offset_nm"][good],
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

        for ax in axes[len(page_results) :]:
            ax.axis("off")

        axes[0].legend(fontsize=7, loc="best")

        fig.suptitle(
            "Telluric-line wavelength offsets w.r.t. HITRAN",
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
    fig, ax = plt.subplots(figsize=(13, 6))

    for res in all_results:
        if res["instrument"] == "SOXS":
            r_orig = res["results_original"]
            r_shift = res["results_shifted"]

            good_orig = clipped_good_mask(r_orig)
            good_shift = clipped_good_mask(r_shift)

            ax.plot(
                r_orig["line_nm"][good_orig],
                r_orig["offset_nm"][good_orig],
                marker="x",
                linestyle="none",
                ms=3,
                color=COLOR_SOXS,
                alpha=0.35,
            )

            ax.plot(
                r_shift["line_nm"][good_shift],
                r_shift["offset_nm"][good_shift],
                marker="o",
                linestyle="none",
                ms=3,
                color=COLOR_SOXS_SHIFTED,
                alpha=0.45,
            )

        else:
            r = res["results_xshooter"]
            good = clipped_good_mask(r)

            ax.plot(
                r["line_nm"][good],
                r["offset_nm"][good],
                marker="o",
                linestyle="none",
                ms=3,
                color=COLOR_XSHOOTER,
                alpha=0.35,
            )

    ax.axhline(0.0, color="black", ls="--", lw=0.8)

    handles = [
        plt.Line2D([], [], color=COLOR_SOXS, marker="x", ls="none", label="SOXS original"),
        plt.Line2D([], [], color=COLOR_SOXS_SHIFTED, marker="o", ls="none", label="SOXS shifted"),
        plt.Line2D([], [], color=COLOR_XSHOOTER, marker="o", ls="none", label="XSHOOTER"),
    ]

    ax.legend(handles=handles, frameon=False)
    ax.set_xlabel("HITRAN air wavelength [nm]")
    ax.set_ylabel("Measured minimum − HITRAN [nm]")
    ax.set_title("All telluric-line offsets")
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


# =========================
# MAIN
# =========================
def main():
    tell_nm_air = load_hitran_air_nm(TELL_FILE)
    print(f"Loaded {tell_nm_air.size} HITRAN air lines.")

    hitran_ref_nm = nearest_hitran_reference_line(tell_nm_air)
    print(f"Reference HITRAN line for SOXS shift: {hitran_ref_nm:.6f} nm")

    print("Searching SOXS spectra...")
    soxs = find_soxs_spectra()

    print("Searching XSHOOTER spectra...")
    xshooter = find_xshooter_spectra()

    print(f"Found {len(soxs)} SOXS spectra.")
    print(f"Found {len(xshooter)} XSHOOTER spectra.")

    print("Applying SOXS shifts...")
    soxs_shifted = apply_soxs_shift(soxs, hitran_ref_nm)

    print("\n=== SOXS SHIFTS ===")
    for s in soxs_shifted:
        print(
            f"{s['date']}  "
            f"min={s['shift_reference_min_nm']:.6f} nm  "
            f"HITRAN={s['shift_reference_hitran_nm']:.6f} nm  "
            f"shift={s['shift_nm']:+.6f} nm  "
            f"{s['file'].name}"
        )

    print("\nMeasuring offsets for all HITRAN lines...")
    all_results = build_analysis_results(soxs_shifted, xshooter, tell_nm_air)

    if not all_results:
        print("No results to plot.")
        return

    plot_offset_grid(all_results)
    plot_summary_offset_vs_wavelength(all_results)

    print("Done.")


if __name__ == "__main__":
    main()
