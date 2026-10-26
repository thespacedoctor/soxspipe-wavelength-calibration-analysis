import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import soxs_telluric_wavelength_offsets as m

STEP_NM = 0.06
LINE_SIGMA_NM = 0.06
LINE_SPACING_NM = 3.0
WAVE_RANGE_NM = (1700.0, 1760.0)


def make_wave(wmin=WAVE_RANGE_NM[0], wmax=WAVE_RANGE_NM[1]):
    return np.arange(wmin, wmax, STEP_NM)


def make_lines(wmin=WAVE_RANGE_NM[0], wmax=WAVE_RANGE_NM[1]):
    return np.arange(wmin + 1.5, wmax - 1.5, LINE_SPACING_NM) + 0.013


def make_spectrum(lines, offset=0.0, depth=0.5, noise=0.003, seed=1, wave=None, sigma=LINE_SIGMA_NM):
    # CONTINUUM OF 1 WITH GAUSSIAN ABSORPTION LINES AT LINE + OFFSET.
    wave = make_wave() if wave is None else wave
    rng = np.random.default_rng(seed)
    flux = np.ones_like(wave)
    depths = np.broadcast_to(depth, np.shape(lines))
    for line, d in zip(lines, depths):
        flux -= d * np.exp(-0.5 * ((wave - (line + offset)) / sigma) ** 2)
    return wave, flux + rng.normal(0.0, noise, wave.size)


def make_sloped_line_spectrum(line, slope, fractionalDepth, level=1.0, sigma=0.03):
    # SLOPED CONTINUUM WITH ONE NARROW LINE ON A GRID SAMPLE (NOISE-FREE).
    wave = make_wave()
    continuum = level + slope * (wave - line)
    return wave, continuum * (1.0 - fractionalDepth * np.exp(-0.5 * ((wave - line) / sigma) ** 2))


def make_offsets_frame(lines, offsets, ok=None, clipped=None):
    ok = np.ones(len(lines), dtype=bool) if ok is None else ok
    clipped = np.zeros(len(lines), dtype=bool) if clipped is None else clipped
    return pd.DataFrame(
        {
            "line_nm": lines,
            "wmin_nm": lines + offsets,
            "offset_nm": offsets,
            "depth": 0.5,
            "core_flux": 0.5,
            "continuum": 1.0,
            "core_ratio": 0.5,
            "continuum_mode": "two_sided",
            "ok": ok,
            "reason": np.where(ok, "ok", "too_shallow"),
            "parabola_fallback": False,
            "clipped": clipped,
        }
    )


# =========================
# find_local_minimum_near_line
# =========================
def test_recovers_injected_offset_within_tolerance():
    line = 1730.013
    wave, flux = make_spectrum([line], offset=0.037)

    result = m.find_local_minimum_near_line(wave, flux, line)

    assert result["ok"]
    assert result["reason"] == "ok"
    assert result["offset_nm"] == pytest.approx(0.037, abs=0.01)


def test_returns_too_few_points_for_empty_window_without_crashing():
    wave, flux = make_spectrum([1730.013])

    result = m.find_local_minimum_near_line(wave, flux, 1900.0)

    assert not result["ok"]
    assert result["reason"] == "too_few_points"
    assert np.isnan(result["offset_nm"])


def test_keeps_too_few_points_reason_when_parabola_refine_enabled():
    wave = np.array([1730.00, 1730.06, 1730.12])
    flux = np.array([1.0, 0.5, 1.0])

    result = m.find_local_minimum_near_line(wave, flux, 1730.06, parabola_refine=True)

    assert result["reason"] == "too_few_points"


def test_flags_edge_minimum_when_lowest_sample_is_at_window_edge():
    line = 1730.013
    wave = make_wave()
    flux = 1.0 - 0.5 * (wave - line) / 10.0 + 0.0 * wave

    result = m.find_local_minimum_near_line(wave, flux, line)

    assert not result["ok"]
    assert result["reason"] == "edge_minimum"


def test_flags_too_shallow_line():
    line = 1730.013
    wave, flux = make_spectrum([line], depth=0.02, noise=0.0)

    result = m.find_local_minimum_near_line(wave, flux, line, min_depth=0.05)

    assert not result["ok"]
    assert result["reason"] == "too_shallow"
    assert result["depth"] == pytest.approx(0.02, abs=0.005)


def test_flags_saturated_core():
    line = 1730.013
    wave, flux = make_spectrum([line], depth=0.98, noise=0.0)

    result = m.find_local_minimum_near_line(wave, flux, line, max_core_ratio=0.1)

    assert not result["ok"]
    assert result["reason"] == "saturated"
    assert result["core_ratio"] < 0.1


def test_accepts_discrete_minimum_when_parabola_refine_disabled():
    line = 1730.013
    wave, flux = make_spectrum([line], offset=0.02)

    result = m.find_local_minimum_near_line(wave, flux, line, parabola_refine=False)

    assert result["ok"]
    assert result["wmin_nm"] in wave


def test_reports_fractional_depth_core_flux_and_continuum_for_accepted_line():
    line = 1730.0
    wave, flux = make_sloped_line_spectrum(line, slope=0.0, fractionalDepth=0.5)

    result = m.find_local_minimum_near_line(wave, flux, line)

    assert result["depth"] == pytest.approx(0.5, abs=0.02)
    assert result["core_flux"] == pytest.approx(0.5, abs=0.02)
    assert result["continuum"] == pytest.approx(1.0, abs=0.02)
    assert result["core_ratio"] == pytest.approx(0.5, abs=0.02)


def test_accepts_line_on_steep_sloped_continuum_with_correct_fractional_depth():
    line = 1730.0
    wave, flux = make_sloped_line_spectrum(line, slope=0.4, fractionalDepth=0.4)

    result = m.find_local_minimum_near_line(wave, flux, line)

    assert result["ok"]
    assert result["depth"] == pytest.approx(0.4, abs=0.03)


def test_rejects_shallow_line_on_slope_that_a_window_percentile_would_accept():
    line = 1730.0
    wave, flux = make_sloped_line_spectrum(line, slope=0.1, fractionalDepth=0.04)

    result = m.find_local_minimum_near_line(wave, flux, line, min_depth=0.05)

    assert result["reason"] == "too_shallow"
    assert result["depth"] == pytest.approx(0.04, abs=0.01)


@pytest.mark.parametrize("scale", [0.5, 1.0, 2.0])
def test_saturation_classification_is_invariant_to_flux_scale(scale):
    line = 1730.0
    wave, deepFlux = make_sloped_line_spectrum(line, slope=0.0, fractionalDepth=0.95)
    _, midFlux = make_sloped_line_spectrum(line, slope=0.0, fractionalDepth=0.5)

    deep = m.find_local_minimum_near_line(wave, scale * deepFlux, line)
    mid = m.find_local_minimum_near_line(wave, scale * midFlux, line)

    assert deep["reason"] == "saturated"
    assert mid["reason"] == "ok"
    assert deep["core_ratio"] == pytest.approx(0.05, abs=0.01)
    assert deep["depth"] == pytest.approx(0.95, abs=0.01)


def test_blended_neighbour_in_sideband_does_not_wreck_the_continuum():
    line = 1730.0
    wave, flux = make_sloped_line_spectrum(line, slope=0.0, fractionalDepth=0.7)
    # DEEP NEIGHBOUR 0.18 NM AWAY: OUTSIDE THE 0.1 NM CORE (2.5 * 0.04) BUT INSIDE THE SIDEBAND.
    flux = flux - 0.4 * np.exp(-0.5 * ((wave - (line + 0.18)) / 0.02) ** 2)

    result = m.find_local_minimum_near_line(wave, flux, line, line_sigma_nm=0.04)

    assert result["ok"]
    assert result["continuum"] == pytest.approx(1.0, abs=0.02)
    assert result["depth"] == pytest.approx(0.7, abs=0.03)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_offset_line_wing_does_not_leak_into_the_sidebands(seed):
    # THE LINE SITS 0.08 NM FROM THE HITRAN POSITION; THE CORE MUST FOLLOW THE MEASURED MINIMUM.
    line = 1730.0
    centre = line + 0.08
    wave = make_wave()
    wave = wave + (centre - wave[np.argmin(np.abs(wave - centre))])
    wave, flux = make_spectrum([centre], depth=0.5, noise=0.003, seed=seed, wave=wave, sigma=0.045)

    result = m.find_local_minimum_near_line(wave, flux, line)

    assert result["ok"]
    assert result["offset_nm"] == pytest.approx(0.08, abs=0.01)
    assert result["depth"] == pytest.approx(0.5, abs=0.02)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_continuum_window_follows_a_minimum_far_from_the_hitran_position(seed):
    # THE LINE SITS -0.15 NM FROM THE HITRAN POSITION: THE SIDEBAND WINDOW MUST FOLLOW THE MINIMUM.
    line = 1730.0
    centre = line - 0.15
    wave = make_wave()
    wave = wave + (centre - wave[np.argmin(np.abs(wave - centre))])
    wave, flux = make_spectrum([centre], depth=0.5, noise=0.003, seed=seed, wave=wave, sigma=0.045)

    result = m.find_local_minimum_near_line(wave, flux, line)

    assert result["ok"]
    assert result["offset_nm"] == pytest.approx(-0.15, abs=0.01)
    assert result["continuum"] == pytest.approx(1.0, abs=0.02)
    assert result["depth"] == pytest.approx(0.5, abs=0.02)


def test_reports_no_continuum_when_one_sideband_is_missing():
    line = 1730.0
    wave, flux = make_sloped_line_spectrum(line, slope=0.0, fractionalDepth=0.5)
    keep = wave >= line - 0.1

    result = m.find_local_minimum_near_line(wave[keep], flux[keep], line, sideband_width_nm=0.1)

    assert result["reason"] == "no_continuum"
    assert np.isnan(result["continuum"])
    assert result["continuum_mode"] == "none"
    assert not result["ok"]


def test_reports_no_continuum_when_too_few_sideband_points():
    line = 1730.0
    wave, flux = make_sloped_line_spectrum(line, slope=0.0, fractionalDepth=0.5)

    result = m.find_local_minimum_near_line(wave, flux, line, sideband_width_nm=0.1)

    assert result["reason"] == "no_continuum"


def test_reports_no_continuum_when_fitted_continuum_is_not_positive():
    line = 1730.0
    wave, _ = make_sloped_line_spectrum(line, slope=0.0, fractionalDepth=0.5)
    flux = -1.0 - 0.5 * np.exp(-0.5 * ((wave - line) / 0.03) ** 2)

    result = m.find_local_minimum_near_line(wave, flux, line)

    assert result["reason"] == "no_continuum"


def test_sidebands_are_independent_of_the_minimum_search_window():
    line = 1730.0
    wave, flux = make_sloped_line_spectrum(line, slope=0.0, fractionalDepth=0.5)

    result = m.find_local_minimum_near_line(wave, flux, line, half_window_nm=0.15)

    assert result["ok"]
    assert result["continuum"] == pytest.approx(1.0, abs=0.02)


def test_continuum_at_minimum_evaluates_the_sideband_line_at_the_minimum_and_ignores_the_core():
    wave = np.arange(-0.25, 0.2501, 0.06)
    flux = 1.0 + 0.2 * wave
    flux[np.abs(wave) <= 0.1] = 0.1

    continuum, mode = m.continuum_at_minimum(wave + 1730.0, flux, 1730.0, core_half_width_nm=0.1)

    assert continuum == pytest.approx(1.0)
    assert mode == "two_sided"


def test_continuum_at_minimum_is_unavailable_without_enough_sideband_points():
    wave = np.arange(0.0, 0.1501, 0.06)

    continuum, mode = m.continuum_at_minimum(wave + 1730.0, np.ones_like(wave), 1730.0, core_half_width_nm=0.1)

    assert np.isnan(continuum)
    assert mode == "none"


def make_gaussian_dip(wave, centre, depth, sigma, level=1.0):
    return level * (1.0 - depth * np.exp(-0.5 * ((wave - centre) / sigma) ** 2))


@pytest.mark.parametrize("sigma", [0.12, 0.14])
def test_broad_saturated_line_gets_unbiased_continuum_and_is_saturated(sigma):
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.95, sigma)

    result = m.find_local_minimum_near_line(wave, flux, line, line_sigma_nm=sigma)

    assert result["continuum"] == pytest.approx(1.0, abs=0.03)
    assert result["reason"] == "saturated"


def test_broad_shallow_line_keeps_correct_fractional_depth():
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.3, 0.14)

    result = m.find_local_minimum_near_line(wave, flux, line, line_sigma_nm=0.14)

    assert result["ok"]
    assert result["depth"] == pytest.approx(0.3, abs=0.03)
    assert result["continuum"] == pytest.approx(1.0, abs=0.03)


def test_core_width_scales_with_the_line_sigma():
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.7, 0.14)

    narrowCore = m.find_local_minimum_near_line(wave, flux, line, line_sigma_nm=0.04)
    broadCore = m.find_local_minimum_near_line(wave, flux, line, line_sigma_nm=0.14)

    assert abs(broadCore["continuum"] - 1.0) < abs(narrowCore["continuum"] - 1.0)


def test_neighbour_hitran_line_is_masked_out_of_the_sidebands():
    line = 1730.0
    wave = make_wave()
    # THE MASK IS 1 SIGMA (0.04 NM) WIDE, SO THE NEIGHBOUR IS NARROW ENOUGH FOR THE MASK TO COVER ITS CORE.
    flux = make_gaussian_dip(wave, line, 0.7, 0.03) - 0.6 * np.exp(-0.5 * ((wave - (line + 0.19)) / 0.02) ** 2)

    masked = m.find_local_minimum_near_line(
        wave, flux, line, line_sigma_nm=0.04, other_lines_nm=np.array([line + 0.19])
    )
    unmasked = m.find_local_minimum_near_line(wave, flux, line, line_sigma_nm=0.04)

    assert masked["continuum"] == pytest.approx(1.0, abs=0.02)
    assert abs(unmasked["continuum"] - 1.0) > abs(masked["continuum"] - 1.0)


def test_neighbour_line_one_and_a_half_sigma_beyond_the_core_does_not_empty_a_sideband():
    # LINE SIGMA 0.1: CORE 0.25 NM, SIDEBAND 0.25-0.5 NM. THE NEIGHBOUR AT +0.55 NM WOULD HAVE MASKED 0.30-0.80 NM
    # WITH A CORE-WIDE MASK, LEAVING THE RIGHT SIDEBAND EMPTY. THE 1 SIGMA MASK ONLY REMOVES 0.45-0.65 NM.
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.5, 0.04)

    result = m.find_local_minimum_near_line(wave, flux, line, other_lines_nm=np.array([line + 0.55]))

    assert result["ok"]
    assert result["continuum_mode"] == "two_sided"
    assert result["depth"] == pytest.approx(0.5, abs=0.02)


def test_one_sided_continuum_gives_correct_depth_on_flat_continuum():
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.5, 0.04)
    keep = wave >= line - 0.1

    result = m.find_local_minimum_near_line(wave[keep], flux[keep], line, line_sigma_nm=0.04)

    assert result["ok"]
    assert result["continuum_mode"] == "one_sided"
    assert result["continuum"] == pytest.approx(1.0, abs=0.01)
    assert result["depth"] == pytest.approx(0.5, abs=0.02)


def test_one_sided_continuum_works_on_the_left_side_too():
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.5, 0.04)
    keep = wave <= line + 0.1

    result = m.find_local_minimum_near_line(wave[keep], flux[keep], line, line_sigma_nm=0.04)

    assert result["continuum_mode"] == "one_sided"
    assert result["depth"] == pytest.approx(0.5, abs=0.02)


def test_one_sided_continuum_applies_when_the_short_side_has_a_single_point():
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.5, 0.04)
    # KEEP ONE POINT (1729.88, 0.12 NM FROM THE LINE) ON THE LEFT SIDEBAND.
    keep = wave >= line - 0.13

    result = m.find_local_minimum_near_line(wave[keep], flux[keep], line, line_sigma_nm=0.04)

    assert result["continuum_mode"] == "one_sided"


def test_two_sided_mode_is_recorded_when_both_sidebands_have_points():
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.5, 0.04)

    result = m.find_local_minimum_near_line(wave, flux, line, line_sigma_nm=0.04)

    assert result["continuum_mode"] == "two_sided"


def test_one_sided_continuum_is_a_robust_median():
    values = np.array([1.0, 1.02, 0.98, 0.4, 1.01, 0.99])

    assert m.one_sided_continuum(values) == pytest.approx(1.0, abs=0.02)
    assert np.mean(values) < 0.95


def test_neighbour_displaced_by_the_spectrum_shift_is_masked_when_the_mask_offset_is_given():
    # EVERY LINE IN THE DATA SITS 0.12 NM REDWARD OF ITS HITRAN POSITION, INCLUDING THE NEIGHBOUR AT +0.35 NM.
    line, shift = 1730.0, 0.12
    wave = make_wave()
    wave = wave + (line + shift - wave[np.argmin(np.abs(wave - (line + shift)))])
    flux = make_gaussian_dip(wave, line + shift, 0.5, 0.04) - 0.6 * np.exp(
        -0.5 * ((wave - (line + 0.35 + shift)) / 0.02) ** 2
    )
    lines, allLines = np.array([line]), np.array([line, line + 0.35])

    aligned = m.measure_offsets_for_spectrum(wave, flux, lines, mask_lines_nm=allLines, mask_offset_nm=shift)
    unaligned = m.measure_offsets_for_spectrum(wave, flux, lines, mask_lines_nm=allLines)

    assert aligned["continuum"].iloc[0] == pytest.approx(1.0, abs=0.02)
    assert abs(unaligned["continuum"].iloc[0] - 1.0) > abs(aligned["continuum"].iloc[0] - 1.0)


def test_one_sided_continuum_bias_on_a_sloped_continuum_stays_within_the_documented_bound():
    # A FLAT ONE-SIDED CONTINUUM IS BIASED BY ABOUT SLOPE * (DISTANCE OF THE SIDEBAND FROM THE MINIMUM).
    line, slope, coreHalf, width = 1730.0, 0.4, 0.1, 0.25
    wave, flux = make_sloped_line_spectrum(line, slope=slope, fractionalDepth=0.5)
    keep = wave >= line - 0.1

    result = m.find_local_minimum_near_line(wave[keep], flux[keep], line, line_sigma_nm=0.04)

    bias = result["continuum"] - 1.0
    assert result["continuum_mode"] == "one_sided"
    assert 0.0 < bias <= slope * (coreHalf + width)
    assert bias == pytest.approx(slope * (coreHalf + width / 2), abs=0.04)


def test_nan_pixel_in_the_sideband_still_gives_a_fit():
    wave = np.arange(-0.4, 0.4001, 0.06) + 1730.0
    flux = 1.0 + 0.2 * (wave - 1730.0)
    flux[np.abs(wave - 1730.0) <= 0.1] = 0.1
    flux[3] = np.nan

    continuum, mode = m.continuum_at_minimum(wave, flux, 1730.0, core_half_width_nm=0.1)

    assert continuum == pytest.approx(1.0)
    assert mode == "two_sided"


def test_infinite_pixel_in_the_sideband_is_ignored():
    wave = np.arange(-0.4, 0.4001, 0.06) + 1730.0
    flux = np.ones_like(wave)
    flux[2] = np.inf

    continuum, _ = m.continuum_at_minimum(wave, flux, 1730.0, core_half_width_nm=0.1)

    assert continuum == pytest.approx(1.0, abs=1e-9)


def test_tiny_flux_scale_gives_the_same_classification_as_unit_scale():
    line = 1730.0
    wave = make_wave()
    flux = make_gaussian_dip(wave, line, 0.6, 0.04) - 0.4 * np.exp(-0.5 * ((wave - (line + 0.18)) / 0.02) ** 2)

    unit = m.find_local_minimum_near_line(wave, flux, line, line_sigma_nm=0.04)
    tiny = m.find_local_minimum_near_line(wave, 1e-16 * flux, line, line_sigma_nm=0.04)

    assert tiny["reason"] == unit["reason"]
    assert tiny["depth"] == pytest.approx(unit["depth"], abs=1e-6)
    assert tiny["continuum"] == pytest.approx(1e-16 * unit["continuum"], rel=1e-6)


def test_few_noisy_sideband_points_give_negligible_mean_continuum_bias():
    # THREE POINTS PER SIDE; THE ONE-SIDED REJECTION MUST NOT PUSH THE CONTINUUM UP SYSTEMATICALLY.
    x = np.array([-0.34, -0.30, -0.26, 0.26, 0.30, 0.34])
    biases = []
    for seed in range(300):
        y = 1.0 + np.random.default_rng(seed).normal(0.0, 0.02, x.size)
        continuum, _ = m.continuum_at_minimum(x + 1730.0, y, 1730.0, core_half_width_nm=0.25)
        biases.append(continuum - 1.0)

    assert abs(np.mean(biases)) < 0.01


def test_fit_keeps_first_fit_when_rejection_would_leave_one_point_per_side():
    # THE DIPPED POINT ON THE SPARSE RIGHT SIDE IS REJECTED, WHICH WOULD LEAVE ONLY ONE POINT THERE.
    wave = np.array([-0.50, -0.46, -0.42, -0.38, -0.34, -0.30, -0.26, 0.26, 0.30]) + 1730.0
    flux = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6])

    continuum, _ = m.continuum_at_minimum(wave, flux, 1730.0, core_half_width_nm=0.25, sideband_width_nm=0.3)

    assert continuum < 0.95


# =========================
# measure_offsets_for_spectrum
# =========================
def test_measure_offsets_returns_one_row_per_line_in_coverage():
    lines = make_lines()
    wave, flux = make_spectrum(lines, offset=0.04)
    outside = np.array([1500.0, 2000.0])

    df = m.measure_offsets_for_spectrum(wave, flux, np.concatenate([outside, lines]))

    assert len(df) == lines.size
    expected = {"line_nm", "wmin_nm", "offset_nm", "depth", "core_flux", "continuum", "core_ratio", "ok", "reason"}
    assert expected <= set(df.columns)
    assert np.nanmedian(df.loc[df["ok"], "offset_nm"]) == pytest.approx(0.04, abs=0.01)


def test_measure_offsets_handles_empty_spectrum():
    df = m.measure_offsets_for_spectrum(np.array([]), np.array([]), make_lines())

    assert len(df) == 0
    assert "offset_nm" in df.columns


# =========================
# STATS AND CLIPPING
# =========================
def test_robust_stats_ignores_outliers_via_mad():
    values = np.concatenate([np.random.default_rng(3).normal(0.02, 0.03, 500), [5.0, -5.0]])

    stats = m.robust_stats(values)

    assert stats["n"] == 502
    assert stats["median"] == pytest.approx(0.02, abs=0.01)
    assert stats["sigma"] == pytest.approx(0.03, abs=0.01)


def test_robust_stats_returns_nan_for_empty_input():
    stats = m.robust_stats(np.array([]))

    assert stats["n"] == 0
    assert np.isnan(stats["median"])


def test_flag_outliers_keeps_smooth_trend_and_flags_isolated_deviation():
    lines = np.linspace(1700, 1800, 100)
    offsets = np.linspace(-0.05, 0.05, 100) + np.random.default_rng(2).normal(0, 0.005, 100)
    offsets[40] += 0.2
    df = make_offsets_frame(lines, offsets)

    outliers = m.flag_outliers(df)

    assert outliers[40]
    assert outliers.sum() <= 3


def test_flag_outliers_does_not_drop_edge_lines_or_flag_rejected_rows():
    lines = np.linspace(1700, 1800, 50)
    offsets = np.random.default_rng(4).normal(0.0, 0.01, 50)
    ok = np.ones(50, dtype=bool)
    ok[10] = False
    df = make_offsets_frame(lines, offsets, ok=ok)

    outliers = m.flag_outliers(df)

    assert not outliers[0] and not outliers[-1]
    assert not outliers[10]


def test_flag_outliers_falls_back_to_global_clipping_for_few_lines():
    lines = np.array([1700.0, 1750.0, 1800.0, 1850.0, 1900.0])
    df = make_offsets_frame(lines, np.array([0.0, 0.01, -0.01, 0.0, 0.5]))

    outliers = m.flag_outliers(df)

    assert outliers.tolist() == [False, False, False, False, True]


def test_good_offsets_excludes_rejected_and_clipped_rows():
    df = make_offsets_frame(np.array([1.0, 2.0, 3.0]), np.array([0.1, 0.2, 0.3]), ok=np.array([True, True, False]))
    df["clipped"] = [False, True, False]

    assert m.good_offsets(df).tolist() == [0.1]


# =========================
# GLOBAL SHIFT
# =========================
def test_global_shift_recovers_injected_uniform_shift():
    lines = make_lines()
    wave, flux = make_spectrum(lines, offset=0.08)

    result = m.derive_global_shift(wave, flux, lines)

    assert result["shift_nm"] == pytest.approx(-0.08, abs=0.01)
    assert result["n_lines"] >= 0.8 * lines.size
    assert result["robust_sigma_nm"] < 0.02


def test_global_shift_unaffected_by_flux_spike_at_one_line():
    lines = make_lines()
    wave, flux = make_spectrum(lines, offset=0.08)
    clean = m.derive_global_shift(wave, flux, lines)["shift_nm"]
    spikeLine = lines[len(lines) // 2]
    spike = np.abs(wave - spikeLine) < 0.2
    flux = np.where(spike, 1.3, flux)

    spiked = m.derive_global_shift(wave, flux, lines)["shift_nm"]

    assert spiked == pytest.approx(clean, abs=0.005)
    assert spiked == pytest.approx(-0.08, abs=0.01)


def test_global_shift_makes_shifted_offsets_centred_on_zero():
    lines = make_lines()
    wave, flux = make_spectrum(lines, offset=-0.07)

    shift = m.derive_global_shift(wave, flux, lines)["shift_nm"]
    df = m.measure_offsets_for_spectrum(wave + shift, flux, lines)

    assert np.nanmedian(df.loc[df["ok"], "offset_nm"]) == pytest.approx(0.0, abs=0.005)


def test_global_shift_raises_when_no_line_is_measurable():
    wave, flux = make_spectrum([], offset=0.0)

    with pytest.raises(ValueError):
        m.derive_global_shift(wave, flux, make_lines())


# =========================
# REFERENCE VALIDATION
# =========================
def test_validate_reference_lines_drops_line_missing_in_most_spectra():
    lines = np.array([1710.0, 1720.0, 1730.0])
    good = make_offsets_frame(lines, np.array([0.0, 0.01, 0.0]))
    bad = make_offsets_frame(lines, np.array([0.0, 0.01, 0.0]), ok=np.array([True, True, False]))
    bad2 = make_offsets_frame(lines, np.array([0.0, 0.01, 0.0]), ok=np.array([True, True, False]))

    table = m.validate_reference_lines([good, bad, bad2], lines, min_fraction=0.5, max_offset_nm=0.05)

    assert table["is_valid"].tolist() == [True, True, False]
    assert table.loc[table["line_nm"] == 1720.0, "median_offset_nm"].iloc[0] == pytest.approx(0.01)


def test_validate_reference_lines_rejects_line_with_large_xshooter_offset():
    lines = np.array([1710.0, 1720.0])
    frames = [make_offsets_frame(lines, np.array([0.0, 0.2])) for _ in range(3)]

    table = m.validate_reference_lines(frames, lines, min_fraction=0.5, max_offset_nm=0.05)

    assert table["is_valid"].tolist() == [True, False]


def test_validate_reference_lines_is_unavailable_with_too_few_spectra():
    lines = np.array([1710.0, 1720.0])
    frames = [make_offsets_frame(lines, np.zeros(2)) for _ in range(2)]

    table = m.validate_reference_lines(frames, lines, min_fraction=0.5, max_offset_nm=0.05, min_spectra=3)

    assert not table["validation_available"].any()
    assert not table["is_valid"].any()
    assert table["median_offset_nm"].isna().all()


def test_validate_reference_lines_is_available_at_minimum_spectrum_count():
    lines = np.array([1710.0, 1720.0])
    frames = [make_offsets_frame(lines, np.zeros(2)) for _ in range(3)]

    table = m.validate_reference_lines(frames, lines, min_fraction=0.5, max_offset_nm=0.05, min_spectra=3)

    assert table["validation_available"].all()
    assert table["is_valid"].all()


def test_validate_reference_lines_ignores_clipped_measurements():
    lines = np.array([1710.0, 1720.0])
    clean = make_offsets_frame(lines, np.array([0.0, 0.01]))
    clipped = make_offsets_frame(lines, np.array([0.0, 0.04]), clipped=np.array([False, True]))

    table = m.validate_reference_lines([clean, clipped, clipped], lines, min_fraction=0.5, max_offset_nm=0.05)

    assert table["is_valid"].tolist() == [True, False]
    assert table.loc[table["line_nm"] == 1720.0, "median_offset_nm"].iloc[0] == pytest.approx(0.01)


def test_select_measurement_lines_uses_validated_lines_when_available():
    lines = np.array([1710.0, 1720.0, 1730.0])
    table = pd.DataFrame({"line_nm": lines, "is_valid": [True, False, True], "validation_available": True})

    assert m.select_measurement_lines(table).tolist() == [1710.0, 1730.0]


def test_select_measurement_lines_falls_back_to_all_lines_when_validation_unavailable():
    lines = np.array([1710.0, 1720.0, 1730.0])
    table = pd.DataFrame({"line_nm": lines, "is_valid": False, "validation_available": False})

    assert m.select_measurement_lines(table).tolist() == lines.tolist()


def test_differential_offsets_subtract_reference_median():
    lines = np.array([1710.0, 1720.0, 1730.0])
    df = make_offsets_frame(lines, np.array([0.10, 0.05, 0.07]))
    reference = pd.DataFrame({"line_nm": lines, "median_offset_nm": [0.01, 0.02, np.nan]})

    diff = m.differential_offsets(df, reference)

    assert diff.iloc[0] == pytest.approx(0.09)
    assert diff.iloc[1] == pytest.approx(0.03)
    assert np.isnan(diff.iloc[2])


# =========================
# TABLES AND DIAGNOSTICS
# =========================
def test_reject_counts_group_rows_by_reason():
    df = make_offsets_frame(np.arange(4.0), np.zeros(4), ok=np.array([True, True, False, False]))
    df.loc[3, "reason"] = "saturated"

    counts = m.reject_counts(df)

    assert counts["ok"] == 2
    assert counts["too_shallow"] == 1
    assert counts["saturated"] == 1
    assert counts["edge_minimum"] == 0


def test_write_diagnostic_cutouts_creates_multi_page_pdf(tmp_path):
    lines = make_lines()
    wave, flux = make_spectrum(lines, offset=0.0, depth=np.linspace(0.01, 0.5, lines.size))
    df = m.add_clip_flags(m.measure_offsets_for_spectrum(wave, flux, lines))
    out = tmp_path / "cutouts.pdf"

    nPanels = m.write_diagnostic_cutouts(wave, flux, df, out, panels_per_page=4)

    assert nPanels == int((~df["ok"] | df["clipped"]).sum())
    assert nPanels > 4
    assert len(re.findall(rb"/Type /Page(?!s)", out.read_bytes())) == int(np.ceil(nPanels / 4))


def test_write_diagnostic_cutouts_skips_file_when_nothing_to_show(tmp_path):
    lines = make_lines()
    wave, flux = make_spectrum(lines, depth=0.5)
    df = m.add_clip_flags(m.measure_offsets_for_spectrum(wave, flux, lines))
    df["ok"] = True
    df["clipped"] = False
    out = tmp_path / "cutouts.pdf"

    assert m.write_diagnostic_cutouts(wave, flux, df, out) == 0
    assert not out.exists()


def test_build_summary_row_reports_shift_stats_and_reject_counts():
    lines = make_lines()
    wave, flux = make_spectrum(lines, offset=0.08)
    shiftInfo = m.derive_global_shift(wave, flux, lines)
    original = m.add_clip_flags(m.measure_offsets_for_spectrum(wave, flux, lines))
    shifted = m.add_clip_flags(m.measure_offsets_for_spectrum(wave + shiftInfo["shift_nm"], flux, lines))

    sigmaInfo = {"sigma_nm": 0.12, "n_fits": 25, "n_failed_fits": 3, "sigma_mad_nm": 0.015}
    row = m.build_summary_row("test.fits", "SOXS", original, shifted, shiftInfo, sigma_info=sigmaInfo)

    assert row["spectrum"] == "test.fits"
    assert row["shift_nm"] == pytest.approx(-0.08, abs=0.01)
    assert row["median_original_nm"] == pytest.approx(0.08, abs=0.01)
    assert row["median_shifted_nm"] == pytest.approx(0.0, abs=0.01)
    assert row["n_ok"] == int(original["ok"].sum())
    assert "reason_edge_minimum" in row
    assert row["line_sigma_nm"] == 0.12
    assert row["n_width_fits"] == 25
    assert row["n_width_failed"] == 3
    assert row["line_sigma_mad_nm"] == 0.015
    assert row["mode_two_sided"] + row["mode_one_sided"] <= row["n_measured"]
    assert row["mode_one_sided"] == int((original["continuum_mode"] == "one_sided").sum())


# =========================
# PIPELINE (INTEGRATION)
# =========================
def make_pipeline_inputs(tmp_path):
    lines = make_lines()
    # THE LAST LINE IS ABSENT FROM THE X-SHOOTER SPECTRA, SO IT MUST NOT BE VALIDATED.
    xshLines = lines[:-1]
    xsh = [
        {
            "wave": make_wave(),
            "flux": make_spectrum(xshLines, offset=0.0, seed=s)[1],
            "file": tmp_path / f"x{s}" / "x.fits",
        }
        for s in range(3)
    ]
    wave, flux = make_spectrum(lines, offset=0.08, seed=10)
    soxs = [{"wave": wave, "flux": flux, "file": tmp_path / "soxs.fits", "date": "2026-01-01"}]
    return lines, xsh, soxs


def test_pipeline_validates_lines_and_builds_shifted_results_with_differentials(tmp_path):
    lines, xsh, soxs = make_pipeline_inputs(tmp_path)
    xshResults, _xshSigmas = m.measure_spectra(xsh, lines)
    validation = m.validate_reference_lines(xshResults, lines)
    validLines = validation.loc[validation["is_valid"], "line_nm"].to_numpy()

    results, failures = m.build_analysis_results(soxs, validLines, validation)

    assert failures == []
    assert validLines.size == lines.size - 1
    assert len(results) == 1
    assert results[0]["shift_nm"] == pytest.approx(-0.08, abs=0.01)
    shiftedStats = m.robust_stats(m.good_offsets(results[0]["results_shifted"]))
    assert shiftedStats["median"] == pytest.approx(0.0, abs=0.01)


def test_pipeline_records_spectrum_without_measurable_lines_as_failed(tmp_path):
    lines, _xsh, soxs = make_pipeline_inputs(tmp_path)
    soxs[0]["flux"] = np.ones_like(soxs[0]["flux"])
    validation = m.validate_reference_lines([], lines)

    results, failures = m.build_analysis_results(soxs, lines, validation)

    assert results == []
    assert len(failures) == 1
    assert failures[0]["file"] == soxs[0]["file"]
    assert "accepted line offsets" in failures[0]["error"]


def test_shifted_offsets_use_refine_window_so_boundary_contaminants_are_ignored():
    lines = make_lines()
    wave, flux = make_spectrum(lines, offset=0.08, depth=0.5, seed=5)
    # A DEEPER CONTAMINANT 0.17 NM FROM 40% OF THE SHIFTED LINES: INSIDE +/-0.25 NM BUT OUTSIDE +/-0.15 NM.
    for line in np.concatenate([lines[::5], lines[1::5]]):
        flux -= 0.6 * np.exp(-0.5 * ((wave - (line + 0.08 + 0.17)) / LINE_SIGMA_NM) ** 2)
    soxs = [{"wave": wave, "flux": flux, "file": Path("soxs.fits"), "date": "2026-01-01"}]
    validation = m.validate_reference_lines([], lines)

    results, _ = m.build_analysis_results(soxs, lines, validation)

    shifted = results[0]["results_shifted"]
    shiftedOffsets = m.good_offsets(shifted)
    assert np.median(shiftedOffsets) == pytest.approx(0.0, abs=0.01)
    assert np.abs(shiftedOffsets).max() < m.REFINE_HALF_WINDOW_NM


def test_differential_is_nan_when_validation_unavailable(tmp_path):
    lines, _xsh, soxs = make_pipeline_inputs(tmp_path)
    validation = m.validate_reference_lines([], lines)

    results, _ = m.build_analysis_results(soxs, lines, validation)

    assert results[0]["differential_shifted"].isna().all()


def test_pipeline_writes_lines_summary_and_reference_csvs(tmp_path, monkeypatch):
    lines, xsh, soxs = make_pipeline_inputs(tmp_path)
    xshResults, xshSigmas = m.measure_spectra(xsh, lines)
    validation = m.validate_reference_lines(xshResults, lines)
    results, _ = m.build_analysis_results(soxs, lines[:-1], validation)
    monkeypatch.setattr(m, "OUTDIR", tmp_path / "out")
    monkeypatch.setattr(m, "OUTFILE_PREFIX", "TEST")

    linesTable, summaryTable = m.assemble_tables(results, xsh, xshResults, None, xshSigmas)
    m.write_tables(linesTable, summaryTable, validation)

    written = pd.read_csv(tmp_path / "out" / "TEST_summary.csv")
    assert set(written["instrument"]) == {"SOXS", "XSHOOTER"}
    assert (tmp_path / "out" / "TEST_lines.csv").exists()
    assert (tmp_path / "out" / "TEST_reference_lines.csv").exists()
    assert set(linesTable["version"]) == {"original", "shifted"}
    assert written.loc[written["instrument"] == "SOXS", "shift_nm"].iloc[0] == pytest.approx(-0.08, abs=0.01)


def test_outfile_prefix_derives_from_soxs_match_string():
    assert m.OUTFILE_PREFIX.startswith(m.SOXS_MATCH_STRING)


# =========================
# LOADERS
# =========================
def write_bad_fits(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("this is not a FITS file")


def test_find_soxs_spectra_records_unreadable_file_and_continues(tmp_path, monkeypatch):
    bad = tmp_path / "20260101T000000_NIR_SOXS_EXTRACTED_MERGED.fits"
    write_bad_fits(bad)
    monkeypatch.setattr(m, "SOXS_ROOT", tmp_path)

    spectra, failures = m.find_soxs_spectra()

    assert spectra == []
    assert [f["file"] for f in failures] == [bad]
    assert failures[0]["error"]


def test_find_xshooter_spectra_records_unreadable_file_and_continues(tmp_path, monkeypatch):
    bad = tmp_path / "obs1" / m.XSHOOTER_FILENAME
    write_bad_fits(bad)
    monkeypatch.setattr(m, "XSHOOTER_ROOT", tmp_path)

    spectra, failures = m.find_xshooter_spectra()

    assert spectra == []
    assert [f["file"] for f in failures] == [bad]


def test_loaders_do_not_swallow_unexpected_exceptions(tmp_path, monkeypatch):
    write_bad_fits(tmp_path / "obs1" / m.XSHOOTER_FILENAME)
    write_bad_fits(tmp_path / "a_EXTRACTED_MERGED.fits")
    monkeypatch.setattr(m, "XSHOOTER_ROOT", tmp_path)
    monkeypatch.setattr(m, "SOXS_ROOT", tmp_path)

    def explode(path):
        raise RuntimeError("bug")

    monkeypatch.setattr(m, "read_xshooter_spectrum", explode)
    monkeypatch.setattr(m, "read_soxs_spectrum", explode)

    with pytest.raises(RuntimeError):
        m.find_xshooter_spectra()
    with pytest.raises(RuntimeError):
        m.find_soxs_spectra()


def test_build_failed_row_marks_status_and_error():
    row = m.build_failed_row("SOXS", {"file": Path("dir/a.fits"), "error": "boom"})

    assert row["spectrum"] == "dir/a.fits"
    assert row["instrument"] == "SOXS"
    assert row["status"] == "failed"
    assert row["error"] == "boom"


def test_summary_table_includes_failed_files(tmp_path):
    lines, xsh, soxs = make_pipeline_inputs(tmp_path)
    xshResults, xshSigmas = m.measure_spectra(xsh, lines)
    validation = m.validate_reference_lines(xshResults, lines)
    results, _ = m.build_analysis_results(soxs, lines, validation)
    failures = [{"instrument": "XSHOOTER", "file": Path("obs9/x.fits"), "error": "bad header"}]

    _, summary = m.assemble_tables(results, xsh, xshResults, failures, xshSigmas)

    failed = summary.loc[summary["status"] == "failed"]
    assert failed["spectrum"].tolist() == ["obs9/x.fits"]
    assert set(summary.loc[summary["status"] == "ok", "instrument"]) == {"SOXS", "XSHOOTER"}


# =========================
# TREND AND PARABOLA
# =========================
def test_local_median_trend_excludes_the_line_itself():
    lines = np.array([0.0, 1.0, 2.0, 3.0])
    offsets = np.array([0.1, 0.1, 0.9, 0.9])

    trend = m.local_median_trend(lines, offsets, window_nm=20.0)

    assert trend[3] == pytest.approx(0.1)
    assert trend[0] == pytest.approx(0.9)


def test_local_median_trend_uses_global_median_with_fewer_than_three_neighbours():
    lines = np.array([0.0, 1.0, 2.0, 500.0])
    offsets = np.array([0.0, 0.0, 5.0, 1.0])

    trend = m.local_median_trend(lines, offsets, window_nm=20.0)

    assert np.allclose(trend, np.median(offsets))


def test_flag_outliers_catches_outlier_that_dominates_its_own_small_neighbourhood():
    lines = np.linspace(1700.0, 1800.0, 12)
    offsets = np.full(12, 0.01) + np.random.default_rng(8).normal(0, 0.002, 12)
    offsets[5] = 0.3

    assert m.flag_outliers(make_offsets_frame(lines, offsets)).tolist() == [i == 5 for i in range(12)]


def test_parabola_minimum_returns_nan_for_rank_deficient_fit_without_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        x0, _, _, _ = m.fit_local_parabola_minimum(np.array([1.0, 1.0, 1.0]), np.array([1.0, 2.0, 3.0]))

    assert np.isnan(x0)


def test_parabola_minimum_returns_nan_on_linalg_error(monkeypatch):
    def fail(*args, **kwargs):
        raise np.linalg.LinAlgError("no convergence")

    monkeypatch.setattr(m.np, "polyfit", fail)

    x0, _, _, _ = m.fit_local_parabola_minimum(np.array([0.0, 1.0, 2.0]), np.array([1.0, 0.0, 1.0]))

    assert np.isnan(x0)


def test_parabola_minimum_does_not_swallow_unexpected_errors(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("bug")

    monkeypatch.setattr(m.np, "polyfit", fail)

    with pytest.raises(RuntimeError):
        m.fit_local_parabola_minimum(np.array([0.0, 1.0, 2.0]), np.array([1.0, 0.0, 1.0]))


def test_records_parabola_fallback_when_fit_fails(monkeypatch):
    line = 1730.013
    wave, flux = make_spectrum([line], noise=0.0)
    monkeypatch.setattr(m, "fit_local_parabola_minimum", lambda x, y: (np.nan, np.nan, None, None))

    result = m.find_local_minimum_near_line(wave, flux, line)

    assert result["ok"]
    assert result["parabola_fallback"]
    assert result["wmin_nm"] in wave


def test_parabola_fallback_is_false_for_normal_and_unrefined_lines():
    line = 1730.013
    wave, flux = make_spectrum([line], noise=0.0)

    refined = m.find_local_minimum_near_line(wave, flux, line)
    discrete = m.find_local_minimum_near_line(wave, flux, line, parabola_refine=False)

    assert not refined["parabola_fallback"]
    assert not discrete["parabola_fallback"]


# =========================
# LINE WIDTH
# =========================
@pytest.mark.parametrize("sigma", [0.08, 0.12, 0.2])
def test_estimate_line_sigma_recovers_injected_sigma_within_ten_percent(sigma):
    lines = make_lines()
    wave, flux = make_spectrum(lines, depth=0.5, noise=0.003, sigma=sigma, seed=7)

    info = m.estimate_line_sigma(wave, flux, lines)

    assert info["sigma_nm"] == pytest.approx(sigma, rel=0.1)
    assert info["n_fits"] >= m.MIN_WIDTH_LINES
    assert info["n_failed_fits"] >= 0
    assert 0.0 <= info["sigma_mad_nm"] < 0.02


def test_estimate_line_sigma_falls_back_to_default_with_too_few_fits():
    lines = make_lines()[:3]
    wave, flux = make_spectrum(lines, depth=0.5, sigma=0.15)

    info = m.estimate_line_sigma(wave, flux, lines)

    assert info["sigma_nm"] == m.DEFAULT_LINE_SIGMA_NM
    assert info["n_fits"] < m.MIN_WIDTH_LINES


def test_estimate_line_sigma_ignores_lines_with_a_close_neighbour():
    lines = make_lines()
    wave, flux = make_spectrum(lines, depth=0.5, sigma=0.12, noise=0.003)
    crowded = np.concatenate([lines, lines + 0.3])

    isolated = m.estimate_line_sigma(wave, flux, lines)
    withNeighbours = m.estimate_line_sigma(wave, flux, lines, all_lines_nm=crowded)

    assert withNeighbours["n_fits"] == 0
    assert isolated["n_fits"] > 0


def test_estimate_line_sigma_ignores_lines_outside_the_depth_range():
    lines = make_lines()
    wave, flux = make_spectrum(lines, depth=0.95, sigma=0.12, noise=0.003)

    info = m.estimate_line_sigma(wave, flux, lines)

    assert info["n_fits"] == 0
    assert info["sigma_nm"] == m.DEFAULT_LINE_SIGMA_NM


def test_measure_spectra_returns_a_line_sigma_per_spectrum(tmp_path):
    lines, xsh, _ = make_pipeline_inputs(tmp_path)

    results, sigmas = m.measure_spectra(xsh, lines)

    assert len(results) == len(sigmas) == len(xsh)
    assert all(info["sigma_nm"] > 0 for info in sigmas)


def test_estimate_line_sigma_excludes_lines_with_a_one_sided_continuum():
    lines = make_lines()
    wave, flux = make_spectrum(lines, depth=0.5, sigma=0.12, noise=0.003)
    # REMOVE THE RIGHT SIDEBAND OF EVERY LINE SO THAT EACH HAS ONLY A ONE-SIDED CONTINUUM.
    gap = np.zeros(wave.size, dtype=bool)
    for line in lines:
        gap |= (wave > line + 0.1) & (wave < line + 0.6)

    info = m.estimate_line_sigma(wave[~gap], flux[~gap], lines)

    assert info["n_fits"] == 0
    assert info["sigma_nm"] == m.DEFAULT_LINE_SIGMA_NM


def test_estimate_line_sigma_reports_failed_fits():
    lines = make_lines()
    wave, flux = make_spectrum(lines, depth=0.5, sigma=0.12, noise=0.003)

    info = m.estimate_line_sigma(wave, flux, lines, min_lines=1)

    assert info["n_fits"] > 0
    assert info["n_failed_fits"] == 0
