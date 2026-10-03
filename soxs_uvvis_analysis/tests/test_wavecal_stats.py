import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.utils import iers

import wavecal_stats as ws


class TestSigmaClipMask:
    def test_masks_outlier_and_keeps_inliers(self):
        # ARRANGE
        rng = np.random.default_rng(1)
        values = np.append(rng.normal(0.0, 1.0, 50), 100.0)

        # ACT
        mask = ws.sigma_clip_mask(values)

        # ASSERT
        assert not mask[-1]
        assert mask[:-1].sum() >= 45

    def test_masks_nan_values(self):
        values = np.array([1.0, np.nan, 1.1, 0.9])

        mask = ws.sigma_clip_mask(values)

        assert mask.tolist() == [True, False, True, True]

    def test_keeps_all_finite_values_when_mad_is_zero(self):
        values = np.array([2.0, 2.0, 2.0, 2.0, 50.0, np.nan])

        mask = ws.sigma_clip_mask(values)

        assert mask.tolist() == [True, True, True, True, True, False]

    def test_returns_empty_mask_for_empty_input(self):
        mask = ws.sigma_clip_mask(np.array([]))

        assert mask.shape == (0,)
        assert mask.dtype == bool

    def test_does_not_modify_input(self):
        values = np.array([1.0, 2.0, 100.0, np.nan])
        original = values.copy()

        ws.sigma_clip_mask(values)

        np.testing.assert_array_equal(values, original)


SPEED_OF_LIGHT_KMS = 299792.458
RESULT_COLUMNS = [
    "Star",
    "Filename",
    "Order",
    "LineType",
    "LineWavelengthAngstrom",
    "ShiftAngstrom",
    "BaryCorrKms",
]
WAVELENGTH_RANGES = {
    "u": (3500.0, 4500.0),
    "g": (4500.0, 5500.0),
    "r": (5500.0, 7000.0),
    "i": (7000.0, 9000.0),
}
STAR_VELOCITIES = {"StarA": 25.0, "StarB": -40.0}
BARY_CORRS = {"a1": 12.0, "a2": -8.0, "b1": 3.0, "b2": -20.0}
EXPOSURE_STARS = {"a1": "StarA", "a2": "StarA", "b1": "StarB", "b2": "StarB"}
NOISE_SIGMA = 0.003
EXPOSURE_OFFSET_SIGMA = 0.02
N_LINES_PER_ORDER = 25


def make_results(seed=42, noiseSigma=NOISE_SIGMA):
    """BUILD A SYNTHETIC RESULTS FRAME AND RETURN IT WITH THE INJECTED ORDER OFFSETS"""
    rng = np.random.default_rng(seed)
    rows = []
    injectedOffsets = {}
    for filename, star in EXPOSURE_STARS.items():
        for order, (wmin, wmax) in WAVELENGTH_RANGES.items():
            offset = rng.normal(0.0, EXPOSURE_OFFSET_SIGMA)
            injectedOffsets[(filename, order)] = offset
            lineTypes = ["stellar"] + (["telluric"] if order in "ri" else [])
            for lineType in lineTypes:
                wavelengths = rng.uniform(wmin, wmax, N_LINES_PER_ORDER)
                shifts = offset + rng.normal(0.0, noiseSigma, N_LINES_PER_ORDER)
                if lineType == "stellar":
                    vShift = STAR_VELOCITIES[star] - BARY_CORRS[filename]
                    shifts = shifts + wavelengths * vShift / SPEED_OF_LIGHT_KMS
                rows.extend(
                    {
                        "Star": star,
                        "Filename": filename,
                        "Order": order,
                        "LineType": lineType,
                        "LineWavelengthAngstrom": w,
                        "ShiftAngstrom": sh,
                        "BaryCorrKms": BARY_CORRS[filename],
                    }
                    for w, sh in zip(wavelengths, shifts)
                )
    return pd.DataFrame(rows, columns=RESULT_COLUMNS), injectedOffsets


class TestTelluricOrderOffsets:
    def test_recovers_injected_offsets_per_exposure_and_order(self):
        resultsDf, injected = make_results()

        offsetsDf = ws.telluric_order_offsets(resultsDf)

        assert list(offsetsDf.columns) == [
            "Filename",
            "Order",
            "OrderOffsetAngstrom",
            "NTelluric",
        ]
        assert len(offsetsDf) == 8
        for row in offsetsDf.itertuples():
            assert row.OrderOffsetAngstrom == pytest.approx(
                injected[(row.Filename, row.Order)], abs=0.004
            )
            assert 0.8 * N_LINES_PER_ORDER <= row.NTelluric <= N_LINES_PER_ORDER

    def test_only_r_and_i_orders_have_offsets(self):
        resultsDf, _ = make_results()

        offsetsDf = ws.telluric_order_offsets(resultsDf)

        assert set(offsetsDf["Order"]) == {"r", "i"}

    def test_clips_outlier_and_ignores_nan_shifts(self):
        resultsDf, injected = make_results()
        telluricIdx = resultsDf.index[
            (resultsDf["LineType"] == "telluric")
            & (resultsDf["Filename"] == "a1")
            & (resultsDf["Order"] == "r")
        ]
        resultsDf.loc[telluricIdx[0], "ShiftAngstrom"] = 5.0
        resultsDf.loc[telluricIdx[1], "ShiftAngstrom"] = np.nan

        offsetsDf = ws.telluric_order_offsets(resultsDf)

        row = offsetsDf[(offsetsDf["Filename"] == "a1") & (offsetsDf["Order"] == "r")]
        assert row["OrderOffsetAngstrom"].iloc[0] == pytest.approx(
            injected[("a1", "r")], abs=0.004
        )
        assert row["NTelluric"].iloc[0] <= N_LINES_PER_ORDER - 2

    def test_drops_exposure_whose_offset_is_an_outlier_across_exposures(self):
        # ARRANGE: ONE TELLURIC LINE PER EXPOSURE, SO A BAD FIT CANNOT BE CLIPPED INSIDE ITS EXPOSURE
        rng = np.random.default_rng(3)
        shifts = list(rng.normal(0.0, 0.02, 9)) + [1.5]
        resultsDf = pd.DataFrame(
            {
                "Star": "StarA",
                "Filename": [f"e{i}" for i in range(10)],
                "Order": "r",
                "LineType": "telluric",
                "LineWavelengthAngstrom": 6277.7,
                "ShiftAngstrom": shifts,
                "BaryCorrKms": 0.0,
            }
        )

        # ACT
        offsetsDf = ws.telluric_order_offsets(resultsDf)

        # ASSERT
        assert "e9" not in set(offsetsDf["Filename"])
        assert len(offsetsDf) == 9

    def test_keeps_all_exposures_when_too_few_to_clip_across_exposures(self):
        resultsDf = pd.DataFrame(
            {
                "Star": "StarA",
                "Filename": ["e0", "e1", "e2"],
                "Order": "r",
                "LineType": "telluric",
                "LineWavelengthAngstrom": 6277.7,
                "ShiftAngstrom": [0.0, 0.01, 1.5],
                "BaryCorrKms": 0.0,
            }
        )

        offsetsDf = ws.telluric_order_offsets(resultsDf)

        assert len(offsetsDf) == 3

    def test_does_not_mutate_input(self):
        resultsDf, _ = make_results()
        original = resultsDf.copy()

        ws.telluric_order_offsets(resultsDf)

        pd.testing.assert_frame_equal(resultsDf, original)

    def test_returns_empty_frame_with_columns_for_empty_input(self):
        emptyDf = pd.DataFrame(columns=RESULT_COLUMNS)

        offsetsDf = ws.telluric_order_offsets(emptyDf)

        assert offsetsDf.empty
        assert list(offsetsDf.columns) == [
            "Filename",
            "Order",
            "OrderOffsetAngstrom",
            "NTelluric",
        ]


VELOCITY_COLUMNS = ["Star", "VStarKms", "VStarErrKms", "VStarMadStdKms", "NLines"]


def add_star_without_telluric_orders(resultsDf):
    """APPEND A STAR WITH ONLY U AND G STELLAR LINES"""
    extraDf = pd.DataFrame(
        {
            "Star": "StarC",
            "Filename": "c1",
            "Order": ["u", "g", "u", "g"],
            "LineType": "stellar",
            "LineWavelengthAngstrom": [4000.0, 5000.0, 4100.0, 5100.0],
            "ShiftAngstrom": [0.01, 0.02, 0.01, 0.02],
            "BaryCorrKms": 5.0,
        }
    )
    return pd.concat([resultsDf, extraDf], ignore_index=True)


class TestFitStarVelocities:
    def test_recovers_injected_star_velocities(self):
        resultsDf, _ = make_results()
        offsetsDf = ws.telluric_order_offsets(resultsDf)

        velocitiesDf = ws.fit_star_velocities(resultsDf, offsetsDf)

        assert list(velocitiesDf.columns) == VELOCITY_COLUMNS
        assert set(velocitiesDf["Star"]) == set(STAR_VELOCITIES)
        for row in velocitiesDf.itertuples():
            assert row.VStarKms == pytest.approx(STAR_VELOCITIES[row.Star], abs=0.2)
            assert row.VStarErrKms == pytest.approx(
                row.VStarMadStdKms / np.sqrt(row.NLines)
            )
            assert row.VStarMadStdKms > 0

    def test_uses_only_stellar_lines_in_orders_with_telluric_offsets(self):
        resultsDf, _ = make_results()
        offsetsDf = ws.telluric_order_offsets(resultsDf)

        velocitiesDf = ws.fit_star_velocities(resultsDf, offsetsDf)

        # 2 EXPOSURES x 2 ORDERS (R, I) x 25 LINES, MINUS ANY CLIPPED
        assert (velocitiesDf["NLines"] <= 100).all()
        assert (velocitiesDf["NLines"] >= 80).all()

    def test_omits_star_without_usable_orders(self):
        resultsDf, _ = make_results()
        resultsDf = add_star_without_telluric_orders(resultsDf)
        offsetsDf = ws.telluric_order_offsets(resultsDf)

        velocitiesDf = ws.fit_star_velocities(resultsDf, offsetsDf)

        assert "StarC" not in set(velocitiesDf["Star"])

    def test_clips_outlier_and_ignores_nan_shifts(self):
        resultsDf, _ = make_results()
        offsetsDf = ws.telluric_order_offsets(resultsDf)
        stellarIdx = resultsDf.index[
            (resultsDf["LineType"] == "stellar")
            & (resultsDf["Star"] == "StarA")
            & (resultsDf["Order"] == "r")
        ]
        resultsDf.loc[stellarIdx[0], "ShiftAngstrom"] = 3.0
        resultsDf.loc[stellarIdx[1], "ShiftAngstrom"] = np.nan

        velocitiesDf = ws.fit_star_velocities(resultsDf, offsetsDf)

        starA = velocitiesDf[velocitiesDf["Star"] == "StarA"].iloc[0]
        assert starA["VStarKms"] == pytest.approx(STAR_VELOCITIES["StarA"], abs=0.2)

    def test_does_not_mutate_inputs(self):
        resultsDf, _ = make_results()
        offsetsDf = ws.telluric_order_offsets(resultsDf)
        resultsCopy = resultsDf.copy()
        offsetsCopy = offsetsDf.copy()

        ws.fit_star_velocities(resultsDf, offsetsDf)

        pd.testing.assert_frame_equal(resultsDf, resultsCopy)
        pd.testing.assert_frame_equal(offsetsDf, offsetsCopy)

    def test_returns_empty_frame_with_columns_for_empty_input(self):
        emptyDf = pd.DataFrame(columns=RESULT_COLUMNS)
        offsetsDf = ws.telluric_order_offsets(emptyDf)

        velocitiesDf = ws.fit_star_velocities(emptyDf, offsetsDf)

        assert velocitiesDf.empty
        assert list(velocitiesDf.columns) == VELOCITY_COLUMNS


def test_fit_star_velocities_omits_star_whose_shifts_are_all_nan():
    # ARRANGE
    resultsDf, _ = make_results()
    resultsDf.loc[
        (resultsDf["Star"] == "StarB") & (resultsDf["LineType"] == "stellar"),
        "ShiftAngstrom",
    ] = np.nan
    offsetsDf = ws.telluric_order_offsets(resultsDf)

    # ACT
    velocitiesDf = ws.fit_star_velocities(resultsDf, offsetsDf)

    # ASSERT
    assert velocitiesDf["Star"].tolist() == ["StarA"]


class TestVelocityCorrectedResiduals:
    def fit(self, resultsDf):
        offsetsDf = ws.telluric_order_offsets(resultsDf)
        return ws.fit_star_velocities(resultsDf, offsetsDf)

    def test_stellar_residuals_in_u_and_g_match_injected_order_offsets(self):
        resultsDf, injected = make_results()
        velocitiesDf = self.fit(resultsDf)

        residualsDf = ws.velocity_corrected_residuals(resultsDf, velocitiesDf)

        stellarDf = residualsDf[
            (residualsDf["LineType"] == "stellar")
            & residualsDf["Order"].isin(["u", "g"])
        ]
        medians = stellarDf.groupby(["Filename", "Order"])["ResidualAngstrom"].median()
        assert len(medians) == 8
        for (filename, order), median in medians.items():
            assert median == pytest.approx(injected[(filename, order)], abs=0.005)

    def test_telluric_expected_shift_is_zero_and_residual_equals_shift(self):
        resultsDf, _ = make_results()
        velocitiesDf = self.fit(resultsDf)

        residualsDf = ws.velocity_corrected_residuals(resultsDf, velocitiesDf)

        telluricDf = residualsDf[residualsDf["LineType"] == "telluric"]
        assert (telluricDf["ExpectedVelocityShiftAngstrom"] == 0).all()
        np.testing.assert_allclose(
            telluricDf["ResidualAngstrom"], telluricDf["ShiftAngstrom"]
        )

    def test_expected_stellar_shift_follows_velocity_and_barycentric_correction(self):
        resultsDf, _ = make_results()
        velocitiesDf = pd.DataFrame(
            {
                "Star": ["StarA", "StarB"],
                "VStarKms": [25.0, -40.0],
            }
        )

        residualsDf = ws.velocity_corrected_residuals(resultsDf, velocitiesDf)

        row = residualsDf[
            (residualsDf["LineType"] == "stellar") & (residualsDf["Filename"] == "b2")
        ].iloc[0]
        expected = row["LineWavelengthAngstrom"] * (-40.0 - -20.0) / SPEED_OF_LIGHT_KMS
        assert row["ExpectedVelocityShiftAngstrom"] == pytest.approx(expected)
        assert row["ResidualAngstrom"] == pytest.approx(row["ShiftAngstrom"] - expected)

    def test_star_without_velocity_gets_nan_stellar_residuals(self):
        resultsDf, _ = make_results()
        resultsDf = add_star_without_telluric_orders(resultsDf)
        velocitiesDf = self.fit(resultsDf)

        residualsDf = ws.velocity_corrected_residuals(resultsDf, velocitiesDf)

        starCDf = residualsDf[residualsDf["Star"] == "StarC"]
        assert len(starCDf) == 4
        assert starCDf["ResidualAngstrom"].isna().all()

    def test_does_not_mutate_inputs_and_returns_new_frame(self):
        resultsDf, _ = make_results()
        velocitiesDf = self.fit(resultsDf)
        resultsCopy = resultsDf.copy()
        velocitiesCopy = velocitiesDf.copy()

        residualsDf = ws.velocity_corrected_residuals(resultsDf, velocitiesDf)

        pd.testing.assert_frame_equal(resultsDf, resultsCopy)
        pd.testing.assert_frame_equal(velocitiesDf, velocitiesCopy)
        assert residualsDf is not resultsDf
        assert len(residualsDf) == len(resultsDf)

    def test_returns_empty_frame_with_new_columns_for_empty_input(self):
        emptyDf = pd.DataFrame(columns=RESULT_COLUMNS)
        velocitiesDf = self.fit(emptyDf)

        residualsDf = ws.velocity_corrected_residuals(emptyDf, velocitiesDf)

        assert residualsDf.empty
        assert list(residualsDf.columns) == RESULT_COLUMNS + [
            "ExpectedVelocityShiftAngstrom",
            "ResidualAngstrom",
        ]


STATS_COLUMNS = [
    "Order",
    "NLines",
    "NExposures",
    "MedianShiftAngstrom",
    "MeanShiftAngstrom",
    "StdShiftAngstrom",
    "SEMeanAngstrom",
    "InternalScatterStd",
    "InternalScatterMadStd",
    "ExposureScatterStd",
    "NStars",
    "StarScatterStd",
    "NClipped",
]
INTERNAL_SIGMA = 0.005
EXPOSURE_SIGMA = 0.02
N_STATS_EXPOSURES = 30
N_STATS_LINES = 60


def make_residuals(orders=("u",), lineType="stellar", seed=7):
    """BUILD A RESIDUALS FRAME WITH KNOWN INTERNAL AND EXPOSURE-TO-EXPOSURE SCATTER"""
    rng = np.random.default_rng(seed)
    rows = []
    for order in orders:
        for exposure in range(N_STATS_EXPOSURES):
            offset = rng.normal(0.0, EXPOSURE_SIGMA)
            for residual in offset + rng.normal(0.0, INTERNAL_SIGMA, N_STATS_LINES):
                rows.append(
                    {
                        "Filename": f"exp{exposure}",
                        "Order": order,
                        "LineType": lineType,
                        "ResidualAngstrom": residual,
                    }
                )
    return pd.DataFrame(rows)


class TestClippedOrderStats:
    def test_internal_scatter_is_unbiased_with_two_lines_per_exposure(self):
        # ARRANGE
        rng = np.random.default_rng(11)
        sigma = 0.1
        rows = [
            {
                "Filename": f"exp{e}",
                "Order": "u",
                "LineType": "stellar",
                "ResidualAngstrom": v,
            }
            for e in range(400)
            for v in rng.normal(rng.normal(0.0, 0.5), sigma, 2)
        ]

        # ACT
        row = ws.clipped_order_stats(pd.DataFrame(rows)).iloc[0]

        # ASSERT
        assert row["InternalScatterStd"] == pytest.approx(sigma, rel=0.1)
        assert row["InternalScatterMadStd"] == pytest.approx(sigma, rel=0.15)

    def test_internal_scatter_is_nan_with_one_line_per_exposure(self):
        # ARRANGE
        residualsDf = make_residuals()
        oneLineDf = residualsDf.groupby("Filename").head(1)

        # ACT
        row = ws.clipped_order_stats(oneLineDf).iloc[0]

        # ASSERT
        assert np.isnan(row["InternalScatterStd"])
        assert np.isnan(row["InternalScatterMadStd"])
        assert np.isfinite(row["ExposureScatterStd"])

    def test_star_scatter_is_std_of_per_star_medians(self):
        # ARRANGE
        residualsDf = make_residuals()
        starOffsets = {"S1": -0.3, "S2": 0.0, "S3": 0.3}
        starNames = [list(starOffsets)[int(f[3:]) % 3] for f in residualsDf["Filename"]]
        residualsDf = residualsDf.assign(Star=starNames)
        residualsDf = residualsDf.assign(
            ResidualAngstrom=residualsDf["ResidualAngstrom"]
            + residualsDf["Star"].map(starOffsets)
        )

        # ACT
        row = ws.clipped_order_stats(residualsDf).iloc[0]

        # ASSERT
        assert row["NStars"] == 3
        assert row["StarScatterStd"] == pytest.approx(0.3, abs=0.03)

    def test_ignores_rows_without_a_filename(self):
        # ARRANGE
        residualsDf = make_residuals()
        withMissingDf = pd.concat(
            [
                residualsDf,
                pd.DataFrame(
                    [
                        {
                            "Filename": np.nan,
                            "Order": "u",
                            "LineType": "stellar",
                            "ResidualAngstrom": 0.0,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        # ACT
        cleanRow = ws.clipped_order_stats(residualsDf).iloc[0]
        missingRow = ws.clipped_order_stats(withMissingDf).iloc[0]

        # ASSERT
        assert missingRow["NLines"] == cleanRow["NLines"]
        assert missingRow["InternalScatterStd"] == pytest.approx(
            cleanRow["InternalScatterStd"]
        )
        assert missingRow["InternalScatterMadStd"] == pytest.approx(
            cleanRow["InternalScatterMadStd"]
        )

    def test_star_columns_are_nan_without_a_star_column(self):
        row = ws.clipped_order_stats(make_residuals()).iloc[0]

        assert np.isnan(row["NStars"])
        assert np.isnan(row["StarScatterStd"])

    def test_scatter_matches_injected_internal_and_exposure_sigmas(self):
        residualsDf = make_residuals()

        statsDf = ws.clipped_order_stats(residualsDf)

        assert list(statsDf.columns) == STATS_COLUMNS
        row = statsDf.iloc[0]
        assert row["Order"] == "u"
        assert row["NExposures"] == N_STATS_EXPOSURES
        assert row["InternalScatterStd"] == pytest.approx(INTERNAL_SIGMA, rel=0.15)
        assert row["InternalScatterMadStd"] == pytest.approx(INTERNAL_SIGMA, rel=0.15)
        assert row["ExposureScatterStd"] == pytest.approx(EXPOSURE_SIGMA, rel=0.35)

    def test_reports_median_mean_std_and_standard_error_of_clipped_values(self):
        residualsDf = make_residuals()

        row = ws.clipped_order_stats(residualsDf).iloc[0]

        values = residualsDf["ResidualAngstrom"].to_numpy()
        keptValues = values[ws.sigma_clip_mask(values, maxIters=1)]
        assert row["NLines"] == keptValues.size
        assert row["MedianShiftAngstrom"] == pytest.approx(np.median(keptValues))
        assert row["MeanShiftAngstrom"] == pytest.approx(keptValues.mean())
        assert row["StdShiftAngstrom"] == pytest.approx(keptValues.std(ddof=1))
        # THE STANDARD ERROR COMES FROM THE PER-EXPOSURE MEDIANS, NOT FROM INDEPENDENT LINES
        assert row["SEMeanAngstrom"] == pytest.approx(
            row["ExposureScatterStd"] / np.sqrt(row["NExposures"])
        )
        assert row["NClipped"] == values.size - keptValues.size

    def test_clips_injected_outlier_without_moving_median_or_mean(self):
        cleanDf = make_residuals()
        outlierDf = cleanDf.copy()
        outlierDf.loc[0, "ResidualAngstrom"] = 10.0

        cleanRow = ws.clipped_order_stats(cleanDf).iloc[0]
        outlierRow = ws.clipped_order_stats(outlierDf).iloc[0]

        assert outlierRow["NClipped"] >= 1
        assert outlierRow["MedianShiftAngstrom"] == pytest.approx(
            cleanRow["MedianShiftAngstrom"], abs=0.001
        )
        assert outlierRow["MeanShiftAngstrom"] == pytest.approx(
            cleanRow["MeanShiftAngstrom"], abs=0.001
        )

    def test_ignores_nan_values(self):
        residualsDf = make_residuals()
        residualsDf.loc[:4, "ResidualAngstrom"] = np.nan

        row = ws.clipped_order_stats(residualsDf).iloc[0]

        assert row["NLines"] <= len(residualsDf) - 5
        assert np.isfinite(row["MedianShiftAngstrom"])

    def test_orders_rows_u_g_r_i_then_others(self):
        residualsDf = make_residuals(orders=("i", "x", "r", "g", "u"))

        statsDf = ws.clipped_order_stats(residualsDf)

        assert statsDf["Order"].tolist() == ["u", "g", "r", "i", "x"]

    def test_anchor_by_order_keeps_only_anchor_line_type_and_records_it(self):
        stellarDf = make_residuals(orders=("u", "r"), lineType="stellar")
        telluricDf = make_residuals(orders=("u", "r"), lineType="telluric", seed=8)
        telluricDf["ResidualAngstrom"] += 1.0
        combinedDf = pd.concat([stellarDf, telluricDf], ignore_index=True)

        statsDf = ws.clipped_order_stats(
            combinedDf, anchorByOrder={"u": "stellar", "r": "telluric"}
        )

        byOrder = statsDf.set_index("Order")
        assert byOrder.loc["u", "Anchor"] == "stellar"
        assert byOrder.loc["r", "Anchor"] == "telluric"
        assert abs(byOrder.loc["u", "MedianShiftAngstrom"]) < 0.02
        assert abs(byOrder.loc["r", "MedianShiftAngstrom"] - 1.0) < 0.02
        assert byOrder.loc["u", "NLines"] <= N_STATS_EXPOSURES * N_STATS_LINES

    def test_orders_missing_from_anchor_map_keep_all_line_types(self):
        residualsDf = make_residuals(orders=("g",))

        statsDf = ws.clipped_order_stats(residualsDf, anchorByOrder={"u": "stellar"})

        assert statsDf.iloc[0]["Anchor"] == "all"
        assert statsDf.iloc[0]["NLines"] > 0

    def test_uses_requested_value_column(self):
        residualsDf = make_residuals().rename(
            columns={"ResidualAngstrom": "ShiftAngstrom"}
        )

        row = ws.clipped_order_stats(residualsDf, valueColumn="ShiftAngstrom").iloc[0]

        assert row["NLines"] > 0

    def test_exposure_scatter_is_nan_for_single_exposure(self):
        residualsDf = make_residuals()
        residualsDf = residualsDf[residualsDf["Filename"] == "exp0"]

        row = ws.clipped_order_stats(residualsDf).iloc[0]

        assert row["NExposures"] == 1
        assert np.isnan(row["ExposureScatterStd"])
        assert np.isfinite(row["InternalScatterStd"])

    def test_does_not_mutate_input(self):
        residualsDf = make_residuals()
        original = residualsDf.copy()

        ws.clipped_order_stats(residualsDf, anchorByOrder={"u": "stellar"})

        pd.testing.assert_frame_equal(residualsDf, original)

    def test_returns_empty_frame_with_columns_for_empty_input(self):
        emptyDf = pd.DataFrame(
            columns=["Filename", "Order", "LineType", "ResidualAngstrom"]
        )

        statsDf = ws.clipped_order_stats(emptyDf)

        assert statsDf.empty
        assert list(statsDf.columns) == STATS_COLUMNS

    def test_empty_input_includes_anchor_column_when_anchors_given(self):
        emptyDf = pd.DataFrame(
            columns=["Filename", "Order", "LineType", "ResidualAngstrom"]
        )

        statsDf = ws.clipped_order_stats(emptyDf, anchorByOrder={"u": "stellar"})

        assert statsDf.empty
        assert list(statsDf.columns) == ["Order", "Anchor"] + STATS_COLUMNS[1:]


# ECLIPTIC TARGET AT RA=0, DEC=0 AND DATES NEAR THE JUNE AND DECEMBER SOLSTICES
JUNE_MJD = 60116.0
DECEMBER_MJD = 60299.0
MAX_EARTH_VELOCITY_KMS = 35.0


@pytest.fixture
def offlineIers():
    """STOP ASTROPY FROM DOWNLOADING EARTH-ORIENTATION DATA DURING THE TESTS"""
    with iers.conf.set_temp("auto_download", False):
        yield


def make_header(ra=0.0, dec=0.0, mjd=JUNE_MJD):
    header = fits.Header()
    header["RA"] = ra
    header["DEC"] = dec
    header["MJD-OBS"] = mjd
    return header


class TestBarycentricCorrectionKms:
    def test_returns_finite_float_within_earth_orbital_speed(self, offlineIers):
        header = make_header()

        correction = ws.barycentric_correction_kms(header)

        assert isinstance(correction, float)
        assert np.isfinite(correction)
        assert abs(correction) < MAX_EARTH_VELOCITY_KMS

    def test_changes_sign_between_dates_six_months_apart(self, offlineIers):
        juneCorrection = ws.barycentric_correction_kms(make_header(mjd=JUNE_MJD))
        decemberCorrection = ws.barycentric_correction_kms(
            make_header(mjd=DECEMBER_MJD)
        )

        assert juneCorrection * decemberCorrection < 0
        assert abs(juneCorrection) > 20.0
        assert abs(decemberCorrection) > 20.0

    def test_accepts_explicit_location(self, offlineIers):
        import astropy.units as u
        from astropy.coordinates import EarthLocation

        header = make_header()
        otherSite = EarthLocation.from_geodetic(
            lon=-155.47 * u.deg, lat=19.82 * u.deg, height=4200 * u.m
        )

        defaultCorrection = ws.barycentric_correction_kms(header)
        otherCorrection = ws.barycentric_correction_kms(header, location=otherSite)

        assert otherCorrection != defaultCorrection
        # EARTH ROTATION ADDS UP TO ~0.46 KM/S AT EACH SITE, SO TWO SITES DIFFER BY < 1 KM/S
        assert otherCorrection == pytest.approx(defaultCorrection, abs=1.0)

    @pytest.mark.parametrize("missingKey", ["RA", "DEC", "MJD-OBS"])
    def test_raises_value_error_naming_missing_key(self, missingKey):
        header = make_header()
        del header[missingKey]

        with pytest.raises(ValueError, match=missingKey):
            ws.barycentric_correction_kms(header)
