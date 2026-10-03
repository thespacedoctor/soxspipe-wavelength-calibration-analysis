"""*statistics helpers for the UV-VIS wavelength-calibration checks*"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Mapping

    import numpy as np
    import pandas as pd
    from astropy.coordinates import EarthLocation

SPEED_OF_LIGHT_KMS = 299792.458
MAD_TO_STD = 1.4826
# NTT (LA SILLA) SITE, HOST TELESCOPE OF SOXS
NTT_LONGITUDE_DEG = -70.7345
NTT_LATITUDE_DEG = -29.2584
NTT_HEIGHT_M = 2375.0
VELOCITY_COLUMNS = ["Star", "VStarKms", "VStarErrKms", "VStarMadStdKms", "NLines"]
ORDER_SEQUENCE = ["u", "g", "r", "i"]
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
OFFSET_COLUMNS = ["Filename", "Order", "OrderOffsetAngstrom", "NTelluric"]
# FEWEST EXPOSURES NEEDED BEFORE PER-EXPOSURE OFFSETS ARE CLIPPED AGAINST EACH OTHER
MIN_EXPOSURES_FOR_CLIPPING = 5


def sigma_clip_mask(
    values: np.ndarray | pd.Series | list, nSigma: float = 3.0, maxIters: int = 5
) -> np.ndarray:
    """*iteratively clip values about the median using a MAD-based standard deviation*

    **Key Arguments:**

    - ``values`` -- array-like of values to clip. NaNs are always masked out.
    - ``nSigma`` -- clipping threshold in units of the MAD-based standard deviation. Default *3.0*
    - ``maxIters`` -- maximum number of clipping iterations. Default *5*

    **Return:**

    - ``mask`` -- boolean array, True for the values that are kept

    **Usage:**

    ```
    mask = sigma_clip_mask(values, nSigma=3.0)
    ```
    """
    import numpy as np

    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    for _ in range(maxIters):
        kept = values[mask]
        if kept.size == 0:
            break
        median = np.median(kept)
        madStd = MAD_TO_STD * np.median(np.abs(kept - median))
        # A ZERO MAD MEANS NO SCALE IS AVAILABLE SO KEEP EVERYTHING
        if madStd == 0:
            break
        newMask = mask & (np.abs(values - median) <= nSigma * madStd)
        if newMask.sum() == mask.sum():
            break
        mask = newMask
    return mask


def _clipped_values(values: np.ndarray | pd.Series, nSigma: float = 3.0) -> np.ndarray:
    """*return the finite, sigma-clipped values as a float array*"""
    import numpy as np

    values = np.asarray(values, dtype=float)
    return values[sigma_clip_mask(values, nSigma=nSigma)]


def telluric_order_offsets(resultsDf: pd.DataFrame) -> pd.DataFrame:
    """*measure the per-exposure, per-order offset from the telluric lines*

    **Key Arguments:**

    - ``resultsDf`` -- per-line results DataFrame (needs ``Filename``, ``Order``, ``LineType`` and ``ShiftAngstrom``)

    **Return:**

    - ``offsetsDf`` -- DataFrame with columns ``Filename``, ``Order``, ``OrderOffsetAngstrom`` (clipped median telluric shift) and ``NTelluric`` (number of lines kept). When an order has at least ``MIN_EXPOSURES_FOR_CLIPPING`` exposures, an exposure whose offset is an outlier against the other exposures is dropped. This catches bad fits in orders with only one telluric line per exposure.

    **Usage:**

    ```
    offsetsDf = telluric_order_offsets(resultsDf)
    ```
    """
    import numpy as np
    import pandas as pd

    telluricDf = resultsDf[resultsDf["LineType"] == "telluric"]
    rows = []
    for (filename, order), groupDf in telluricDf.groupby(["Filename", "Order"]):
        kept = _clipped_values(groupDf["ShiftAngstrom"])
        if kept.size == 0:
            continue
        rows.append(
            {
                "Filename": filename,
                "Order": order,
                "OrderOffsetAngstrom": float(np.median(kept)),
                "NTelluric": int(kept.size),
            }
        )
    offsetsDf = pd.DataFrame(rows, columns=OFFSET_COLUMNS)

    # A SINGLE BAD TELLURIC FIT CANNOT BE CLIPPED INSIDE ITS OWN EXPOSURE, SO CLIP ACROSS EXPOSURES
    keepMask = np.ones(len(offsetsDf), dtype=bool)
    for _, groupDf in offsetsDf.groupby("Order"):
        if len(groupDf) >= MIN_EXPOSURES_FOR_CLIPPING:
            positions = offsetsDf.index.get_indexer(groupDf.index)
            # ONE PASS ONLY: ITERATING ON A FEW EXPOSURES SHRINKS THE MAD AND REJECTS GOOD ONES
            keepMask[positions] = sigma_clip_mask(
                groupDf["OrderOffsetAngstrom"], maxIters=1
            )
    return offsetsDf[keepMask].reset_index(drop=True)


def _mad_std(values: np.ndarray) -> float:
    """*MAD-based standard deviation of an array (NaN if empty)*"""
    import numpy as np

    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.nan
    return float(MAD_TO_STD * np.median(np.abs(values - np.median(values))))


def fit_star_velocities(
    resultsDf: pd.DataFrame, orderOffsetsDf: pd.DataFrame
) -> pd.DataFrame:
    """*fit the barycentric radial velocity of each star from its telluric-corrected stellar lines*

    **Key Arguments:**

    - ``resultsDf`` -- per-line results DataFrame
    - ``orderOffsetsDf`` -- output of ``telluric_order_offsets``. Only stellar lines in an exposure and order with a telluric offset are used.

    **Return:**

    - ``velocitiesDf`` -- DataFrame with columns ``Star``, ``VStarKms`` (clipped median), ``VStarErrKms`` (MAD-std / sqrt(N); a lower limit, because lines in one exposure share an offset), ``VStarMadStdKms`` and ``NLines``. Stars with no usable lines are omitted.

    **Usage:**

    ```
    velocitiesDf = fit_star_velocities(resultsDf, orderOffsetsDf)
    ```
    """
    import numpy as np
    import pandas as pd

    stellarDf = resultsDf[resultsDf["LineType"] == "stellar"].merge(
        orderOffsetsDf[["Filename", "Order", "OrderOffsetAngstrom"]],
        on=["Filename", "Order"],
        how="inner",
    )
    correctedShift = stellarDf["ShiftAngstrom"] - stellarDf["OrderOffsetAngstrom"]
    vTopo = SPEED_OF_LIGHT_KMS * correctedShift / stellarDf["LineWavelengthAngstrom"]
    stellarDf = stellarDf.assign(VStarKms=vTopo + stellarDf["BaryCorrKms"])
    rows = []
    for star, groupDf in stellarDf.groupby("Star"):
        kept = _clipped_values(groupDf["VStarKms"])
        if kept.size == 0:
            continue
        madStd = _mad_std(kept)
        rows.append(
            {
                "Star": star,
                "VStarKms": float(np.median(kept)),
                "VStarErrKms": madStd / np.sqrt(kept.size),
                "VStarMadStdKms": madStd,
                "NLines": int(kept.size),
            }
        )
    return pd.DataFrame(rows, columns=VELOCITY_COLUMNS)


def velocity_corrected_residuals(
    resultsDf: pd.DataFrame, starVelocitiesDf: pd.DataFrame
) -> pd.DataFrame:
    """*subtract the expected stellar velocity shift from every line shift*

    **Key Arguments:**

    - ``resultsDf`` -- per-line results DataFrame
    - ``starVelocitiesDf`` -- output of ``fit_star_velocities`` (needs ``Star`` and ``VStarKms``)

    **Return:**

    - ``residualsDf`` -- a new DataFrame with the extra columns ``ExpectedVelocityShiftAngstrom`` (wavelength * (vStar - baryCorr) / c for stellar lines, 0 for telluric lines) and ``ResidualAngstrom`` (shift minus expected shift). Stellar lines of stars without a velocity get NaN.

    **Usage:**

    ```
    residualsDf = velocity_corrected_residuals(resultsDf, starVelocitiesDf)
    ```
    """
    vStarByStar = starVelocitiesDf.set_index("Star")["VStarKms"]
    vStar = resultsDf["Star"].map(vStarByStar).astype(float)
    stellarExpected = (
        resultsDf["LineWavelengthAngstrom"].astype(float)
        * (vStar - resultsDf["BaryCorrKms"].astype(float))
        / SPEED_OF_LIGHT_KMS
    )
    isTelluric = resultsDf["LineType"] == "telluric"
    expected = stellarExpected.mask(isTelluric, 0.0)
    return resultsDf.assign(
        ExpectedVelocityShiftAngstrom=expected,
        ResidualAngstrom=resultsDf["ShiftAngstrom"].astype(float) - expected,
    )


def _std_ddof1(values: np.ndarray) -> float:
    """*sample standard deviation (ddof=1), NaN for fewer than two values*"""
    import numpy as np

    values = np.asarray(values, dtype=float)
    return float(np.std(values, ddof=1)) if values.size > 1 else np.nan


def _internal_scatter(keptDf: pd.DataFrame, valueColumn: str) -> tuple[float, float]:
    """*scatter of the lines about their own exposure's median, corrected for the fitted medians*"""
    import numpy as np

    groupSizes = (
        keptDf.groupby("Filename")[valueColumn].transform("size").to_numpy(dtype=float)
    )
    deviations = (
        keptDf[valueColumn]
        - keptDf.groupby("Filename")[valueColumn].transform("median")
    ).to_numpy(dtype=float)
    # EACH EXPOSURE MEDIAN USES ONE DEGREE OF FREEDOM; EXPOSURES WITH ONE LINE ADD NO INFORMATION
    degreesOfFreedom = deviations.size - keptDf["Filename"].nunique()
    if degreesOfFreedom <= 0:
        return np.nan, np.nan
    internalStd = float(np.sqrt(np.sum(deviations**2) / degreesOfFreedom))
    multiLine = groupSizes > 1
    scaledDeviations = deviations[multiLine] * np.sqrt(
        groupSizes[multiLine] / (groupSizes[multiLine] - 1)
    )
    internalMadStd = float(MAD_TO_STD * np.median(np.abs(scaledDeviations)))
    return internalStd, internalMadStd


def _star_scatter(keptDf: pd.DataFrame, valueColumn: str) -> tuple[float, float]:
    """*number of stars and the scatter of the per-star medians (NaN without a Star column)*"""
    import numpy as np

    if "Star" not in keptDf.columns:
        return np.nan, np.nan
    starMedians = keptDf.groupby("Star")[valueColumn].median().to_numpy(dtype=float)
    return float(starMedians.size), _std_ddof1(starMedians)


def _order_stats_row(order: str, groupDf: pd.DataFrame, valueColumn: str) -> dict:
    """*compute the clipped statistics for the lines of a single order*"""
    import numpy as np

    # A LINE WITHOUT AN EXPOSURE CANNOT BE PLACED IN THE PER-EXPOSURE STATISTICS
    groupDf = groupDf.dropna(subset=["Filename"])
    values = groupDf[valueColumn].to_numpy(dtype=float)
    # ONE 3-SIGMA MAD PASS, THE SAME AS THE CROSS-EXPOSURE CLIP IN telluric_order_offsets
    keepMask = sigma_clip_mask(values, maxIters=1)
    keptDf = groupDf[keepMask]
    kept = values[keepMask]
    exposureMedians = (
        keptDf.groupby("Filename")[valueColumn].median().to_numpy(dtype=float)
    )
    exposureStd = _std_ddof1(exposureMedians)
    internalStd, internalMadStd = _internal_scatter(keptDf, valueColumn)
    nStars, starStd = _star_scatter(keptDf, valueColumn)
    return {
        "Order": order,
        "NLines": int(kept.size),
        "NExposures": int(exposureMedians.size),
        "MedianShiftAngstrom": float(np.median(kept)) if kept.size else np.nan,
        "MeanShiftAngstrom": float(np.mean(kept)) if kept.size else np.nan,
        "StdShiftAngstrom": _std_ddof1(kept),
        # LINES IN ONE EXPOSURE SHARE ITS OFFSET, SO THE ERROR COMES FROM THE PER-EXPOSURE MEDIANS
        "SEMeanAngstrom": exposureStd / np.sqrt(exposureMedians.size),
        "InternalScatterStd": internalStd,
        "InternalScatterMadStd": internalMadStd,
        "ExposureScatterStd": exposureStd,
        "NStars": nStars,
        "StarScatterStd": starStd,
        "NClipped": int(np.isfinite(values).sum() - kept.size),
    }


def _select_anchor_lines(resultsDf: pd.DataFrame, anchorByOrder: dict) -> pd.DataFrame:
    """*keep only the anchor line type for each order that has one, and tag rows with an Anchor column*"""
    anchors = resultsDf["Order"].map(anchorByOrder).fillna("all")
    keepMask = (anchors == "all") | (anchors == resultsDf["LineType"])
    return resultsDf.assign(Anchor=anchors)[keepMask]


def clipped_order_stats(
    resultsDf: pd.DataFrame,
    valueColumn: str = "ResidualAngstrom",
    anchorByOrder: dict | None = None,
) -> pd.DataFrame:
    """*sigma-clipped summary statistics of a value column for each order*

    **Key Arguments:**

    - ``resultsDf`` -- per-line DataFrame with ``Filename``, ``Order``, ``LineType`` and ``valueColumn``
    - ``valueColumn`` -- column to summarise. Default *ResidualAngstrom*
    - ``anchorByOrder`` -- optional dict mapping order to the line type that anchors it, e.g. ``{"u": "stellar", "r": "telluric"}``. Only that line type is used for the order and it is recorded in an ``Anchor`` column. Orders missing from the dict use all line types and get the anchor *all*.

    **Return:**

    - ``statsDf`` -- one row per order (u, g, r, i, then any others) with ``NLines``, ``NExposures``, ``MedianShiftAngstrom``, ``MeanShiftAngstrom``, ``StdShiftAngstrom``, ``SEMeanAngstrom``, ``InternalScatterStd``, ``InternalScatterMadStd`` (scatter about each exposure's own median, corrected for the fitted medians; NaN when every exposure has one line), ``ExposureScatterStd`` (scatter of the per-exposure medians), ``NStars`` and ``StarScatterStd`` (scatter of the per-star medians; NaN without a ``Star`` column) and ``NClipped``. ``SEMeanAngstrom`` is ``ExposureScatterStd / sqrt(NExposures)``. Shifts are measured minus reference, so the correction to apply to an order is minus its median (or mean) shift.

    **Usage:**

    ```
    statsDf = clipped_order_stats(residualsDf, anchorByOrder={"u": "stellar", "r": "telluric"})
    ```
    """
    import pandas as pd

    columns = list(STATS_COLUMNS)
    if anchorByOrder is not None:
        columns.insert(1, "Anchor")
        selectedDf = _select_anchor_lines(resultsDf, anchorByOrder)
    else:
        selectedDf = resultsDf
    rows = []
    for order, groupDf in selectedDf.groupby("Order"):
        row = _order_stats_row(order, groupDf, valueColumn)
        if anchorByOrder is not None:
            row["Anchor"] = groupDf["Anchor"].iloc[0]
        rows.append(row)
    statsDf = pd.DataFrame(rows, columns=columns)
    rank = statsDf["Order"].map(
        lambda order: (
            ORDER_SEQUENCE.index(order)
            if order in ORDER_SEQUENCE
            else len(ORDER_SEQUENCE)
        )
    )
    return (
        statsDf.assign(OrderRank=rank)
        .sort_values(["OrderRank", "Order"], kind="stable")
        .drop(columns="OrderRank")
        .reset_index(drop=True)
    )


def barycentric_correction_kms(
    header: Mapping[str, Any], location: EarthLocation | None = None
) -> float:
    """*compute the barycentric radial-velocity correction for one exposure*

    **Key Arguments:**

    - ``header`` -- FITS primary header (or dict-like) with ``RA`` and ``DEC`` in degrees and ``MJD-OBS``
    - ``location`` -- astropy ``EarthLocation`` of the telescope. Default *None* (the NTT at La Silla)

    **Return:**

    - ``correctionKms`` -- the velocity (km/s) to add to a topocentric velocity to get a barycentric velocity. ``MJD-OBS`` is the exposure start; for exposures under an hour the error is below 0.05 km/s

    **Usage:**

    ```
    correctionKms = barycentric_correction_kms(fits.getheader(path))
    ```
    """
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time

    missingKeys = [key for key in ("RA", "DEC", "MJD-OBS") if key not in header]
    if missingKeys:
        raise ValueError(
            f"header is missing the keyword(s) needed for the barycentric correction: {', '.join(missingKeys)}"
        )

    if location is None:
        location = EarthLocation.from_geodetic(
            lon=NTT_LONGITUDE_DEG * u.deg,
            lat=NTT_LATITUDE_DEG * u.deg,
            height=NTT_HEIGHT_M * u.m,
        )

    targetCoord = SkyCoord(
        ra=float(header["RA"]) * u.deg, dec=float(header["DEC"]) * u.deg
    )
    obsTime = Time(float(header["MJD-OBS"]), format="mjd")
    correction = targetCoord.radial_velocity_correction(
        "barycentric", obstime=obsTime, location=location
    )
    return float(correction.to(u.km / u.s).value)
