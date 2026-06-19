"""Temporal LST products built server-side in Earth Engine."""

from __future__ import annotations

from collections.abc import Mapping

import ee


def filter_season(collection: ee.ImageCollection, start_month: int, end_month: int) -> ee.ImageCollection:
    """Filter an image collection to an inclusive month range."""
    return collection.filter(ee.Filter.calendarRange(start_month, end_month, "month"))


def filter_year(collection: ee.ImageCollection, year: int) -> ee.ImageCollection:
    """Filter an image collection to a calendar year."""
    return collection.filter(ee.Filter.calendarRange(year, year, "year"))


def target_season_collection(collection: ee.ImageCollection, config: Mapping[str, object]) -> ee.ImageCollection:
    return filter_season(
        filter_year(collection, int(config["target_year"])),
        int(config["season_start_month"]),
        int(config["season_end_month"]),
    )


def baseline_season_collection(collection: ee.ImageCollection, config: Mapping[str, object]) -> ee.ImageCollection:
    return filter_season(
        collection.filter(
            ee.Filter.calendarRange(
                int(config["baseline_start_year"]), int(config["baseline_end_year"]), "year"
            )
        ),
        int(config["season_start_month"]),
        int(config["season_end_month"]),
    )


def build_target_composite(collection: ee.ImageCollection, config: Mapping[str, object]) -> dict[str, ee.Image]:
    """Build target-year seasonal LST products."""
    target = target_season_collection(collection, config)
    lst = target.select("LST_C")
    return {
        "01_lst_median_target_C.tif": lst.median().rename("lst_median_target_C").toFloat(),
        "02_lst_mean_target_C.tif": lst.mean().rename("lst_mean_target_C").toFloat(),
        "03_lst_p10_target_C.tif": lst.reduce(ee.Reducer.percentile([10])).rename("lst_p10_target_C").toFloat(),
        "04_lst_p90_target_C.tif": lst.reduce(ee.Reducer.percentile([90])).rename("lst_p90_target_C").toFloat(),
        "05_lst_std_target_C.tif": lst.reduce(ee.Reducer.stdDev()).rename("lst_std_target_C").toFloat(),
        "06_lst_valid_observations.tif": target.select("valid_observation")
        .sum()
        .rename("lst_valid_observations")
        .toUint16(),
        "07_lst_uncertainty_K.tif": target.select("ST_uncertainty_K")
        .median()
        .rename("lst_uncertainty_K")
        .toFloat(),
    }


def build_baseline_products(collection: ee.ImageCollection, config: Mapping[str, object]) -> dict[str, ee.Image]:
    """Build baseline seasonal climatology products."""
    baseline = baseline_season_collection(collection, config)
    lst = baseline.select("LST_C")
    return {
        "08_lst_baseline_median_C.tif": lst.median().rename("lst_baseline_median_C").toFloat(),
        "baseline_mean_C": lst.mean().rename("lst_baseline_mean_C").toFloat(),
        "baseline_std_C": lst.reduce(ee.Reducer.stdDev()).rename("lst_baseline_std_C").toFloat(),
        "baseline_valid_observations": baseline.select("valid_observation")
        .sum()
        .rename("lst_baseline_valid_observations")
        .toUint16(),
    }


def build_anomaly_products(
    target_products: dict[str, ee.Image],
    baseline_products: dict[str, ee.Image],
    *,
    moderate_z: float,
    strong_z: float,
) -> dict[str, ee.Image]:
    """Build thermal anomaly, z-score and categorical stress anomaly class."""
    target_median = target_products["01_lst_median_target_C.tif"]
    baseline_median = baseline_products["08_lst_baseline_median_C.tif"]
    baseline_mean = baseline_products["baseline_mean_C"]
    baseline_std = baseline_products["baseline_std_C"]
    nonzero_std = baseline_std.neq(0)

    anomaly = target_median.subtract(baseline_median).rename("lst_anomaly_C").toFloat()
    zscore = (
        target_median.subtract(baseline_mean)
        .divide(baseline_std)
        .updateMask(nonzero_std)
        .rename("lst_zscore")
        .toFloat()
    )

    valid = zscore.mask()
    stress = (
        ee.Image.constant(0)
        .where(valid, 1)
        .where(zscore.gte(moderate_z), 2)
        .where(zscore.gte(strong_z), 3)
        .rename("thermal_stress_anomaly_class")
        .toUint8()
    )
    return {
        "09_lst_anomaly_C.tif": anomaly,
        "10_lst_zscore.tif": zscore,
        "11_thermal_stress_class.tif": stress,
    }


def build_all_raster_products(collection: ee.ImageCollection, config: Mapping[str, object]) -> dict[str, ee.Image]:
    """Build all required raster products keyed by output filename."""
    target = build_target_composite(collection, config)
    baseline = build_baseline_products(collection, config)
    anomaly = build_anomaly_products(
        target,
        baseline,
        moderate_z=float(config["zscore_moderate_threshold"]),
        strong_z=float(config["zscore_strong_threshold"]),
    )
    return {
        **target,
        "08_lst_baseline_median_C.tif": baseline["08_lst_baseline_median_C.tif"],
        **anomaly,
    }

