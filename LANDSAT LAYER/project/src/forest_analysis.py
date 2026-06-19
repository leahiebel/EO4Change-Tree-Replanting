"""Forest-control LST statistics and yearly time-series tables."""

from __future__ import annotations

from collections.abc import Mapping

import ee

from .composites import filter_season

STAT_BANDS = ["LST_C", "ST_uncertainty_K", "valid_observation"]


def _statistics_reducer() -> ee.Reducer:
    return (
        ee.Reducer.median()
        .combine(ee.Reducer.mean(), sharedInputs=True)
        .combine(ee.Reducer.stdDev(), sharedInputs=True)
        .combine(ee.Reducer.count(), sharedInputs=True)
    )


def _prefix_properties(feature: ee.Feature, prefix: str, preserve: list[str]) -> ee.Feature:
    keys = feature.propertyNames()
    rename_from = keys.filter(ee.Filter.inList("item", preserve).Not())
    rename_to = rename_from.map(lambda key: ee.String(prefix).cat(ee.String(key)))
    return feature.select(rename_from, rename_to).copyProperties(feature, preserve)


def forest_control_statistics(
    target_collection: ee.ImageCollection,
    forest_fc: ee.FeatureCollection,
    control_fc: ee.FeatureCollection | None,
    config: Mapping[str, object],
) -> ee.FeatureCollection:
    """Calculate forest and optional paired control LST statistics.

    Pairing is based on a shared ``feature_id`` property. When controls are
    missing this returns forest-only statistics and labels them exploratory.
    """
    scale = int(config["export_scale_m"])
    crs = str(config["output_crs"])
    composite = target_collection.select(STAT_BANDS).median()
    reducer = _statistics_reducer()

    forest_stats = composite.reduceRegions(
        collection=forest_fc,
        reducer=reducer,
        scale=scale,
        crs=crs,
        tileScale=4,
    ).map(lambda f: ee.Feature(f).set("analysis_role", "forest"))

    if control_fc is None:
        return forest_stats.map(
            lambda f: ee.Feature(f).set(
                "control_status",
                "missing_control_exploratory_no_cooling_difference",
                "forest_cooling_C",
                None,
            )
        )

    control_stats = composite.reduceRegions(
        collection=control_fc,
        reducer=reducer,
        scale=scale,
        crs=crs,
        tileScale=4,
    ).map(lambda f: _prefix_properties(ee.Feature(f).set("analysis_role", "control"), "control_", ["feature_id"]))

    join = ee.Join.inner()
    paired = join.apply(
        forest_stats,
        control_stats,
        ee.Filter.equals(leftField="feature_id", rightField="feature_id"),
    )

    def _cooling(joined: ee.Feature) -> ee.Feature:
        forest = ee.Feature(joined.get("primary"))
        control = ee.Feature(joined.get("secondary"))
        cooling = ee.Number(forest.get("LST_C_median")).subtract(
            ee.Number(control.get("control_LST_C_median"))
        )
        return forest.copyProperties(control).set(
            "forest_cooling_C",
            cooling,
            "interpretation_note",
            "negative means forest polygon is cooler than paired control; not causal unless controls are comparable",
        )

    return ee.FeatureCollection(paired.map(_cooling))


def yearly_lst_timeseries(
    collection: ee.ImageCollection,
    aoi: ee.Geometry,
    forest_fc: ee.FeatureCollection,
    control_fc: ee.FeatureCollection | None,
    config: Mapping[str, object],
) -> ee.FeatureCollection:
    """Build a tidy yearly seasonal median LST table for AOI, forest and controls."""
    years = ee.List.sequence(int(config["baseline_start_year"]), int(config["target_year"]))
    scale = int(config["export_scale_m"])
    crs = str(config["output_crs"])

    def _stats_for_regions(year: ee.Number, regions: ee.FeatureCollection, role: str) -> ee.FeatureCollection:
        seasonal = filter_season(
            collection.filter(ee.Filter.calendarRange(year, year, "year")),
            int(config["season_start_month"]),
            int(config["season_end_month"]),
        )
        image = seasonal.select("LST_C").median().rename("seasonal_median_lst_C")
        return image.reduceRegions(
            collection=regions,
            reducer=ee.Reducer.median().combine(ee.Reducer.count(), sharedInputs=True),
            scale=scale,
            crs=crs,
            tileScale=4,
        ).map(lambda f: ee.Feature(f).set("year", year, "analysis_role", role))

    aoi_fc = ee.FeatureCollection([ee.Feature(aoi, {"feature_id": "aoi", "site_name": "complete_aoi"})])

    def _for_year(year: ee.Number) -> ee.FeatureCollection:
        parts = _stats_for_regions(year, aoi_fc, "aoi").merge(
            _stats_for_regions(year, forest_fc, "forest")
        )
        if control_fc is not None:
            parts = parts.merge(_stats_for_regions(year, control_fc, "control"))
        return parts

    return ee.FeatureCollection(years.map(_for_year)).flatten()


def difference_in_differences(
    collection: ee.ImageCollection,
    forest_fc: ee.FeatureCollection,
    control_fc: ee.FeatureCollection,
    config: Mapping[str, object],
) -> ee.FeatureCollection:
    """Estimate optional pre/post planting difference-in-differences where possible."""
    scale = int(config["export_scale_m"])
    crs = str(config["output_crs"])

    def _median_for(feature: ee.Feature, start_year: ee.Number, end_year: ee.Number) -> ee.Number:
        seasonal = filter_season(
            collection.filter(ee.Filter.calendarRange(start_year, end_year, "year")),
            int(config["season_start_month"]),
            int(config["season_end_month"]),
        )
        value = seasonal.select("LST_C").median().reduceRegion(
            reducer=ee.Reducer.median(),
            geometry=feature.geometry(),
            scale=scale,
            crs=crs,
            bestEffort=True,
            tileScale=4,
        ).get("LST_C")
        return ee.Number(value)

    def _per_forest(forest: ee.Feature) -> ee.Feature:
        planting_year = ee.Number(forest.get("planting_year"))
        control = ee.Feature(control_fc.filter(ee.Filter.eq("feature_id", forest.get("feature_id"))).first())
        pre_start = ee.Number(config["baseline_start_year"])
        pre_end = planting_year.subtract(1)
        post_start = planting_year
        post_end = ee.Number(config["target_year"])

        pre_forest = _median_for(forest, pre_start, pre_end)
        post_forest = _median_for(forest, post_start, post_end)
        pre_control = _median_for(control, pre_start, pre_end)
        post_control = _median_for(control, post_start, post_end)
        effect = post_forest.subtract(pre_forest).subtract(post_control.subtract(pre_control))

        return ee.Feature(forest.geometry(), {
            "feature_id": forest.get("feature_id"),
            "planting_year": planting_year,
            "pre_forest_lst_C": pre_forest,
            "post_forest_lst_C": post_forest,
            "pre_control_lst_C": pre_control,
            "post_control_lst_C": post_control,
            "effect_C": effect,
            "interpretation_note": "exploratory DiD; requires comparable controls and stable observation support",
        })

    return ee.FeatureCollection(forest_fc.filter(ee.Filter.notNull(["planting_year"])).map(_per_forest))

