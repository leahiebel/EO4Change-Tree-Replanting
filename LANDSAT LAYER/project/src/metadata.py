"""Processing metadata generation."""

from __future__ import annotations

import datetime as dt
from collections.abc import Mapping
from pathlib import Path

import ee
import geopandas as gpd

from .landsat_collection import LANDSAT_8_COLLECTION, LANDSAT_9_COLLECTION
from .quality_mask import ST_CDIST_SCALE_KM, ST_QA_SCALE_K
from .temperature import ST_B10_OFFSET_K, ST_B10_SCALE
from .web_layers import write_json

LIMITATIONS = [
    "LST is surface temperature, not 2-m air temperature.",
    "The output is distributed on a 30-m grid.",
    "The native TIRS thermal sampling is approximately 100 m.",
    "Resampling to 30 m does not create independent 30-m thermal information.",
    "Cloud proximity and emissivity uncertainty can affect LST.",
    "Young or narrow plantations may be mixed with surrounding land cover.",
    "This is a surface-temperature anomaly indicator, not a direct measurement of plant physiological stress.",
]


def _dates_from_collection(collection: ee.ImageCollection) -> list[str]:
    millis = collection.aggregate_array("system:time_start").getInfo()
    return [dt.datetime.utcfromtimestamp(ms / 1000).date().isoformat() for ms in millis]


def aoi_area_km2(aoi_gdf: gpd.GeoDataFrame) -> float:
    """Calculate AOI area in square kilometres using an equal-area projection."""
    return float(aoi_gdf.to_crs("EPSG:3035").area.sum() / 1_000_000)


def valid_observation_percentage(
    valid_observation_image: ee.Image,
    aoi: ee.Geometry,
    scale: int,
    crs: str,
) -> float | None:
    """Calculate percentage of AOI pixels with at least one valid observation."""
    valid = valid_observation_image.gt(0).rename("valid")
    stats = valid.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi,
        scale=scale,
        crs=crs,
        bestEffort=True,
        tileScale=4,
        maxPixels=1_000_000_000,
    ).getInfo()
    value = stats.get("valid")
    return None if value is None else float(value) * 100.0


def build_processing_metadata(
    *,
    config: Mapping[str, object],
    aoi_gdf: gpd.GeoDataFrame,
    aoi: ee.Geometry,
    landsat8_raw: ee.ImageCollection,
    landsat9_raw: ee.ImageCollection,
    valid_collection: ee.ImageCollection,
    valid_observation_image: ee.Image,
    export_tasks: list[Mapping[str, object]],
) -> dict[str, object]:
    """Build metadata JSON using only small Earth Engine summaries."""
    l8_count = int(landsat8_raw.size().getInfo())
    l9_count = int(landsat9_raw.size().getInfo())
    valid_count = int(valid_collection.size().getInfo())
    actual_dates = _dates_from_collection(valid_collection)

    return {
        "processing_date_utc": dt.datetime.now(dt.UTC).isoformat(),
        "earth_engine_collections": [LANDSAT_8_COLLECTION, LANDSAT_9_COLLECTION],
        "landsat_satellites_used": ["Landsat 8", "Landsat 9"],
        "requested_date_range": {"start": config["start_date"], "end": config["end_date"]},
        "actual_acquisition_dates_used": actual_dates,
        "seasonal_window": {
            "start_month": config["season_start_month"],
            "end_month": config["season_end_month"],
        },
        "baseline_period": {
            "start_year": config["baseline_start_year"],
            "end_year": config["baseline_end_year"],
        },
        "target_year": config["target_year"],
        "aoi_area_km2": aoi_area_km2(aoi_gdf),
        "export_crs": config["output_crs"],
        "export_scale_m": config["export_scale_m"],
        "st_b10_scale_offset": {"scale": ST_B10_SCALE, "offset_kelvin": ST_B10_OFFSET_K},
        "qa_mask_rules": {
            "qa_pixel_excluded": [
                "fill",
                "dilated_cloud",
                "high_confidence_cirrus",
                "cloud",
                "cloud_shadow",
                "snow_ice",
                "water_when_configured",
            ],
            "qa_radsat": "exclude all saturated or terrain-occluded pixels",
            "st_cdist_scale_km": ST_CDIST_SCALE_KM,
            "st_qa_scale_k": ST_QA_SCALE_K,
        },
        "cloud_distance_threshold_km": config["min_cloud_distance_km"],
        "uncertainty_threshold_k": config["max_st_uncertainty_k"],
        "number_landsat8_images": l8_count,
        "number_landsat9_images": l9_count,
        "number_valid_images": valid_count,
        "percentage_aoi_with_valid_observations": valid_observation_percentage(
            valid_observation_image,
            aoi,
            int(config["export_scale_m"]),
            str(config["output_crs"]),
        ),
        "export_tasks": export_tasks,
        "known_limitations": LIMITATIONS,
    }


def write_processing_metadata(path: str | Path, metadata: Mapping[str, object]) -> None:
    write_json(path, metadata)

