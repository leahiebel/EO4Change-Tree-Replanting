"""Landsat 8/9 Collection 2 L2SP collection construction."""

from __future__ import annotations

from collections.abc import Mapping

import ee

from .quality_mask import mask_landsat_lst
from .temperature import add_lst_celsius

LANDSAT_8_COLLECTION = "LANDSAT/LC08/C02/T1_L2"
LANDSAT_9_COLLECTION = "LANDSAT/LC09/C02/T1_L2"

REQUIRED_BANDS = ["ST_B10", "ST_QA", "ST_CDIST", "QA_PIXEL", "QA_RADSAT"]
PRESERVED_PROPERTIES = [
    "system:time_start",
    "SPACECRAFT_ID",
    "LANDSAT_PRODUCT_ID",
    "PROCESSING_LEVEL",
    "CLOUD_COVER",
    "WRS_PATH",
    "WRS_ROW",
]


def _filtered_collection(collection_id: str, aoi: ee.Geometry, config: Mapping[str, object]) -> ee.ImageCollection:
    return (
        ee.ImageCollection(collection_id)
        .filterBounds(aoi)
        .filterDate(str(config["start_date"]), str(config["end_date"]))
        .filter(ee.Filter.lte("CLOUD_COVER", float(config["max_scene_cloud_cover"])))
        .filter(ee.Filter.eq("PROCESSING_LEVEL", "L2SP"))
        .select(REQUIRED_BANDS, REQUIRED_BANDS)
    )


def get_landsat_l2sp_collections(
    aoi: ee.Geometry, config: Mapping[str, object]
) -> tuple[ee.ImageCollection, ee.ImageCollection]:
    """Return filtered Landsat 8 and Landsat 9 L2SP collections."""
    return (
        _filtered_collection(LANDSAT_8_COLLECTION, aoi, config),
        _filtered_collection(LANDSAT_9_COLLECTION, aoi, config),
    )


def preprocess_collection(collection: ee.ImageCollection, config: Mapping[str, object]) -> ee.ImageCollection:
    """Apply identical masking and temperature conversion to a collection."""

    def _preprocess(image: ee.Image) -> ee.Image:
        masked = mask_landsat_lst(
            image,
            mask_water=bool(config["mask_water"]),
            min_cloud_distance_km=float(config["min_cloud_distance_km"]),
            max_st_uncertainty_k=float(config["max_st_uncertainty_k"]),
        )
        return add_lst_celsius(masked).copyProperties(image, PRESERVED_PROPERTIES)

    return collection.map(_preprocess)


def build_merged_landsat_lst_collection(
    aoi: ee.Geometry, config: Mapping[str, object]
) -> tuple[ee.ImageCollection, ee.ImageCollection, ee.ImageCollection]:
    """Filter, preprocess and merge Landsat 8 and Landsat 9 for LST."""
    landsat8, landsat9 = get_landsat_l2sp_collections(aoi, config)
    merged = preprocess_collection(landsat8, config).merge(preprocess_collection(landsat9, config))
    return landsat8, landsat9, merged.sort("system:time_start")

