"""Landsat 8/9 Land Surface Temperature layer for the EO4Change pipeline.

Distilled from the LANDSAT LAYER sub-project:
  quality_mask.py + temperature.py + landsat_collection.py + composites.py

Public API: build_lst_layer(aoi, cfg) -> dict[str, ee.Image]
"""
from __future__ import annotations

from collections.abc import Mapping

import ee

# --- ST_B10 conversion constants (USGS Collection 2 L2 product guide) --------
_ST_SCALE  = 0.00341802
_ST_OFFSET = 149.0          # Kelvin offset
_CELSIUS   = 273.15

_L8 = "LANDSAT/LC08/C02/T1_L2"
_L9 = "LANDSAT/LC09/C02/T1_L2"


# --- QA masking --------------------------------------------------------------

def _mask_landsat(
    image: ee.Image,
    *,
    mask_water: bool,
    min_cloud_km: float,
    max_uncert_k: float,
) -> ee.Image:
    """Mask fill, cloud, shadow, snow, radsat, cloud-proximity and high-uncertainty pixels."""
    qa = image.select("QA_PIXEL")
    bad = (
        qa.bitwiseAnd(1 << 0).neq(0)           # fill
        .Or(qa.bitwiseAnd(1 << 1).neq(0))      # dilated cloud
        .Or(qa.bitwiseAnd(1 << 2).neq(0))      # cirrus
        .Or(qa.bitwiseAnd(1 << 3).neq(0))      # cloud
        .Or(qa.bitwiseAnd(1 << 4).neq(0))      # cloud shadow
        .Or(qa.bitwiseAnd(1 << 5).neq(0))      # snow/ice
        .Or(image.select("QA_RADSAT").neq(0))  # radiometric saturation
    )
    if mask_water:
        bad = bad.Or(qa.bitwiseAnd(1 << 7).neq(0))
    # ST_CDIST scale = 0.01 km/DN; ST_QA scale = 0.01 K/DN
    bad = (bad
           .Or(image.select("ST_CDIST").multiply(0.01).lt(min_cloud_km))
           .Or(image.select("ST_QA").multiply(0.01).gt(max_uncert_k)))
    return image.updateMask(bad.Not())


# --- Temperature conversion --------------------------------------------------

def _add_lst_celsius(image: ee.Image) -> ee.Image:
    """Add LST_C band (degrees Celsius) from ST_B10 digital number."""
    lst = (
        image.select("ST_B10")
        .multiply(_ST_SCALE)
        .add(_ST_OFFSET - _CELSIUS)
        .rename("LST_C")
        .toFloat()
    )
    return image.addBands(lst, overwrite=True)


# --- Collection builder ------------------------------------------------------

def _build_collection(collection_id: str, aoi: ee.Geometry, cfg: Mapping) -> ee.ImageCollection:
    return (
        ee.ImageCollection(collection_id)
        .filterBounds(aoi)
        .filterDate(str(cfg["start_date"]), str(cfg["end_date"]))
        .filter(ee.Filter.lte("CLOUD_COVER", float(cfg["max_scene_cloud_cover"])))
        .filter(ee.Filter.eq("PROCESSING_LEVEL", "L2SP"))
        .select(["ST_B10", "ST_QA", "ST_CDIST", "QA_PIXEL", "QA_RADSAT"])
    )


# --- Public API --------------------------------------------------------------

def build_lst_layer(aoi: ee.Geometry, cfg: Mapping) -> dict[str, ee.Image]:
    """Build Landsat LST raster products ready to add as folium map layers.

    Parameters
    ----------
    aoi : ee.Geometry
        Study area (same geometry used by the S2 pipeline).
    cfg : Mapping
        The ``landsat:`` block from the pipeline config YAML.

    Returns
    -------
    dict
        ``lst_median``  — target-season median LST in °C (float)
        ``lst_anomaly`` — target median minus baseline median in °C (float)
        ``lst_stress``  — thermal stress class 0–3 (uint8):
                          0 = no valid data, 1 = near-normal,
                          2 = moderate warm anomaly, 3 = strong warm anomaly
    """
    mask_kw = dict(
        mask_water=bool(cfg.get("mask_water", True)),
        min_cloud_km=float(cfg.get("min_cloud_distance_km", 1.0)),
        max_uncert_k=float(cfg.get("max_st_uncertainty_k", 3.0)),
    )

    def _preprocess(col_id: str) -> ee.ImageCollection:
        return _build_collection(col_id, aoi, cfg).map(
            lambda img: _add_lst_celsius(_mask_landsat(img, **mask_kw))
        )

    collection = _preprocess(_L8).merge(_preprocess(_L9)).sort("system:time_start")

    sm = int(cfg["season_start_month"])
    em = int(cfg["season_end_month"])

    target = (
        collection
        .filter(ee.Filter.calendarRange(int(cfg["target_year"]), int(cfg["target_year"]), "year"))
        .filter(ee.Filter.calendarRange(sm, em, "month"))
    )
    baseline = (
        collection
        .filter(ee.Filter.calendarRange(
            int(cfg["baseline_start_year"]), int(cfg["baseline_end_year"]), "year"
        ))
        .filter(ee.Filter.calendarRange(sm, em, "month"))
    )

    lst_target   = target.select("LST_C").median().rename("lst_median_target_C").clip(aoi)
    lst_baseline = baseline.select("LST_C").median()
    lst_anomaly  = lst_target.subtract(lst_baseline).rename("lst_anomaly_C").clip(aoi)

    baseline_mean = baseline.select("LST_C").mean()
    baseline_std  = baseline.select("LST_C").reduce(ee.Reducer.stdDev()).rename("LST_C")
    zscore = (
        lst_target.subtract(baseline_mean)
        .divide(baseline_std)
        .updateMask(baseline_std.neq(0))
        .rename("lst_zscore")
    )

    mz = float(cfg.get("zscore_moderate_threshold", 1.5))
    sz = float(cfg.get("zscore_strong_threshold", 2.5))
    stress = (
        ee.Image.constant(0)
        .where(zscore.mask(), 1)
        .where(zscore.gte(mz), 2)
        .where(zscore.gte(sz), 3)
        .rename("thermal_stress_class")
        .toUint8()
        .clip(aoi)
    )

    return {
        "lst_median":  lst_target,
        "lst_anomaly": lst_anomaly,
        "lst_stress":  stress,
    }
