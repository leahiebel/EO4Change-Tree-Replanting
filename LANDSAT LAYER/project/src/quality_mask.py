"""Reusable Landsat LST quality masking."""

from __future__ import annotations

import ee

ST_QA_SCALE_K = 0.01
ST_CDIST_SCALE_KM = 0.01

QA_FILL_BIT = 0
QA_DILATED_CLOUD_BIT = 1
QA_CIRRUS_BIT = 2
QA_CLOUD_BIT = 3
QA_CLOUD_SHADOW_BIT = 4
QA_SNOW_BIT = 5
QA_WATER_BIT = 7


def _bit_is_zero(image: ee.Image, bit: int) -> ee.Image:
    return image.bitwiseAnd(1 << bit).eq(0)


def qa_pixel_is_clear(value: int, mask_water: bool = True) -> bool:
    """Pure-Python QA_PIXEL bit check used by unit tests."""
    blocked_bits = [
        QA_FILL_BIT,
        QA_DILATED_CLOUD_BIT,
        QA_CIRRUS_BIT,
        QA_CLOUD_BIT,
        QA_CLOUD_SHADOW_BIT,
        QA_SNOW_BIT,
    ]
    if mask_water:
        blocked_bits.append(QA_WATER_BIT)
    return all((value & (1 << bit)) == 0 for bit in blocked_bits)


def qa_radsat_is_clear(value: int) -> bool:
    """Return true when QA_RADSAT has no saturation or terrain-occlusion flag."""
    return value == 0


def mask_landsat_lst(
    image: ee.Image,
    *,
    mask_water: bool,
    min_cloud_distance_km: float,
    max_st_uncertainty_k: float,
) -> ee.Image:
    """Mask Landsat Collection 2 L2SP pixels unsuitable for LST analysis.

    Mask criteria:
    - QA_PIXEL excludes fill, dilated cloud, high-confidence cirrus, cloud,
      cloud shadow, snow/ice and optionally water.
    - QA_RADSAT excludes radiometrically saturated and terrain-occluded pixels.
    - ST_CDIST excludes pixels too close to clouds after applying the 0.01 km
      scale factor.
    - ST_QA excludes pixels whose surface-temperature uncertainty exceeds the
      configured Kelvin threshold after applying the 0.01 K scale factor.
    """
    qa_pixel = image.select("QA_PIXEL")
    qa_radsat = image.select("QA_RADSAT")
    uncertainty = image.select("ST_QA").multiply(ST_QA_SCALE_K).rename("ST_uncertainty_K")
    cloud_distance = image.select("ST_CDIST").multiply(ST_CDIST_SCALE_KM).rename(
        "cloud_distance_km"
    )

    mask = (
        _bit_is_zero(qa_pixel, QA_FILL_BIT)
        .And(_bit_is_zero(qa_pixel, QA_DILATED_CLOUD_BIT))
        .And(_bit_is_zero(qa_pixel, QA_CIRRUS_BIT))
        .And(_bit_is_zero(qa_pixel, QA_CLOUD_BIT))
        .And(_bit_is_zero(qa_pixel, QA_CLOUD_SHADOW_BIT))
        .And(_bit_is_zero(qa_pixel, QA_SNOW_BIT))
        .And(qa_radsat.eq(0))
        .And(cloud_distance.gte(min_cloud_distance_km))
        .And(uncertainty.lte(max_st_uncertainty_k))
    )
    if mask_water:
        mask = mask.And(_bit_is_zero(qa_pixel, QA_WATER_BIT))

    valid = ee.Image.constant(1).rename("valid_observation").uint8().updateMask(mask)
    return image.addBands([uncertainty.toFloat(), cloud_distance.toFloat(), valid], overwrite=True).updateMask(mask)

