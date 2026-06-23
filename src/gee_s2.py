"""
EO4Change Group 4 — Sentinel-2 collection builder, spectral indices, composites.
"""
import ee
from gee_config import (aoi, CS_BAND, CS_THRESH, MAX_CLOUD,
                         SPECTRAL_BANDS, BIOPHYS_VARS, RGB_BANDS, ALL_BANDS, COMPOSITE_BANDS)


def build_s2_col(ws: str, we: str) -> ee.ImageCollection:
    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(aoi)
            .filterDate(ws, we)
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", MAX_CLOUD)))
    cs = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
    return (s2.linkCollection(cs, [CS_BAND])
              .map(lambda img: img.updateMask(img.select(CS_BAND).gte(CS_THRESH))))


def add_spectral_indices(img: ee.Image) -> ee.Image:
    ndvi = img.normalizedDifference(["B8",  "B4" ]).rename("NDVI")
    ndre = img.normalizedDifference(["B8A", "B5" ]).rename("NDRE")
    ndmi = img.normalizedDifference(["B8",  "B11"]).rename("NDMI")
    ndwi = img.normalizedDifference(["B3",  "B8" ]).rename("NDWI")
    nbr  = img.normalizedDifference(["B8",  "B12"]).rename("NBR")
    bsi  = (img.select("B11").add(img.select("B4"))
               .subtract(img.select("B8").add(img.select("B2")))
               .divide(img.select("B11").add(img.select("B4"))
                          .add(img.select("B8")).add(img.select("B2")))
               ).rename("BSI")
    return img.addBands([ndvi, ndre, ndmi, bsi, ndwi, nbr])


def add_biophys_proxies(img: ee.Image) -> ee.Image:
    ndvi  = img.select("NDVI")
    laie  = ndvi.multiply(3.618).subtract(0.118).clamp(0, 8).rename("laie")
    fapar = ndvi.multiply(1.136).subtract(0.040).clamp(0, 1).rename("fapar")
    fcover = laie.multiply(-0.5).exp().multiply(-1).add(1).clamp(0, 1).rename("fcover")
    return img.addBands([laie, fapar, fcover])


def add_requested_biophys(img: ee.Image) -> ee.Image:
    return add_biophys_proxies(img)


_PLACEHOLDER = (ee.Image.constant([0] * len(COMPOSITE_BANDS))
                  .rename(COMPOSITE_BANDS)
                  .updateMask(ee.Image(0)))


def make_composite(ws: str, we: str) -> ee.Image:
    col = (build_s2_col(ws, we)
             .map(add_spectral_indices)
             .map(add_requested_biophys)
             .map(lambda img: img.select(COMPOSITE_BANDS)))
    return ee.Image(
        ee.Algorithms.If(
            col.size().gt(0),
            col.median().select(COMPOSITE_BANDS).clip(aoi),
            _PLACEHOLDER
        )
    )


def make_seasonal_composite(year_start: int, year_end: int,
                            month_start: int, month_end: int) -> ee.Image:
    cols = []
    for y in range(year_start, year_end + 1):
        ws = f"{y}-{month_start:02d}-01"
        we = f"{y}-{month_end:02d}-01"
        cols.append(build_s2_col(ws, we)
                      .map(add_spectral_indices)
                      .map(add_requested_biophys)
                      .map(lambda img: img.select(COMPOSITE_BANDS)))
    merged = cols[0]
    for c in cols[1:]:
        merged = merged.merge(c)
    return merged.median().select(COMPOSITE_BANDS).clip(aoi)
