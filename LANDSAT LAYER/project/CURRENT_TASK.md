# Current Task: Landsat LST With Earth Engine Python API

## Decision

For this project, the simplest operational path is to use the Google Earth Engine Python API. Python will access Landsat datasets directly in the cloud, filter the Portuguese domain, convert surface temperature, and produce the final raster without downloading all original Landsat scenes locally.

This current task is Landsat-only.

## Dataset

Use Landsat 8 and Landsat 9 Collection 2 Tier 1 Level-2:

```python
LANDSAT/LC08/C02/T1_L2
LANDSAT/LC09/C02/T1_L2
```

Meaning:

- `LC08` / `LC09`: Landsat 8 and Landsat 9;
- `C02`: Collection 2;
- `T1`: Tier 1 scenes with geometric quality suitable for time series;
- `L2`: Level-2, containing Surface Reflectance and, for `L2SP` products, Surface Temperature.

The thermal surface-temperature band is:

```text
ST_B10
```

`ST_B10` is not raw radiance. It is a Level-2 Land Surface Temperature product stored in Earth Engine as a scaled integer value.

Conversion:

```text
LST [K] = ST_B10 * 0.00341802 + 149.0
LST [degC] = ST_B10 * 0.00341802 + 149.0 - 273.15
```

USGS requires applying the scale factor to convert the stored values to temperature.

The final product is distributed on a 30 m grid, but the TIRS thermal information is originally acquired at approximately 100 m and resampled to 30 m. It must not be presented as independent thermal information at every 30 m pixel.

## Installation

```bash
python -m pip install earthengine-api geopandas
```

Requirements:

1. Google account.
2. Access to Google Earth Engine.
3. Google Cloud project configured for Earth Engine.

## QGIS AOI Export

In QGIS:

```text
Right click polygon
Export
Save Features As
```

Set:

```text
Format: GeoJSON
CRS: EPSG:4326
File: portugal_domain.geojson
```

## Complete Python Example

This example:

- opens the Portuguese polygon;
- accesses Landsat 8 and 9;
- keeps only `L2SP` products with surface temperature;
- masks clouds, cloud shadows, snow, saturation and water;
- converts `ST_B10` to degrees Celsius;
- creates a summer median;
- exports a Cloud Optimized GeoTIFF;
- generates a tile URL for an HTML map prototype.

```python
from pathlib import Path

import ee
import geopandas as gpd


# ============================================================
# CONFIGURATION
# ============================================================

EARTH_ENGINE_PROJECT = "YOUR_GOOGLE_CLOUD_PROJECT_ID"

AOI_PATH = Path(
    "/Users/nicolocaron/Desktop/portugal_domain.geojson"
)

START_DATE = "2025-06-01"
END_DATE = "2025-10-01"

OUTPUT_CRS = "EPSG:3763"  # ETRS89 / Portugal TM06
EXPORT_SCALE_M = 30

DRIVE_FOLDER = "reforestation_project"
OUTPUT_NAME = "portugal_LST_summer_2025"


# ============================================================
# 1. AUTHENTICATE AND INITIALIZE EARTH ENGINE
# ============================================================

# Normally required only the first time, or when credentials expire.
ee.Authenticate()

ee.Initialize(
    project=EARTH_ENGINE_PROJECT
)


# ============================================================
# 2. READ THE PORTUGUESE DOMAIN
# ============================================================

if not AOI_PATH.exists():
    raise FileNotFoundError(
        f"AOI file not found: {AOI_PATH}"
    )

aoi_gdf = gpd.read_file(AOI_PATH)

if aoi_gdf.empty:
    raise ValueError(
        "The GeoJSON does not contain any geometry."
    )

if aoi_gdf.crs is None:
    raise ValueError(
        "The AOI does not have a defined coordinate system."
    )

# Earth Engine geometries should be provided in longitude/latitude.
aoi_gdf = aoi_gdf.to_crs("EPSG:4326")

# Combine all polygons into one geometry.
if hasattr(aoi_gdf.geometry, "union_all"):
    combined_geometry = aoi_gdf.geometry.union_all()
else:
    combined_geometry = aoi_gdf.unary_union

if combined_geometry.is_empty:
    raise ValueError(
        "The resulting AOI geometry is empty."
    )

aoi = ee.Geometry(
    combined_geometry.__geo_interface__
)


# ============================================================
# 3. LANDSAT QUALITY MASK AND LST CONVERSION
# ============================================================

def prepare_landsat_lst(image: ee.Image) -> ee.Image:
    """
    Mask invalid Landsat pixels and convert ST_B10
    from scaled digital numbers to Celsius.
    """

    qa_pixel = image.select("QA_PIXEL")

    # QA_PIXEL bits:
    # 0 = fill
    # 1 = dilated cloud
    # 2 = cirrus
    # 3 = cloud
    # 4 = cloud shadow
    # 5 = snow
    # 7 = water

    invalid_mask = (
        (1 << 0)
        | (1 << 1)
        | (1 << 2)
        | (1 << 3)
        | (1 << 4)
        | (1 << 5)
    )

    clear_sky = qa_pixel.bitwiseAnd(
        invalid_mask
    ).eq(0)

    # Exclude water for the forest analysis.
    non_water = qa_pixel.bitwiseAnd(
        1 << 7
    ).eq(0)

    # Exclude radiometrically saturated pixels.
    not_saturated = image.select(
        "QA_RADSAT"
    ).eq(0)

    valid_mask = (
        clear_sky
        .And(non_water)
        .And(not_saturated)
    )

    # Convert ST_B10 DN to Kelvin and then Celsius.
    lst_celsius = (
        image
        .select("ST_B10")
        .multiply(0.00341802)
        .add(149.0)
        .subtract(273.15)
        .rename("LST_C")
    )

    return (
        lst_celsius
        .updateMask(valid_mask)
        .copyProperties(
            image,
            [
                "system:time_start",
                "LANDSAT_PRODUCT_ID",
                "SPACECRAFT_ID",
                "PROCESSING_LEVEL",
                "CLOUD_COVER",
            ],
        )
    )


# ============================================================
# 4. LOAD LANDSAT 8 AND LANDSAT 9
# ============================================================

landsat_8 = (
    ee.ImageCollection(
        "LANDSAT/LC08/C02/T1_L2"
    )
    .filterBounds(aoi)
    .filterDate(START_DATE, END_DATE)
    .filter(
        ee.Filter.eq(
            "PROCESSING_LEVEL",
            "L2SP",
        )
    )
)

landsat_9 = (
    ee.ImageCollection(
        "LANDSAT/LC09/C02/T1_L2"
    )
    .filterBounds(aoi)
    .filterDate(START_DATE, END_DATE)
    .filter(
        ee.Filter.eq(
            "PROCESSING_LEVEL",
            "L2SP",
        )
    )
)


# ============================================================
# 5. MERGE AND PROCESS THE COLLECTIONS
# ============================================================

landsat_lst = (
    landsat_8
    .merge(landsat_9)
    .map(prepare_landsat_lst)
)

scene_count = landsat_lst.size().getInfo()

print(
    f"Valid Landsat scenes found: {scene_count}"
)

if scene_count == 0:
    raise RuntimeError(
        "No valid Landsat L2SP scenes were found "
        "for the requested domain and dates."
    )


# ============================================================
# 6. CREATE THE SEASONAL MEDIAN LST
# ============================================================

lst_median = (
    landsat_lst
    .select("LST_C")
    .median()
    .clip(aoi)
    .rename("LST_median_C")
)


# ============================================================
# 7. BASIC DOMAIN STATISTICS
# ============================================================

statistics = lst_median.reduceRegion(
    reducer=(
        ee.Reducer.mean()
        .combine(
            reducer2=ee.Reducer.median(),
            sharedInputs=True,
        )
        .combine(
            reducer2=ee.Reducer.stdDev(),
            sharedInputs=True,
        )
        .combine(
            reducer2=ee.Reducer.minMax(),
            sharedInputs=True,
        )
    ),
    geometry=aoi,
    scale=EXPORT_SCALE_M,
    maxPixels=1e13,
    tileScale=4,
)

print(
    "LST statistics [degC]:",
    statistics.getInfo(),
)


# ============================================================
# 8. EXPORT AS CLOUD OPTIMIZED GEOTIFF
# ============================================================

NODATA_VALUE = -9999

export_image = lst_median.unmask(
    value=NODATA_VALUE,
    sameFootprint=False,
)

export_task = ee.batch.Export.image.toDrive(
    image=export_image,
    description=OUTPUT_NAME,
    folder=DRIVE_FOLDER,
    fileNamePrefix=OUTPUT_NAME,
    region=aoi,
    scale=EXPORT_SCALE_M,
    crs=OUTPUT_CRS,
    maxPixels=1e13,
    fileFormat="GeoTIFF",
    formatOptions={
        "cloudOptimized": True,
        "noData": NODATA_VALUE,
    },
)

export_task.start()

print(
    "Export started."
)

print(
    f"Earth Engine task ID: {export_task.id}"
)


# ============================================================
# 9. CREATE A TILE URL FOR AN HTML MAP PROTOTYPE
# ============================================================

visualisation = {
    "min": 15,
    "max": 45,
    "palette": [
        "313695",
        "4575b4",
        "74add1",
        "abd9e9",
        "e0f3f8",
        "ffffbf",
        "fee090",
        "fdae61",
        "f46d43",
        "d73027",
        "a50026",
    ],
}

map_information = lst_median.getMapId(
    visualisation
)

tile_url = map_information[
    "tile_fetcher"
].url_format

print(
    "Earth Engine tile URL:"
)

print(tile_url)
```

## Execution Model

These lines:

```python
landsat_8 = ee.ImageCollection(
    "LANDSAT/LC08/C02/T1_L2"
)
```

do not download the dataset. They create a server-side reference to the Earth Engine collection.

Earth Engine executes work when calling:

```python
.getInfo()
```

or when starting:

```python
ee.batch.Export.image.toDrive(...)
```

For a large Portuguese domain, processing should remain server-side and only final products should be exported.

## Leaflet Prototype

The printed `tile_url` can be used in Leaflet:

```javascript
const lstLayer = L.tileLayer(
  "PASTE_TILE_URL_HERE",
  {
    opacity: 0.75,
    attribution:
      "USGS Landsat 8-9, processed with Google Earth Engine"
  }
);

lstLayer.addTo(map);
```

`getMapId()` generates a tile endpoint suitable for prototyping.

For the operational version:

```text
Earth Engine
Cloud Optimized GeoTIFF
Cloud Storage
Tile server, for example TiTiler
Leaflet HTML
```

## USGS M2M Alternative

Use the USGS Machine-to-Machine API only if the goal is to download original Landsat products.

Recommended API by task:

| Need | Recommended API |
| --- | --- |
| Process LST, masks, composites and statistics | Earth Engine Python API |
| Download original Landsat scenes | USGS M2M API |
| Create final HTML layer | Earth Engine COG export plus tile service |

For the operational product, start from the Earth Engine Python API because it avoids local download, extraction, mosaicking and manual management of many Landsat scenes.

## References To Verify During Implementation

- Google Earth Engine Landsat 8 Level 2 Collection 2 Tier 1 dataset documentation.
- USGS Landsat Collection 2 Surface Temperature documentation.
- USGS Landsat 8-9 OLI/TIRS Collection 2 Level-2 Science Products documentation.
- Google Earth Engine Python installation documentation.
- Google Earth Engine authentication and initialization documentation.
- Google Earth Engine `Export.image.toDrive` documentation.
- Google Earth Engine `ee.Image.getMapId` documentation.
- USGS Machine-to-Machine API documentation.
