# EO4C Objective: Landsat Land Surface Temperature Workflow

## Scope Decision

This project will work only with Landsat 8 and Landsat 9 Collection 2 Level-2 Science Products for Land Surface Temperature.

Authoritative Earth Engine collections:

- `LANDSAT/LC08/C02/T1_L2`
- `LANDSAT/LC09/C02/T1_L2`

Only scenes with:

- `PROCESSING_LEVEL == "L2SP"`

may contribute to Land Surface Temperature products, because `ST_B10` is masked for `L2SR` products.

Sentinel-2 is not used for the LST workflow. Sentinel-2 NDVI, NDRE and NDMI may be added later as separate vegetation or moisture layers to interpret the Landsat LST products, but they are not part of the thermal-temperature calculation.

## Project Objective

Build a complete, modular, production-oriented Python workflow to generate a Land Surface Temperature layer for a large reforestation-monitoring domain in mainland Portugal.

The final raster products must be suitable for:

- scientific analysis;
- forest versus control-area temperature comparison;
- interactive visualisation in a Leaflet-based HTML application.

The workflow must assess:

1. The spatial distribution of Land Surface Temperature, LST, across the Portuguese study domain.
2. Whether forested or reforested areas are cooler than comparable surrounding non-forested areas.
3. Whether unusually high surface temperatures may indicate potential thermal stress in reforested areas.

The workflow must distinguish between:

- absolute Land Surface Temperature, in degrees Celsius;
- temporal LST anomaly, in degrees Celsius;
- forest cooling difference relative to a control area, in degrees Celsius;
- observation confidence and number of valid acquisitions.

LST alone must not be interpreted as proof of physiological plant stress.

## Technology Requirements

Use:

- Python 3.11 or newer;
- Google Earth Engine Python API;
- GeoPandas for local vector input validation;
- Rasterio only for local raster inspection or validation;
- Pandas for output statistics;
- Pathlib for paths;
- Logging instead of print statements where appropriate;
- type hints and docstrings;
- configuration separated from processing logic.

Authentication must support:

```python
ee.Authenticate()
ee.Initialize(project=EARTH_ENGINE_PROJECT)
```

No credentials or service-account secrets may be placed in source code.

## Configuration Inputs

The workflow must expose a YAML configuration with:

```yaml
earth_engine_project: "replace-with-project-id"

aoi_path: "/path/to/portugal_domain.geojson"
forest_path: "/path/to/reforestation_polygons.geojson"
control_path: "/path/to/control_polygons.geojson"

start_date: "2013-04-01"
end_date: "2026-12-31"

season_start_month: 6
season_end_month: 9

baseline_start_year: 2013
baseline_end_year: 2019
target_year: 2025

max_scene_cloud_cover: 70

min_cloud_distance_km: 1.0
max_st_uncertainty_k: 2.0

mask_water: true

output_crs: "EPSG:3763"
export_scale_m: 30
nodata_value: -9999

output_directory: "/path/to/output"

gcs_bucket: null
```

`control_path` may be optional. If no control polygons are provided, scientific forest-cooling analysis must be skipped or clearly labelled exploratory.

## Vector Input Handling

Read AOI, forest polygons and optional control polygons using GeoPandas.

Validation must ensure:

- files exist;
- geometries are non-empty;
- geometries are polygon or multipolygon;
- invalid geometries are repaired where safely possible;
- all layers are converted to EPSG:4326 before conversion to Earth Engine geometries;
- forest and control polygons intersect the AOI;
- duplicated geometries are removed;
- informative exceptions are raised when validation fails.

Preserve useful identifiers when present:

- `project_id`
- `site_name`
- `planting_year`
- `species`
- `management_type`

If these fields do not exist, create a stable feature identifier.

## Landsat Collection Filtering

For both Landsat 8 and Landsat 9:

1. Filter by AOI.
2. Filter by date.
3. Filter `CLOUD_COVER` to the configured maximum.
4. Filter `PROCESSING_LEVEL == "L2SP"`.

Retain at least:

- `ST_B10`
- `ST_QA`
- `ST_CDIST`
- `QA_PIXEL`
- `QA_RADSAT`

Preserve image properties:

- `system:time_start`
- `SPACECRAFT_ID`
- `LANDSAT_PRODUCT_ID`
- `PROCESSING_LEVEL`
- `CLOUD_COVER`
- `WRS_PATH`
- `WRS_ROW`

Landsat 8 and Landsat 9 must be merged only after applying identical preprocessing.

## Quality Mask

Implement a reusable:

```python
mask_landsat_lst(image)
```

Use `QA_PIXEL` to exclude:

- fill pixels;
- dilated clouds;
- high-confidence cirrus;
- clouds;
- cloud shadows;
- snow and ice;
- water when `mask_water` is true.

Use `QA_RADSAT` to exclude radiometrically saturated or terrain-occluded pixels.

Use `ST_CDIST` with its scale factor to exclude pixels closer to clouds than:

```text
min_cloud_distance_km
```

Use `ST_QA` with its scale factor to exclude pixels whose estimated surface-temperature uncertainty exceeds:

```text
max_st_uncertainty_k
```

Masking criteria must be configurable and documented.

Do not apply arbitrary smoothing before calculating composites.

## Temperature Conversion

Convert `ST_B10` digital numbers to Kelvin using:

```text
LST_K = ST_B10 * 0.00341802 + 149.0
```

Convert Kelvin to Celsius using:

```text
LST_C = LST_K - 273.15
```

Rename the output band:

```text
LST_C
```

Retain quality bands:

- `ST_uncertainty_K`
- `cloud_distance_km`
- `valid_observation`

Do not clamp scientific LST rasters to display ranges. Min/max limits are visualisation-only.

## Temporal Products

### Per-Acquisition Collection

For every valid acquisition, retain:

- `LST_C`;
- ST uncertainty;
- acquisition date;
- satellite identifier;
- valid-pixel mask.

### Target-Year Seasonal Composite

For the configured seasonal months in `target_year`, calculate:

- median LST;
- mean LST;
- 10th percentile;
- 90th percentile;
- standard deviation;
- number of valid observations;
- median ST uncertainty.

Use the median as the primary visual layer.

### Baseline Climatology

For the baseline period and same seasonal months, calculate:

- baseline median LST;
- baseline mean LST;
- baseline standard deviation;
- baseline valid-observation count.

### Thermal Anomaly

Calculate:

```text
LST_anomaly_C =
target_year seasonal median LST
-
baseline seasonal median LST
```

### Yearly Time Series

Optionally calculate yearly seasonal median LST over:

- complete AOI;
- every forest polygon;
- every control polygon.

Export as tidy CSV.

## Forest Cooling Analysis

For each forest or reforestation polygon, calculate:

- median forest LST;
- mean forest LST;
- standard deviation;
- valid pixel count;
- valid acquisition count;
- median uncertainty.

For the corresponding control polygon, calculate the same statistics.

Calculate:

```text
forest_cooling_C =
median_LST_forest
-
median_LST_control
```

Interpretation:

- `forest_cooling_C < 0`: forest is cooler than the control.
- `forest_cooling_C > 0`: forest is warmer than the control.

Do not describe the difference as a causal forest effect unless control areas are demonstrably comparable.

If `planting_year` is available, also calculate:

- pre-planting seasonal LST;
- post-planting seasonal LST;
- corresponding control change;
- difference-in-differences estimate.

Use:

```text
effect_C =
(post_forest - pre_forest)
-
(post_control - pre_control)
```

Report the number of valid observations supporting each estimate.

## Thermal-Stress Anomaly Indicator

Create a cautious anomaly layer based on relative anomalies, not on a universal Celsius threshold.

Calculate:

```text
LST_z =
(target_LST - baseline_mean_LST)
/
baseline_standard_deviation
```

Handle zero baseline standard deviation explicitly.

Create categorical values:

- `0`: no valid data;
- `1`: near-normal temperature;
- `2`: moderately warm anomaly;
- `3`: strong warm anomaly.

Z-score thresholds must be configurable.

Output band:

```text
thermal_stress_anomaly_class
```

Metadata warning:

```text
This is a surface-temperature anomaly indicator, not a direct measurement of plant physiological stress.
```

## Required Scientific Outputs

Export rasters as Cloud Optimized GeoTIFFs:

- `01_lst_median_target_C.tif`
- `02_lst_mean_target_C.tif`
- `03_lst_p10_target_C.tif`
- `04_lst_p90_target_C.tif`
- `05_lst_std_target_C.tif`
- `06_lst_valid_observations.tif`
- `07_lst_uncertainty_K.tif`
- `08_lst_baseline_median_C.tif`
- `09_lst_anomaly_C.tif`
- `10_lst_zscore.tif`
- `11_thermal_stress_class.tif`

Export settings:

- scale: 30 m;
- CRS: configured output CRS, default `EPSG:3763`;
- region: AOI;
- NoData: `-9999`;
- `cloudOptimized: true`;
- continuous rasters as Float32;
- categorical rasters as integer where possible.

Also export:

- `forest_control_statistics.csv`
- `yearly_lst_timeseries.csv`
- `processing_metadata.json`
- `layer_manifest.json`
- `style_config.json`

## Metadata Requirements

`processing_metadata.json` must include:

- processing date;
- Earth Engine collections;
- Landsat satellites used;
- requested date range;
- actual acquisition dates used;
- seasonal window;
- baseline period;
- AOI area;
- export CRS;
- export scale;
- scale and offset applied to `ST_B10`;
- QA mask rules;
- cloud-distance threshold;
- uncertainty threshold;
- number of Landsat 8 images;
- number of Landsat 9 images;
- number of valid images;
- percentage of AOI with valid observations;
- known limitations.

Known limitations must explicitly state:

1. LST is surface temperature, not 2-m air temperature.
2. The output is distributed on a 30-m grid.
3. The native TIRS thermal sampling is approximately 100 m.
4. Resampling to 30 m does not create independent 30-m thermal information.
5. Cloud proximity and emissivity uncertainty can affect LST.
6. Young or narrow plantations may be mixed with surrounding land cover.

## Web-Layer Architecture

Final products must support a Leaflet-based HTML map.

### Prototype Mode

Create Earth Engine visualisation tile URLs using `getMapId` for:

- target median LST;
- LST anomaly;
- thermal-stress class.

Return XYZ URL templates and visualisation parameters.

Do not expose private credentials in HTML.

Prototype mode is for development and testing only.

### Operational Mode

Use Cloud Optimized GeoTIFFs as persistent source files.

Support export to Google Cloud Storage when `gcs_bucket` is configured.

Generate a manifest compatible with a COG tile service such as TiTiler:

```json
{
  "id": "lst_median_2025",
  "title": "Median land surface temperature, summer 2025",
  "type": "raster",
  "unit": "°C",
  "cog_url": "https://storage.../01_lst_median_target_C.tif",
  "tile_url": "https://tile-server/.../{z}/{x}/{y}.png",
  "legend": {},
  "bounds": [],
  "min_zoom": 5,
  "max_zoom": 14,
  "opacity": 0.75,
  "attribution": "USGS Landsat 8-9 Collection 2 Level-2"
}
```

Leaflet example must include:

- OpenStreetMap basemap;
- layer control;
- LST layer;
- anomaly layer;
- stress-class layer;
- opacity control;
- legend;
- scale bar;
- AOI boundary;
- forest polygon boundary;
- mouse-click or popup value query when supported;
- loading and error messages;
- attribution.

Do not embed GeoTIFFs in HTML.

## Visualisation Styles

`style_config.json` must define display-only styles.

Target LST:

```yaml
min: 15
max: 45
palette:
  - "#313695"
  - "#4575b4"
  - "#74add1"
  - "#abd9e9"
  - "#e0f3f8"
  - "#ffffbf"
  - "#fee090"
  - "#fdae61"
  - "#f46d43"
  - "#d73027"
  - "#a50026"
```

Anomaly:

```yaml
min: -5
max: 5
palette: diverging centred on 0
```

Categorical stress:

- `0`: transparent/no data;
- `1`: normal;
- `2`: moderate warm anomaly;
- `3`: strong warm anomaly.

Styles must not alter scientific raster values.

## Required Project Structure

```text
project/
├── config/
│   └── config.yaml
├── data/
│   ├── vectors/
│   └── outputs/
├── src/
│   ├── __init__.py
│   ├── authentication.py
│   ├── vector_io.py
│   ├── landsat_collection.py
│   ├── quality_mask.py
│   ├── temperature.py
│   ├── composites.py
│   ├── forest_analysis.py
│   ├── exports.py
│   ├── web_layers.py
│   └── metadata.py
├── web/
│   ├── index.html
│   ├── app.js
│   └── styles.css
├── tests/
│   ├── test_vector_io.py
│   ├── test_temperature_conversion.py
│   ├── test_qa_mask.py
│   └── test_metadata.py
├── requirements.txt
├── README.md
└── run_pipeline.py
```

## Error Handling Requirements

The workflow must detect and clearly report:

- missing AOI files;
- invalid geometry;
- no Landsat scenes found;
- scenes found but no L2SP temperature data;
- zero valid pixels after masking;
- missing control polygons;
- failed Earth Engine export task;
- Google Cloud Storage permission errors;
- invalid date ranges;
- insufficient baseline observations;
- division by zero when baseline standard deviation is zero.

Do not silently continue when a scientifically important output cannot be calculated.

## Performance Requirements

The AOI may cover approximately 1,300 square kilometres.

All raster processing must occur server-side in Earth Engine.

Do not call `getInfo()` on large images, image collections or full pixel arrays.

Use `getInfo()` only for:

- small metadata dictionaries;
- task status;
- summary statistics.

Use Earth Engine batch exports for raster outputs.

Avoid downloading original Landsat scenes locally.

## Acceptance Tests

The final implementation is accepted only if:

1. Landsat 8 and Landsat 9 are merged correctly.
2. Only L2SP scenes contribute to LST.
3. `ST_B10` is converted correctly to Celsius.
4. Clouds, shadows, snow, fill and saturated pixels are masked.
5. The target composite contains valid values inside the AOI.
6. The COG has expected CRS, bounds, NoData and band count.
7. Output raster values are not unintentionally clipped by visualisation limits.
8. The valid-observation raster is produced.
9. The forest-control CSV is produced when control polygons exist.
10. The HTML displays the LST layer through an XYZ tile source.
11. No authentication secret appears in frontend code.
12. The README explains setup, authentication, execution and publication.
13. Metadata clearly describes the 30-m output grid and approximately 100-m native thermal sampling.
