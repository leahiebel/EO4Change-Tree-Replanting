# EO4C Landsat Land Surface Temperature Workflow

This project generates production-oriented Land Surface Temperature (LST) layers for a mainland Portugal reforestation-monitoring domain using Landsat 8/9 Collection 2 Level-2 Science Products in Google Earth Engine.

## Architecture

All raster processing runs server-side in Earth Engine. Local Python is used for configuration, vector validation, export orchestration, metadata, and web-layer manifests.

- `src/vector_io.py`: validates AOI, forest and optional control polygons with GeoPandas.
- `src/landsat_collection.py`: filters Landsat 8/9 L2SP scenes and merges them after identical preprocessing.
- `src/quality_mask.py`: masks fill, clouds, cirrus, shadows, snow, water, saturation, cloud proximity and high ST uncertainty.
- `src/temperature.py`: converts `ST_B10` digital numbers to `LST_C`.
- `src/composites.py`: builds target, baseline, anomaly, z-score and stress-anomaly rasters.
- `src/forest_analysis.py`: exports forest/control statistics and yearly time series.
- `src/exports.py`: starts COG and CSV Earth Engine batch exports.
- `src/web_layers.py`: creates prototype EE tiles and operational COG manifests.

## Installation

```bash
cd /Users/nicolocaron/Documents/GitHub/EO4C/project
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The existing EO4C sandbox can also be used:

```bash
cd /Users/nicolocaron/Documents/GitHub/EO4C
source .venv-eo/bin/activate
python -m pip install -r project/requirements.txt
```

## Configuration

Edit `config/config.yaml` and replace:

- `earth_engine_project`
- `aoi_path`
- `forest_path`
- `control_path`, or keep `null` to skip paired cooling analysis
- `output_directory`
- `gcs_bucket`, when exporting operational COGs to Google Cloud Storage

AOI, forest and control vectors must be polygon or multipolygon layers with a valid CRS. They are converted to EPSG:4326 before Earth Engine conversion.

## Earth Engine Authentication

Interactive user authentication:

```bash
python run_pipeline.py --config config/config.yaml
```

For an already authenticated environment:

```bash
python run_pipeline.py --config config/config.yaml --no-auth
```

The code uses:

```python
ee.Authenticate()
ee.Initialize(project=EARTH_ENGINE_PROJECT)
```

No credentials or service-account secrets are stored in source code.

## Execution

```bash
cd project
source .venv/bin/activate
python run_pipeline.py --config config/config.yaml
```

The script starts Earth Engine batch exports for:

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
- `forest_control_statistics.csv`
- `yearly_lst_timeseries.csv`

Set `wait_for_exports: true` to block and fail locally if an Earth Engine task fails.

## COG Export And Publication

When `gcs_bucket` is configured, COGs are exported to:

```text
gs://<gcs_bucket>/<gcs_prefix>/
```

Make the COGs readable by the tile service account or public, depending on your deployment. `data/outputs/layer_manifest.json` contains TiTiler-style layer definitions. Replace `titiler_base_url` with your tile service endpoint.

When `gcs_bucket` is `null`, exports go to the configured Google Drive folder. Download or publish those COGs separately, then update `layer_manifest.json` with public COG URLs.

## Web Prototype

The pipeline writes:

- `data/outputs/style_config.json`
- `data/outputs/layer_manifest.json`
- `data/outputs/prototype_tile_layers.json`, when enabled
- `web/aoi.geojson`
- `web/forest.geojson`

Serve the project directory over HTTP:

```bash
cd project
python -m http.server 8000
```

Open:

```text
http://localhost:8000/web/index.html
```

The HTML uses OpenStreetMap, layer controls, opacity, legend, scale bar, AOI and forest boundaries, loading/error messages, attribution and optional point-value queries when the tile service supports them.

## Tests

```bash
cd project
source .venv/bin/activate
pytest
```

The tests cover vector validation, ST_B10 Celsius conversion, QA bitmask rules and metadata limitations.

## Scientific Limitations

LST is surface temperature, not 2-m air temperature. Outputs are distributed on a 30-m grid, but native Landsat TIRS thermal sampling is approximately 100 m; resampling does not create independent 30-m thermal information. Cloud proximity and emissivity uncertainty can affect LST. Young or narrow plantations may be mixed with surrounding land cover. The thermal-stress class is a relative surface-temperature anomaly indicator, not a direct measurement of plant physiological stress.

