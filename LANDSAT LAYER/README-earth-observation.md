# EO4C Earth Observation Sandbox

Local virtual environment:

```bash
source .venv-eo/bin/activate
```

Registered Jupyter kernel:

```text
Python (EO4C Earth Observation)
```

Kernel name:

```text
eo4c-earth-observation
```

Start JupyterLab from this workspace:

```bash
.venv-eo/bin/jupyter lab
```

Dependencies are listed in:

```text
requirements-earth-observation.txt
```

Core workflow support:

- Landsat 8/9 LST through Google Earth Engine: `earthengine-api`
- AOI GeoJSON/vector data: `geopandas`, `shapely`, `pyproj`, `fiona`
- GeoTIFF/raster inspection: `rasterio`, `rioxarray`, `xarray`
- Notebook plotting/mapping: `matplotlib`, `folium`, `mapclassify`

Primary project objective:

```text
project/OBJECTIVE.md
```

Current implementation task:

```text
project/CURRENT_TASK.md
```
