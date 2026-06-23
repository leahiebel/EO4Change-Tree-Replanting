# EO4Change Group 4 — Post-Fire Forest Recovery Monitor

Sentinel-2 + Landsat 8/9 pipeline that maps per-pixel forest recovery after wildfires using a Relative Recovery Indicator (RRI) and Land Surface Temperature. Built for the DTU Space EO4Change course, June 2026.

---

## What it does

1. Builds cloud-free seasonal composites from Sentinel-2 for three periods: **pre-fire**, **post-fire (fire scar)**, and **recent**.
2. Computes a per-pixel **RRI** — how much of the pre-fire NDVI has been restored:
   ```
   RRI = (recent_NDVI − postfire_NDVI) / (prefire_NDVI − postfire_NDVI)
   ```
3. Classifies pixels into a **recovery class**: recovering well / recovering weakly / not recovering / outside assessment.
4. Generates time-series plots (RRI, NDVI, NDRE, NBR) and a pixel-count histogram by class.
5. Exports an interactive **HTML map** with all layers, a dynamic legend, and a Chart.js time-series panel.
6. Optionally adds a **Land Surface Temperature** layer (Landsat 8/9 Collection 2 L2SP).

---

## Study sites

| Region | Config | Reference event |
|--------|--------|-----------------|
| Rugballegård Skov, Denmark | `data/config_DK.yaml` | Tree planting April 2022 |
| Pedrógão Grande, Portugal | `data/config_pt.yaml` | Wildfire June 2017 |

---

## Setup

```bash
# Clone and enter the repo
git clone https://github.com/leahiebel/EO4Change-Tree-Replanting.git
cd EO4Change-Tree-Replanting

# Create and activate a virtual environment
python -m venv ../venv
source ../venv/bin/activate          # macOS/Linux
# ../venv/Scripts/activate           # Windows

# Install dependencies
pip install earthengine-api folium pyyaml pandas matplotlib python-dateutil

# Authenticate with Google Earth Engine (one-time)
python -c "import ee; ee.Authenticate()"
```

---

## Running

```bash
# Denmark — default config
python src/gee.py

# Portugal
python src/gee.py --config config_pt.yaml

# Second run is fast: time-series loaded from cached CSV automatically
python src/gee.py --config config_pt.yaml

# Force re-extraction from GEE (discard the cache)
python src/gee.py --config config_pt.yaml --force-ts

# Skip time-series entirely (map only, no plots)
python src/gee.py --config config_pt.yaml --map-only
```

> **Speed tip:** The first full run extracts all time windows from GEE (slow, ~1 min/window).
> Every subsequent run loads the cached `timeseries_{region}.csv` instantly and skips GEE.
> Use `--force-ts` to refresh the cache after changing the temporal config.

---

## Outputs

| File | Description |
|------|-------------|
| `src/timeseries_{region}.csv` | Per-window AOI-mean values for all bands + RRI (cached) |
| `src/timeseries_rri_{region}.png` | RRI time series |
| `src/timeseries_ndvi_{region}.png` | NDVI time series |
| `src/timeseries_ndre_{region}.png` | NDRE time series |
| `src/timeseries_nbr_{region}.png` | NBR time series |
| `src/histogram_recovery_class_{region}.png` | Pixel counts by recovery class |
| `src/map_{region}.html` | Interactive Folium map — open in any browser |

### HTML map

- Toggle any layer → the **legend** in the bottom-right updates automatically.
- Toggle any layer → the **time-series panel** (top-left) shows the matching index over time with a reference-event line (RRI layer → RRI chart; NDVI layer → NDVI chart; etc.).
- Click a **red cluster marker** (failed-pixel diagnostics) → popup with reasoning and an interactive time-series chart for that location.

---

## Config

All pipeline parameters live in a YAML file under `data/`. Key sections:

```yaml
spatial:
  type: geojson
  geojson_path: area_pt.geojson
  region_name: portugal-aoi

temporal:
  start: "2018-01-02"
  end:   "2024-01-01"
  cadence:
    type: fixed
    interval: monthly

comparison:
  prefire_years:  [2017, 2017]   # pre-fire composite window
  prefire_months: [4, 6]
  postfire_years: [2017, 2017]   # fire-scar composite window
  recent_years:   [2022, 2025]   # recent recovery composite
  months:         [5, 9]

classification:
  rri_good_threshold:      0.80  # ≥ this → recovering well
  rri_low_threshold:       0.30  # < this → not recovering
  burn_severity_threshold: 0.10  # dNBR floor to restrict to burned pixels
```

---

## Architecture

```
gee_config.py       CLI args, YAML parsing, GEE init, AOI, date windows
gee_s2.py           Sentinel-2 composite builders, spectral indices, biophys proxies
gee_recovery.py     RRI, recovery_class, guard_flags, forest mask
gee_diagnostics.py  Failed-cluster sampling, TS extraction, diagnostic plots
gee_sensitivity.py  Threshold sensitivity analysis
gee_map.py          Folium map builder, dynamic legend, Chart.js TS panel
landsat_lst.py      Landsat 8/9 LST median, anomaly, and stress class
gee.py              Orchestrates all of the above
```

---

## Spectral indices

| Index | Formula (S2 SR) | Interpretation |
|-------|----------------|----------------|
| NDVI  | (B8−B4)/(B8+B4) | Vegetation greenness |
| NDRE  | (B8A−B5)/(B8A+B5) | Chlorophyll / canopy quality |
| NDMI  | (B8−B11)/(B8+B11) | Canopy water content |
| BSI   | ((B11+B4)−(B8+B2))/((B11+B4)+(B8+B2)) | Bare soil exposure |
| NDWI  | (B3−B8)/(B3+B8) | Open water |
| NBR   | (B8−B12)/(B8+B12) | Burn severity / fire recovery |

---

## Recovery classes

| Class | Color | Meaning |
|-------|-------|---------|
| 0 | Grey | Uncertain — masked, water, or insufficient data |
| 1 | Dark green | Recovering well (RRI ≥ 0.80, NDRE confirms canopy) |
| 2 | Light green | Recovering but weak |
| 3 | Red | Not recovering (RRI < 0.30) |
| 4 | Brown | Outside burn recovery assessment (unburned) |

---

## Citation / acknowledgements

- Sentinel-2 data: ESA Copernicus `COPERNICUS/S2_SR_HARMONIZED` via Google Earth Engine
- Cloud masking: Cloud Score Plus (`GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED`)
- Forest mask: Dynamic World V1 (`GOOGLE/DYNAMICWORLD/V1`)
- LST: Landsat Collection 2 Level-2 (`LANDSAT/LC08/C02/T1_L2`, `LANDSAT/LC09/C02/T1_L2`)
- RRI methodology adapted from Bright et al. (2019), *Remote Sensing of Environment*
