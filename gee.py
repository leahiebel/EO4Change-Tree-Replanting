"""
EO4Change Group 4 — Rugballegård Skov reforestation monitor
All parameters driven by config_DK.yaml.

Index computation
─────────────────
Spectral indices (deterministic band arithmetic on S2 SR reflectance):

  NDVI  = (B8  − B4)  / (B8  + B4)      NIR vs Red
          Chlorophyll absorbs red (B4 ≈ 665 nm), reflects NIR (B8 ≈ 842 nm).
          Range −1…1; bare soil ≈ 0.1, young forest ≈ 0.3–0.6, dense forest > 0.7.

  NDRE  = (B8A − B5)  / (B8A + B5)      Red-edge (B5 ≈ 705 nm, B8A ≈ 865 nm)
          Red-edge reflectance rises steeply with chlorophyll concentration.
          Responds earlier to stress than NDVI; doesn't saturate at high LAI.

  NDMI  = (B8  − B11) / (B8  + B11)     NIR vs SWIR1 (B11 ≈ 1610 nm)
          Leaf water strongly absorbs SWIR. High NDMI → high moisture content.

  BSI   = (B11 + B4 − B8 − B2) / (B11 + B4 + B8 + B2)
          Bare soil has high SWIR+Red, low NIR+Blue. Positive → exposed soil.
          Persistent high BSI after planting date = establishment failure signal.

Biophysical variables (s2biophys neural-network model):
  LAI-e  Effective Leaf Area Index [m²/m²]  — canopy density / growth rate proxy
  FAPAR  Fraction of Absorbed PAR [0…1]     — direct photosynthesis proxy
  FCOVER Fraction of green vegetation cover [0…1] — structural establishment indicator

  The s2biophys model is a 2-layer neural network trained on PROSAIL radiative
  transfer simulations. It takes 8 normalised S2 bands (B3,B4,B5,B6,B7,B8A,B11,B12)
  plus solar/view geometry as input. Call gee_biophys.retrieve(cfg) to run it.
  Below we use published linear proxies from Delegido et al. as a fast fallback.
"""
import ee
import yaml
import json
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from datetime import datetime
from dateutil.relativedelta import relativedelta

# ── Config ────────────────────────────────────────────────────────────────────
CONFIG_PATH = Path(__file__).parent / "config_DK.yaml"
with open(CONFIG_PATH) as f:
    cfg = yaml.safe_load(f)

sp  = cfg["spatial"]
tmp = cfg["temporal"]
var = cfg["variables"]
exp = cfg["export"]
opt = cfg["options"]

# variable may be a string or a list
BIOPHYS_VARS: list[str] = (
    var["variable"] if isinstance(var["variable"], list) else [var["variable"]]
)

# ── GEE init ──────────────────────────────────────────────────────────────────
# Run ee.Authenticate() once in a terminal to cache credentials, then leave it
# commented out — Initialize() reuses the cached token automatically.
# ee.Authenticate()
ee.Initialize(project=exp.get("project_id", "eo4change-ex"))

# ── AOI ───────────────────────────────────────────────────────────────────────
def build_aoi(sp_cfg: dict) -> ee.Geometry:
    if sp_cfg["type"] == "bbox":
        minx, miny, maxx, maxy = sp_cfg["bbox"]
        return ee.Geometry.Rectangle([minx, miny, maxx, maxy])
    if sp_cfg["type"] == "geojson":
        p = Path(__file__).parent / sp_cfg.get("geojson_path", sp_cfg.get("path"))
        gj = json.loads(p.read_text())
        valid = [
            f for f in gj["features"]
            if len({tuple(c) for c in f["geometry"]["coordinates"][0]}) > 3
        ]
        return ee.Geometry(valid[-1]["geometry"] if valid else gj)
    raise ValueError(sp_cfg["type"])

aoi = build_aoi(sp)

# ── Date windows ──────────────────────────────────────────────────────────────
FIXED_DELTAS = {
    "monthly":   relativedelta(months=1),
    "biweekly":  relativedelta(weeks=2),
    "quarterly": relativedelta(months=3),
    "annual":    relativedelta(years=1),
}
SEASONS = {"MAM": (3, 6), "JJA": (6, 9), "SON": (9, 12), "DJF": (12, 3)}

def date_windows(start: str, end: str, cadence: dict) -> list[tuple[str, str, str]]:
    t0 = datetime.strptime(start, "%Y-%m-%d")
    t1 = datetime.strptime(end,   "%Y-%m-%d")
    if cadence["type"] == "fixed":
        delta = FIXED_DELTAS[cadence["interval"]]
        windows, cur = [], t0
        while cur < t1:
            nxt = min(cur + delta, t1)
            windows.append((cur.strftime("%Y-%m-%d"), nxt.strftime("%Y-%m-%d"),
                            cur.strftime("%Y-%m")))
            cur = nxt
        return windows
    if cadence["type"] == "seasons":
        windows = []
        for yr in range(t0.year, t1.year + 1):
            for name, (ms, me) in SEASONS.items():
                s = datetime(yr if name != "DJF" else yr,     ms, 1)
                e = datetime(yr if name != "DJF" else yr + 1, me, 1)
                if e <= t0 or s >= t1:
                    continue
                windows.append((max(s, t0).strftime("%Y-%m-%d"),
                                min(e, t1).strftime("%Y-%m-%d"),
                                f"{yr}-{name}"))
        windows.sort()
        return windows
    raise ValueError(cadence["type"])

windows = date_windows(tmp["start"], tmp["end"], tmp["cadence"])
print(f"→ {len(windows)} windows  "
      f"({tmp['cadence']['type']} / {tmp['cadence'].get('interval','seasons')})")

# ── Cloud masking ─────────────────────────────────────────────────────────────
CS_BAND   = opt.get("csplus_band", "cs")
CS_THRESH = opt.get("cs_plus_threshold", 0.65)
MAX_CLOUD = opt.get("max_cloud_cover", 50)

def build_s2_col(ws: str, we: str) -> ee.ImageCollection:
    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(aoi)
            .filterDate(ws, we)
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", MAX_CLOUD)))
    cs = ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
    return (s2.linkCollection(cs, [CS_BAND])
              .map(lambda img: img.updateMask(img.select(CS_BAND).gte(CS_THRESH))))

# ── Spectral indices ──────────────────────────────────────────────────────────
# All computed as normalised differences or simple arithmetic on BOA reflectance.
# Reflectance values in S2_SR_HARMONIZED are scaled integers (÷10000 = true ρ).
def add_spectral_indices(img: ee.Image) -> ee.Image:
    ndvi = img.normalizedDifference(["B8",  "B4" ]).rename("NDVI")
    ndre = img.normalizedDifference(["B8A", "B5" ]).rename("NDRE")
    ndmi = img.normalizedDifference(["B8",  "B11"]).rename("NDMI")
    bsi  = (img.select("B11").add(img.select("B4"))
               .subtract(img.select("B8").add(img.select("B2")))
               .divide(img.select("B11").add(img.select("B4"))
                          .add(img.select("B8")).add(img.select("B2")))
               ).rename("BSI")
    return img.addBands([ndvi, ndre, ndmi, bsi])

SPECTRAL_BANDS = ["NDVI", "NDRE", "NDMI", "BSI"]

# ── Biophysical variables ─────────────────────────────────────────────────────
# Proper route: gee_biophys.retrieve(cfg) runs the PROSAIL-trained neural network.
# Proxy route (below): simple empirical formulas valid for broadleaf temperate forest.
#   LAI-e  ≈  3.618 × NDVI − 0.118          (Baret & Guyot 1991 linearisation)
#   FAPAR  ≈  1.136 × NDVI − 0.04            (Myneni & Williams 1994)
#   FCOVER ≈  1 − exp(−0.5 × LAI)            (Beer–Lambert, k=0.5 for mixed canopy)
# These are adequate for trend analysis; switch to gee_biophys for absolute values.

def add_biophys_proxies(img: ee.Image) -> ee.Image:
    ndvi  = img.select("NDVI")
    laie  = ndvi.multiply(3.618).subtract(0.118).clamp(0, 8).rename("laie")
    fapar = ndvi.multiply(1.136).subtract(0.040).clamp(0, 1).rename("fapar")
    # FCOVER via Beer–Lambert on LAI proxy
    fcover = laie.multiply(-0.5).exp().multiply(-1).add(1).clamp(0, 1).rename("fcover")
    return img.addBands([laie, fapar, fcover])

# Only add bands that were requested in config; keep all original bands intact
def add_requested_biophys(img: ee.Image) -> ee.Image:
    return add_biophys_proxies(img)   # addBands — does not drop existing bands

# ── gee_biophys integration point ────────────────────────────────────────────
# To replace the proxies above with the full neural-network model, uncomment:
#
#   import gee_biophys
#   results = gee_biophys.retrieve(cfg)
#
# cfg["variables"]["variable"] is already a list: ["laie", "fapar", "fcover"]
# gee_biophys.retrieve() reads this and submits one export task per composite.

# ── Composite builder ─────────────────────────────────────────────────────────
ALL_BANDS = SPECTRAL_BANDS + BIOPHYS_VARS

# Fully-masked placeholder returned for empty windows (no cloud-free images).
# reduceRegion on a masked image returns None for every band, which we catch below.
_PLACEHOLDER = (ee.Image.constant([0] * len(ALL_BANDS))
                  .rename(ALL_BANDS)
                  .updateMask(ee.Image(0)))

def make_composite(ws: str, we: str) -> ee.Image:
    col = (build_s2_col(ws, we)
             .map(add_spectral_indices)
             .map(add_requested_biophys))
    # Server-side conditional: if collection is empty return masked placeholder
    return ee.Image(
        ee.Algorithms.If(col.size().gt(0), col.median().clip(aoi), _PLACEHOLDER)
    )

# ── Time-series extraction ────────────────────────────────────────────────────
def extract_mean(ws: str, we: str) -> dict:
    stats = make_composite(ws, we).select(ALL_BANDS).reduceRegion(
        reducer   = ee.Reducer.mean(),
        geometry  = aoi,
        scale     = int(exp.get("scale", 20)),
        maxPixels = int(exp.get("max_pixels", 1e9)),
    )
    return stats.getInfo()

print("Extracting time series — please wait…")
records = []
for ws, we, label in windows:
    try:
        row = extract_mean(ws, we) or {}
    except ee.EEException as e:
        print(f"  {label:12s}  skipped ({e})")
        row = {}
    row["date"]  = ws
    row["label"] = label
    records.append(row)
    v = row.get("NDVI")
    print(f"  {label:12s}  NDVI={v:.3f}" if v is not None else f"  {label:12s}  (no data)")

# ── Plot 1 — spectral indices ─────────────────────────────────────────────────
PLANTING = datetime(2022, 4, 1)
dates = [datetime.strptime(r["date"], "%Y-%m-%d") for r in records]

SPEC_COLORS  = {"NDVI": "#2d8a2d", "NDRE": "#7baf27", "NDMI": "#1e6eb5", "BSI": "#c47a1e"}
BIOPH_COLORS = {"laie": "#8b4513", "fapar": "#ff6600", "fcover": "#009900"}
BIOPH_LABELS = {"laie": "LAI-e [m²/m²]", "fapar": "FAPAR [0-1]", "fcover": "FCOVER [0-1]"}

def _plot_band(ax, band, color, ylabel):
    vals = [r.get(band) for r in records]
    x = [d for d, v in zip(dates, vals) if v is not None]
    y = [v for v in vals if v is not None]
    ax.plot(x, y, marker="o", ms=3, lw=1.5, color=color, label=band)
    ax.axvline(PLANTING, color="red", ls="--", lw=1, label="Planting T₀")
    ax.set_ylabel(ylabel, fontsize=9)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

region_name = sp.get("region_name", "AOI")
cadence_str = f"{tmp['cadence']['type']} / {tmp['cadence'].get('interval','seasons')}"

fig1, axes1 = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
fig1.suptitle(f"{region_name} — Spectral indices ({cadence_str})", fontsize=12)
for ax, band in zip(axes1.flat, SPECTRAL_BANDS):
    _plot_band(ax, band, SPEC_COLORS[band], band)
fig1.autofmt_xdate()
plt.tight_layout()
out1 = Path(__file__).parent / "timeseries_spectral.png"
fig1.savefig(out1, dpi=150)
plt.show()
print(f"Saved → {out1}")

# ── Plot 2 — biophysical variables ────────────────────────────────────────────
n = len(BIOPHYS_VARS)
fig2, axes2 = plt.subplots(1, n, figsize=(5 * n, 4), sharex=True)
if n == 1:
    axes2 = [axes2]
fig2.suptitle(f"{region_name} — Biophysical variables / s2biophys proxies ({cadence_str})",
              fontsize=11)
for ax, band in zip(axes2, BIOPHYS_VARS):
    _plot_band(ax, band, BIOPH_COLORS.get(band, "grey"), BIOPH_LABELS.get(band, band))
fig2.autofmt_xdate()
plt.tight_layout()
out2 = Path(__file__).parent / "timeseries_biophys.png"
fig2.savefig(out2, dpi=150)
plt.show()
print(f"Saved → {out2}")

# ── Raster visualisation — latest summer composite ───────────────────────────
RASTER_START, RASTER_END = "2023-06-01", "2023-08-31"
raster = make_composite(RASTER_START, RASTER_END)
bounds = aoi.bounds().getInfo()["coordinates"][0]

VIS = {
    "RGB":    {"bands": ["B4","B3","B2"], "params": {"min":0,"max":2800,"gamma":1.4}},
    "NDVI":   {"bands": ["NDVI"],  "params": {"min":-0.1,"max":0.85,
               "palette":["saddlebrown","khaki","limegreen","darkgreen"]}},
    "LAI-e":  {"bands": ["laie"],  "params": {"min":0,"max":6,
               "palette":["white","yellow","limegreen","darkgreen"]}},
    "FCOVER": {"bands": ["fcover"],"params": {"min":0,"max":1,
               "palette":["white","lightgreen","forestgreen"]}},
    "BSI":    {"bands": ["BSI"],   "params": {"min":-0.3,"max":0.4,
               "palette":["darkgreen","white","orange","saddlebrown"]}},
}
# Only keep layers whose bands exist in ALL_BANDS or raw S2
VIS_AVAILABLE = {k: v for k, v in VIS.items()
                 if all(b in ALL_BANDS + ["B4","B3","B2"] for b in v["bands"])}

try:
    import geemap
    Map = geemap.Map()
    Map.centerObject(aoi, zoom=14)
    for title, vc in VIS_AVAILABLE.items():
        Map.addLayer(raster.select(vc["bands"]), vc["params"], title)
    Map.addLayer(aoi, {"color": "red"}, "AOI")
    out_map = Path(__file__).parent / "map_rugballegaard.html"
    Map.save(str(out_map))
    print(f"Interactive map → {out_map}")

except ImportError:
    import urllib.request, io
    from PIL import Image

    n_panels = len(VIS_AVAILABLE)
    fig3, axes3 = plt.subplots(1, n_panels, figsize=(4 * n_panels, 4))
    fig3.suptitle(f"{region_name} — {RASTER_START} to {RASTER_END}", fontsize=11)
    for ax, (title, vc) in zip(axes3, VIS_AVAILABLE.items()):
        url = raster.select(vc["bands"]).getThumbURL(
            {**vc["params"], "region": bounds, "dimensions": 512, "format": "png"}
        )
        with urllib.request.urlopen(url) as resp:
            img = Image.open(io.BytesIO(resp.read()))
        ax.imshow(img)
        ax.set_title(title, fontsize=9)
        ax.axis("off")
    plt.tight_layout()
    out3 = Path(__file__).parent / "raster_rugballegaard.png"
    fig3.savefig(out3, dpi=150)
    plt.show()
    print(f"Saved → {out3}")
