"""
EO4Change Group 4 — Sentinel-2 reforestation / vegetation monitor.
All parameters driven by a YAML config (default: data/config_DK.yaml).

Run for a different region with:
    python src/gee.py --config config_pt.yaml

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

  NDWI  = (B3  − B8)  / (B3  + B8)      Green vs NIR (McFeeters 1996)
          Open water reflects green and absorbs NIR. Positive → standing water /
          ponding / flooding. Distinct from NDMI (which measures *leaf* water).

  NBR   = (B8  − B12) / (B8  + B12)     NIR vs SWIR2 (B12 ≈ 2190 nm)
          Healthy vegetation: high NIR, low SWIR2 → high NBR.
          Burned vegetation (char/dry): low NIR, high SWIR2 → low / negative NBR.
          ΔNBR (pre minus post) is the standard fire-severity metric.

Biophysical variables (s2biophys neural-network model):
  LAI-e  Effective Leaf Area Index [m²/m²]  — canopy density / growth rate proxy
  FAPAR  Fraction of Absorbed PAR [0…1]     — direct photosynthesis proxy
  FCOVER Fraction of green vegetation cover [0…1] — structural establishment indicator

  The s2biophys model is a 2-layer neural network trained on PROSAIL radiative
  transfer simulations. It takes 8 normalised S2 bands (B3,B4,B5,B6,B7,B8A,B11,B12)
  plus solar/view geometry as input. Call gee_biophys.retrieve(cfg) to run it.
  Below we use published linear proxies from Delegido et al. as a fast fallback.
"""
import argparse
import ee
import yaml
import json
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from datetime import datetime
from dateutil.relativedelta import relativedelta

# ── Paths & CLI ───────────────────────────────────────────────────────────────
# Layout: <repo>/src/gee.py  +  <repo>/data/<config>.yaml  +  <repo>/data/<area>.geojson
SRC_DIR  = Path(__file__).resolve().parent
DATA_DIR = SRC_DIR.parent / "data"

_parser = argparse.ArgumentParser(add_help=True, description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
_parser.add_argument("--config", default="config_DK.yaml",
                     help="config filename inside data/ (or absolute path). "
                          "Default: config_DK.yaml")
_parser.add_argument("--map-only", action="store_true",
                     help="Skip the slow per-window time-series extraction "
                          "and plots; only build composites + change image + map.")
_args, _ = _parser.parse_known_args()
MAP_ONLY = _args.map_only

_cfg_arg = Path(_args.config)
CONFIG_PATH = _cfg_arg if _cfg_arg.is_absolute() else DATA_DIR / _cfg_arg
print(f"→ config: {CONFIG_PATH}")

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

# Comparison periods for early-vs-recent change maps (vision M1+M2).
# Defaults match vision.md spec; configs may override.
cmp_cfg = cfg.get("comparison", {}) or {}
EARLY_YEARS  = tuple(cmp_cfg.get("early_years",  [2017, 2018]))
RECENT_YEARS = tuple(cmp_cfg.get("recent_years", [2022, 2025]))
CMP_MONTHS   = tuple(cmp_cfg.get("months",       [5, 9]))

# Establishment status classification thresholds (vision M3). Per vision.md
# these are intentionally empirical and meant to be tuned after inspecting
# the first classified map.
cls_cfg = cfg.get("classification", {}) or {}
T_CHANGE      = float(cls_cfg.get("change_magnitude",     0.05))
T_NDVI_LOW    = float(cls_cfg.get("ndvi_low_threshold",   0.30))
T_BSI_HIGH    = float(cls_cfg.get("bsi_high_threshold",   0.20))
T_WATER       = float(cls_cfg.get("ndwi_water_threshold", 0.20))

# Reference event date for the vertical line on time-series plots.
# Set reference_date: null in the config to hide the line entirely.
REF_DATE  = tmp.get("reference_date",  "2022-04-01")
REF_LABEL = tmp.get("reference_label", "Reference T₀")

# Region label (used in all output filenames — needs to be defined before any
# block that might be skipped by --map-only).
region_name = sp.get("region_name", "AOI")

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
        gj_arg = Path(sp_cfg.get("geojson_path", sp_cfg.get("path")))
        p = gj_arg if gj_arg.is_absolute() else DATA_DIR / gj_arg
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
    ndwi = img.normalizedDifference(["B3",  "B8" ]).rename("NDWI")   # McFeeters: open water
    nbr  = img.normalizedDifference(["B8",  "B12"]).rename("NBR")    # fire damage
    bsi  = (img.select("B11").add(img.select("B4"))
               .subtract(img.select("B8").add(img.select("B2")))
               .divide(img.select("B11").add(img.select("B4"))
                          .add(img.select("B8")).add(img.select("B2")))
               ).rename("BSI")
    return img.addBands([ndvi, ndre, ndmi, bsi, ndwi, nbr])

SPECTRAL_BANDS = ["NDVI", "NDRE", "NDMI", "BSI", "NDWI", "NBR"]

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
# To replace the proxies above with the full neural-network model, set up
# gee_biophys.py (see --biophys-method flag) and uncomment:
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

# Multi-year seasonal composite (vision.md Milestone 1).
# Merges (year_start..year_end) × (month_start..month_end) S2 scenes, applies
# the same Cloud Score Plus mask + index/biophys functions as the monthly
# pipeline, then takes the median. Returned image is lazy — evaluated only
# when added to a Map or passed to getThumbURL/reduceRegion.
def make_seasonal_composite(year_start: int, year_end: int,
                            month_start: int, month_end: int) -> ee.Image:
    cols = []
    for y in range(year_start, year_end + 1):
        ws = f"{y}-{month_start:02d}-01"
        we = f"{y}-{month_end:02d}-01"   # exclusive end (filterDate convention)
        cols.append(build_s2_col(ws, we)
                      .map(add_spectral_indices)
                      .map(add_requested_biophys))
    merged = cols[0]
    for c in cols[1:]:
        merged = merged.merge(c)
    return merged.median().clip(aoi)

# Establishment status map (vision.md Milestone 3).
# Per-pixel categorical classification from the rule logic in vision.md
# section "Explicit rule logic and interpretation layer".
#
# Class values (priority high→low: Uncertain > Weak > Good > Moderate):
#   0 = Uncertain  (any input band masked over the pixel)
#   1 = Good       (NDVI ↗ AND BSI ↘ AND NDMI stable↗ AND not standing water)
#   2 = Moderate   (default: valid pixel that is neither Good nor Weak)
#   3 = Weak       (NDVI stayed low OR BSI stayed high)
def build_establishment_status(early: ee.Image, recent: ee.Image,
                               change: ee.Image) -> ee.Image:
    ndvi_increased = change.select("NDVI_change").gt( T_CHANGE)
    bsi_decreased  = change.select("BSI_change" ).lt(-T_CHANGE)
    ndmi_stable_up = change.select("NDMI_change").gt(-T_CHANGE)
    not_water      = recent.select("NDWI"       ).lt( T_WATER)
    good = ndvi_increased.And(bsi_decreased).And(ndmi_stable_up).And(not_water)

    ndvi_low      = recent.select("NDVI").lt(T_NDVI_LOW)
    bsi_high      = recent.select("BSI" ).gt(T_BSI_HIGH)
    weak = ndvi_low.Or(bsi_high)

    # Start with Moderate (2); promote Good (1); override with Weak (3) per
    # priority; finally mark Uncertain (0) wherever any input is masked.
    cls = ee.Image.constant(2).rename("establishment_status")
    cls = cls.where(good, 1)
    cls = cls.where(weak, 3)

    # Compose validity mask from all input bands used by the rules. .mask()
    # returns 1 where data is present, 0 where masked.
    valid = (recent.select("NDVI").mask()
             .multiply(recent.select("BSI" ).mask())
             .multiply(recent.select("NDWI").mask())
             .multiply(change.select("NDVI_change").mask())
             .multiply(change.select("BSI_change" ).mask())
             .multiply(change.select("NDMI_change").mask()))
    cls = cls.where(valid.eq(0), 0)

    # Re-mask outside the AOI so the categorical raster respects the polygon.
    return cls.clip(aoi).toUint8().rename("establishment_status")

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
if MAP_ONLY:
    print("  (—map-only: skipped, no time-series plots will be produced)")
    records: list[dict] = []
else:
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

# ── CSV export — time-series data ────────────────────────────────────────────
if records:
    import csv as _csv
    _ts_csv = SRC_DIR / f"timeseries_{region_name}.csv"
    _fields = ["date", "label"] + ALL_BANDS
    with open(_ts_csv, "w", newline="") as _f:
        _w = _csv.DictWriter(_f, fieldnames=_fields, extrasaction="ignore")
        _w.writeheader()
        _w.writerows(records)
    print(f"Saved → {_ts_csv}")

# ── Plots — spectral & biophysical time series ───────────────────────────────
# Skipped entirely in --map-only mode.
PLANTING = datetime.strptime(REF_DATE, "%Y-%m-%d") if REF_DATE else None

SPEC_COLORS  = {
    "NDVI": "#2d8a2d", "NDRE": "#7baf27", "NDMI": "#1e6eb5", "BSI": "#c47a1e",
    "NDWI": "#1ab1c4", "NBR":  "#a83232",
}
BIOPH_COLORS = {"laie": "#8b4513", "fapar": "#ff6600", "fcover": "#009900"}
BIOPH_LABELS = {"laie": "LAI-e [m²/m²]", "fapar": "FAPAR [0-1]", "fcover": "FCOVER [0-1]"}

def _plot_band(ax, band, color, ylabel, records, dates):
    vals = [r.get(band) for r in records]
    x = [d for d, v in zip(dates, vals) if v is not None]
    y = [v for v in vals if v is not None]
    ax.plot(x, y, marker="o", ms=3, lw=1.5, color=color, label=band)
    if PLANTING is not None:
        ax.axvline(PLANTING, color="red", ls="--", lw=1, label=REF_LABEL)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

if not MAP_ONLY:
    dates = [datetime.strptime(r["date"], "%Y-%m-%d") for r in records]
    cadence_str = f"{tmp['cadence']['type']} / {tmp['cadence'].get('interval','seasons')}"

    # Plot 1 — spectral indices (2 × 3 grid)
    fig1, axes1 = plt.subplots(2, 3, figsize=(16, 8), sharex=True)
    fig1.suptitle(f"{region_name} — Spectral indices ({cadence_str})", fontsize=12)
    for ax, band in zip(axes1.flat, SPECTRAL_BANDS):
        _plot_band(ax, band, SPEC_COLORS[band], band, records, dates)
    fig1.autofmt_xdate()
    plt.tight_layout()
    out1 = SRC_DIR / f"timeseries_spectral_{region_name}.png"
    fig1.savefig(out1, dpi=150)
    plt.show()
    print(f"Saved → {out1}")

    # Plot 2 — biophysical variables
    n = len(BIOPHYS_VARS)
    fig2, axes2 = plt.subplots(1, n, figsize=(5 * n, 4), sharex=True)
    if n == 1:
        axes2 = [axes2]
    fig2.suptitle(f"{region_name} — Biophysical variables / s2biophys proxies "
                  f"({cadence_str})", fontsize=11)
    for ax, band in zip(axes2, BIOPHYS_VARS):
        _plot_band(ax, band, BIOPH_COLORS.get(band, "grey"),
                   BIOPH_LABELS.get(band, band), records, dates)
    fig2.autofmt_xdate()
    plt.tight_layout()
    out2 = SRC_DIR / f"timeseries_biophys_{region_name}.png"
    fig2.savefig(out2, dpi=150)
    plt.show()
    print(f"Saved → {out2}")

# ── Raster visualisation — latest summer composite ───────────────────────────
RASTER_START, RASTER_END = "2023-06-01", "2023-08-31"
raster = make_composite(RASTER_START, RASTER_END)
bounds = aoi.bounds().getInfo()["coordinates"][0]
# ── Early vs Recent seasonal composites (vision M1) ──────────────────────
print(f"Building early composite (months {CMP_MONTHS[0]}-{CMP_MONTHS[1]-1}, "
      f"years {EARLY_YEARS[0]}-{EARLY_YEARS[1]})…")
early_composite  = make_seasonal_composite(EARLY_YEARS[0],  EARLY_YEARS[1],
                                           CMP_MONTHS[0],    CMP_MONTHS[1])
print(f"Building recent composite (months {CMP_MONTHS[0]}-{CMP_MONTHS[1]-1}, "
      f"years {RECENT_YEARS[0]}-{RECENT_YEARS[1]})…")
recent_composite = make_seasonal_composite(RECENT_YEARS[0], RECENT_YEARS[1],
                                           CMP_MONTHS[0],    CMP_MONTHS[1])

# ── Change image (vision M2) ─────────────────────────────────────
print("Computing change image (recent − early)…")
CHANGE_BANDS = [f"{b}_change" for b in SPECTRAL_BANDS]
change_img = (recent_composite.select(SPECTRAL_BANDS)
                              .subtract(early_composite.select(SPECTRAL_BANDS))
                              .rename(CHANGE_BANDS))

# Diagnostic: print min/max per change band. bestEffort avoids the
# "Image.reduceRegion: Too many pixels" error on large AOIs at 20 m.
try:
    _stats = change_img.reduceRegion(
        reducer   = ee.Reducer.minMax(),
        geometry  = aoi,
        scale     = int(exp.get("scale", 20)),
        maxPixels = int(exp.get("max_pixels", 1e9)),
        bestEffort= True,
    ).getInfo()
    for _b in CHANGE_BANDS:
        _mn, _mx = _stats.get(f"{_b}_min"), _stats.get(f"{_b}_max")
        if _mn is not None and _mx is not None:
            print(f"  {_b:14s}  min={_mn:+.3f}  max={_mx:+.3f}")
except ee.EEException as _e:
    print(f"  (change-band stats unavailable: {_e})")

# Diverging palettes for the change visualisation.
_DIVERGE_RG    = ["#7d2222", "white", "#1f5e1f"]      # red→white→green (+ = good)
_DIVERGE_GR    = ["#1f5e1f", "white", "#7d2222"]      # reversed (+ = bad, BSI)
_DIVERGE_BROWN = ["saddlebrown", "white", "steelblue"] # dry→neutral→wet

VIS_CHANGE = {
    "NDVI_change": {"min":-0.3, "max":0.3, "palette": _DIVERGE_RG},
    "NDRE_change": {"min":-0.3, "max":0.3, "palette": _DIVERGE_RG},
    "NDMI_change": {"min":-0.3, "max":0.3, "palette": _DIVERGE_RG},
    "BSI_change":  {"min":-0.3, "max":0.3, "palette": _DIVERGE_GR},
    "NDWI_change": {"min":-0.3, "max":0.3, "palette": _DIVERGE_BROWN},
    "NBR_change":  {"min":-0.3, "max":0.3, "palette": _DIVERGE_RG},
}

# ── Establishment status map (vision M3) ────────────────────────────────
print("Building establishment status map (M3)…")
establishment_status = build_establishment_status(
    early_composite, recent_composite, change_img,
)

# Class metadata — keep these in sync with build_establishment_status() values.
ESTAB_CLASSES = [
    (0, "Uncertain",              "#888888"),   # grey  — masked / no data
    (1, "Good establishment",     "#1f5e1f"),   # dark green
    (2, "Moderate establishment", "#a8d666"),   # light green
    (3, "Weak / bare-soil risk",  "#c47a1e"),   # orange
]
ESTAB_PALETTE = [c for _, _, c in sorted(ESTAB_CLASSES)]
VIS_ESTAB = {"min": 0, "max": 3, "palette": ESTAB_PALETTE}

# ── Landsat LST layer (optional) ─────────────────────────────────────────────
lst_products = None
_lst_cfg = cfg.get("landsat")
if _lst_cfg:
    from landsat_lst import build_lst_layer
    print("Building Landsat LST layers…")
    try:
        lst_products = build_lst_layer(aoi, _lst_cfg)
    except Exception as _e:
        print(f"  (LST build failed: {_e})")

# ── Per-layer legend content ──────────────────────────────────────────────────
def _safe_id(name: str) -> str:
    return "leg-" + "".join(c if c.isalnum() else "_" for c in name)

def _gradient_legend(title: str, palette: list, vmin: float, vmax: float,
                     unit: str = "", note: str = "") -> str:
    grad = ", ".join(palette)
    mid  = (vmin + vmax) / 2
    return (
        f'<div style="font-weight:600;margin-bottom:5px;">{title}</div>'
        f'<div style="height:12px;background:linear-gradient(to right,{grad});'
        f'border:1px solid #aaa;border-radius:2px;margin-bottom:4px;"></div>'
        f'<div style="display:flex;justify-content:space-between;font-size:10px;margin-bottom:2px;">'
        f'  <span>{vmin:.2g}</span><span>{mid:.2g}</span><span>{vmax:.2g}</span>'
        f'</div>'
        + (f'<div style="font-size:10px;color:#555;">{unit}</div>' if unit else '')
        + (f'<div style="font-size:10px;color:#777;font-style:italic;margin-top:3px;">{note}</div>' if note else '')
    )

def _categorical_legend(title: str, classes: list, note: str = "") -> str:
    rows = "".join(
        f'<div style="display:flex;align-items:center;margin:2px 0;">'
        f'<span style="display:inline-block;width:14px;height:14px;background:{c};'
        f'border:1px solid #555;margin-right:6px;flex-shrink:0;"></span>'
        f'<span style="font-size:11px;">{idx} — {lbl}</span></div>'
        for idx, lbl, c in classes
    )
    return (
        f'<div style="font-weight:600;margin-bottom:5px;">{title}</div>'
        + rows
        + (f'<div style="font-size:10px;color:#777;font-style:italic;margin-top:5px;">{note}</div>' if note else '')
    )

def _text_legend(title: str, note: str = "") -> str:
    return (
        f'<div style="font-weight:600;margin-bottom:4px;">{title}</div>'
        f'<div style="font-size:10px;color:#777;font-style:italic;">{note or "No radiometric scale."}</div>'
    )

_early_label  = f"Early RGB ({EARLY_YEARS[0]}–{EARLY_YEARS[1]})"
_recent_label = f"Recent RGB ({RECENT_YEARS[0]}–{RECENT_YEARS[1]})"

_LEGENDS: dict = {
    "RGB":    _text_legend("True-colour composite (B4/B3/B2)"),
    "NDVI":   _gradient_legend("NDVI", ["saddlebrown","khaki","limegreen","darkgreen"],
                               -0.1, 0.85,
                               note="Vegetation greenness — young forest ≈ 0.3–0.6"),
    "LAI-e":  _gradient_legend("LAI-e (proxy)", ["white","yellow","limegreen","darkgreen"],
                               0, 6, unit="m²/m²",
                               note="Linear proxy: 3.618 × NDVI − 0.118"),
    "FCOVER": _gradient_legend("FCOVER (proxy)", ["white","lightgreen","forestgreen"],
                               0, 1,
                               note="Green vegetation cover fraction"),
    "BSI":    _gradient_legend("BSI", ["darkgreen","white","orange","saddlebrown"],
                               -0.3, 0.4,
                               note="High positive = exposed bare soil"),
    "NDWI":   _gradient_legend("NDWI", ["lightyellow","powderblue","steelblue","navy"],
                               -0.3, 0.5,
                               note="Open water / flooding (McFeeters 1996)"),
    "NBR":    _gradient_legend("NBR", ["#7d2222","#cc6600","#ffe066","#a8d666","#1f5e1f"],
                               -0.3, 0.7,
                               note="Low = burned or stressed vegetation"),
    "Establishment status": _categorical_legend(
        "Establishment status", ESTAB_CLASSES,
        note="Greening signal only — not proof of tree survival"),
    _early_label:  _text_legend(_early_label,
                                f"Months {CMP_MONTHS[0]}–{CMP_MONTHS[1]-1} seasonal median"),
    _recent_label: _text_legend(_recent_label,
                                f"Months {CMP_MONTHS[0]}–{CMP_MONTHS[1]-1} seasonal median"),
    **{band: _gradient_legend(band, params["palette"], params["min"], params["max"],
                              note="Change = recent − early composite")
       for band, params in {
           "NDVI_change": {"min":-0.3,"max":0.3,"palette":["#7d2222","white","#1f5e1f"]},
           "NDRE_change": {"min":-0.3,"max":0.3,"palette":["#7d2222","white","#1f5e1f"]},
           "NDMI_change": {"min":-0.3,"max":0.3,"palette":["#7d2222","white","#1f5e1f"]},
           "BSI_change":  {"min":-0.3,"max":0.3,"palette":["#1f5e1f","white","#7d2222"]},
           "NDWI_change": {"min":-0.3,"max":0.3,"palette":["saddlebrown","white","steelblue"]},
           "NBR_change":  {"min":-0.3,"max":0.3,"palette":["#7d2222","white","#1f5e1f"]},
       }.items()
    },
    **(
        {
            "LST median (Landsat)": _gradient_legend(
                "LST median (Landsat)",
                ["#313695","#4575b4","#abd9e9","#ffffbf","#fdae61","#a50026"],
                15, 45, unit="°C",
                note=f"Target season {_lst_cfg['target_year']} — Landsat 8/9 C2 L2SP",
            ),
            "LST anomaly (Landsat)": _gradient_legend(
                "LST anomaly (Landsat)",
                ["#313695","#4575b4","#e0f3f8","#ffffff","#fee090","#a50026"],
                -5, 5, unit="°C",
                note="Target minus baseline median; blue = cooler than normal",
            ),
            "LST stress class (Landsat)": _categorical_legend(
                "Thermal stress class",
                [(0, "No valid data",          "#888888"),
                 (1, "Near-normal",            "#2ca25f"),
                 (2, "Moderate warm anomaly",  "#feb24c"),
                 (3, "Strong warm anomaly",    "#de2d26")],
            ),
        }
        if _lst_cfg else {}
    ),
}

# Diagnostic: per-class pixel count and percentage over the AOI.
try:
    _hist = establishment_status.reduceRegion(
        reducer   = ee.Reducer.frequencyHistogram(),
        geometry  = aoi,
        scale     = int(exp.get("scale", 20)),
        maxPixels = int(exp.get("max_pixels", 1e9)),
        bestEffort= True,
    ).getInfo()
    _counts = _hist.get("establishment_status") or {}
    # frequencyHistogram returns string-keyed counts; coerce + total
    _counts = {int(float(k)): int(v) for k, v in _counts.items()}
    _total  = sum(_counts.values()) or 1
    print(f"Establishment status histogram ({region_name}):")
    for _idx, _name, _ in ESTAB_CLASSES:
        _n = _counts.get(_idx, 0)
        _pct = 100.0 * _n / _total
        print(f"  {_name:26s}  {_n:>10,d} px  ({_pct:5.1f}%)")
except ee.EEException as _e:
    print(f"  (establishment histogram unavailable: {_e})")

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
    "NDWI":   {"bands": ["NDWI"],  "params": {"min":-0.3,"max":0.5,
               "palette":["lightyellow","powderblue","steelblue","navy"]}},
    "NBR":    {"bands": ["NBR"],   "params": {"min":-0.3,"max":0.7,
               "palette":["#7d2222","#cc6600","#ffe066","#a8d666","#1f5e1f"]}},
}
# Only keep layers whose bands exist in ALL_BANDS or raw S2
VIS_AVAILABLE = {k: v for k, v in VIS.items()
                 if all(b in ALL_BANDS + ["B4","B3","B2"] for b in v["bands"])}

try:
    # Use folium directly. geemap.Map (ipyleaflet) produces HTML that
    # depends on Jupyter widgets and fails outside a notebook
    # ("Class null not found in module @jupyter-widgets/base").
    # geemap.foliumap fails to import on the current xyzservices version.
    # Raw folium is reliable, self-contained, and ships with a working
    # LayerControl widget.
    import folium

    # Centre the map on the AOI bounds (reuse `bounds` already fetched above)
    _lons = [c[0] for c in bounds]
    _lats = [c[1] for c in bounds]
    _center = [(min(_lats) + max(_lats)) / 2, (min(_lons) + max(_lons)) / 2]
    fmap = folium.Map(location=_center, zoom_start=10, tiles="OpenStreetMap",
                      control_scale=True)

    _visible_on_load: list = []

    def _add_ee_layer(ee_image: ee.Image, params: dict, name: str,
                      show: bool = False) -> None:
        if show:
            _visible_on_load.append(name)
        map_id = ee_image.getMapId(params)
        folium.TileLayer(
            tiles=map_id["tile_fetcher"].url_format,
            attr="Google Earth Engine",
            name=name,
            overlay=True,
            control=True,
            show=show,
        ).add_to(fmap)

    # Single-period rasters — only RGB on by default; others toggled from layer panel
    for title, vc in VIS_AVAILABLE.items():
        _add_ee_layer(raster.select(vc["bands"]), vc["params"], title,
                      show=(title == "RGB"))
    # Establishment status (M3) — headline product, visible by default
    _add_ee_layer(establishment_status, VIS_ESTAB, "Establishment status",
                  show=True)
    # Early & recent seasonal composites (RGB) — hidden by default
    _RGB_PARAMS = {"min": 0, "max": 2800, "gamma": 1.4}
    _add_ee_layer(early_composite.select(["B4", "B3", "B2"]),  _RGB_PARAMS,
                  f"Early RGB ({EARLY_YEARS[0]}–{EARLY_YEARS[1]})",  show=False)
    _add_ee_layer(recent_composite.select(["B4", "B3", "B2"]), _RGB_PARAMS,
                  f"Recent RGB ({RECENT_YEARS[0]}–{RECENT_YEARS[1]})", show=False)
    # Per-index change layers — hidden by default; tick in the layer panel
    for _band, _params in VIS_CHANGE.items():
        _add_ee_layer(change_img.select(_band), _params, _band, show=False)

    # Landsat LST layers — hidden by default; only added when landsat: block present
    if lst_products:
        _LST_PALETTE  = ["#313695","#4575b4","#abd9e9","#ffffbf","#fdae61","#a50026"]
        _ANOM_PALETTE = ["#313695","#4575b4","#e0f3f8","#ffffff","#fee090","#a50026"]
        _STRESS_PALETTE = ["#888888","#2ca25f","#feb24c","#de2d26"]
        _add_ee_layer(lst_products["lst_median"],
                      {"min": 15, "max": 45, "palette": _LST_PALETTE},
                      "LST median (Landsat)", show=False)
        _add_ee_layer(lst_products["lst_anomaly"],
                      {"min": -5, "max": 5, "palette": _ANOM_PALETTE},
                      "LST anomaly (Landsat)", show=False)
        _add_ee_layer(lst_products["lst_stress"],
                      {"min": 0, "max": 3, "palette": _STRESS_PALETTE},
                      "LST stress class (Landsat)", show=False)

    # AOI outline as a true vector overlay (clearer than an EE tile mask)
    folium.GeoJson(
        aoi.getInfo(),
        name="AOI",
        style_function=lambda _f: {"color": "red", "weight": 2, "fillOpacity": 0},
    ).add_to(fmap)

    # ── Dynamic per-layer legend ───────────────────────────────────────────────
    # One div per layer; shown/hidden via JS overlayadd / overlayremove events.
    _inner_divs = "\n".join(
        f'<div id="{_safe_id(n)}" '
        f'style="display:{"block" if n in _visible_on_load else "none"};">'
        f'{html}</div>'
        for n, html in _LEGENDS.items()
    )
    _any_visible_init = "none" if _visible_on_load else "block"
    _legend_outer = (
        '<div id="dyn-legend" style="'
        'position:fixed;bottom:30px;right:30px;z-index:9999;'
        'background:rgba(255,255,255,0.95);padding:10px 12px;'
        'border:1px solid #888;border-radius:4px;'
        'font:12px/1.3 system-ui,sans-serif;max-width:260px;'
        'box-shadow:0 2px 6px rgba(0,0,0,0.2);">'
        + _inner_divs
        + f'<div id="dyn-legend-empty" style="display:{_any_visible_init};font-size:11px;color:#888;">'
          f'Toggle a layer to see its legend.</div>'
        + '</div>'
    )
    fmap.get_root().html.add_child(folium.Element(_legend_outer))

    _id_map_js  = json.dumps({n: _safe_id(n) for n in _LEGENDS})
    _map_var    = fmap.get_name()
    # Raw JS only — no <script> wrapper. folium embeds this inside its own
    # <script> block, so adding wrapper tags would create nested <script> tags
    # which cause the browser to close the outer block early and never run
    # the L.map() initialisation. The poll loop handles the case where this
    # code executes before the map variable has been assigned.
    _dyn_script = f"""
(function poll() {{
  var m = window['{_map_var}'];
  if (!m) {{ setTimeout(poll, 50); return; }}

  var ID_MAP = {_id_map_js};

  function refresh() {{
    var any = Object.values(ID_MAP).some(function(id) {{
      var el = document.getElementById(id);
      return el && el.style.display !== 'none';
    }});
    var emp = document.getElementById('dyn-legend-empty');
    if (emp) emp.style.display = any ? 'none' : 'block';
  }}

  function setLayer(name, visible) {{
    var id = ID_MAP[name];
    if (!id) return;
    var el = document.getElementById(id);
    if (el) el.style.display = visible ? 'block' : 'none';
    refresh();
  }}

  m.on('overlayadd',    function(e) {{ setLayer(e.name, true);  }});
  m.on('overlayremove', function(e) {{ setLayer(e.name, false); }});
  refresh();
}})();
"""
    fmap.get_root().script.add_child(folium.Element(_dyn_script))

    # ── Interactive time-series side panel ────────────────────────────────────
    # Chart.js + date adapter + annotation plugin (CDN, loaded in <head>)
    fmap.get_root().header.add_child(folium.Element(
        '<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>\n'
        '<script src="https://cdn.jsdelivr.net/npm/chartjs-adapter-date-fns@3.0.0/dist/chartjs-adapter-date-fns.bundle.min.js"></script>\n'
        '<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@3.0.1/dist/chartjs-plugin-annotation.min.js"></script>'
    ))

    # Embed time-series records as JSON (only the bands we plot)
    _ts_records_clean = [
        {k: r.get(k) for k in ["date", "label"] + ALL_BANDS}
        for r in records
    ]
    _ts_records_json   = json.dumps(_ts_records_clean)
    _ts_layer_vars     = {
        "RGB":                  ["NDVI"],
        "NDVI":                 ["NDVI"],
        "LAI-e":                ["laie"],
        "FCOVER":               ["fcover"],
        "BSI":                  ["BSI"],
        "NDWI":                 ["NDWI"],
        "NBR":                  ["NBR"],
        "Establishment status": ["NDVI", "BSI", "NDMI"],
        _early_label:           ["NDVI"],
        _recent_label:          ["NDVI"],
        "NDVI_change":          ["NDVI"],
        "NDRE_change":          ["NDRE"],
        "NDMI_change":          ["NDMI"],
        "BSI_change":           ["BSI"],
        "NDWI_change":          ["NDWI"],
        "NBR_change":           ["NBR"],
    }
    _ts_layer_vars_json = json.dumps(_ts_layer_vars)
    _ts_colors_json     = json.dumps({**SPEC_COLORS, **BIOPH_COLORS})
    _ts_ref_date        = json.dumps(REF_DATE  or "")
    _ts_ref_label       = json.dumps(REF_LABEL or "")
    _ts_init_json       = json.dumps(_visible_on_load)

    _ts_panel_html = (
        '<div id="ts-panel" style="'
        'position:fixed;left:50px;top:60px;z-index:9998;'
        'background:rgba(255,255,255,0.97);'
        'padding:12px 14px;border:1px solid #888;border-radius:4px;'
        'box-shadow:0 2px 8px rgba(0,0,0,0.25);width:340px;'
        'font:11px/1.4 system-ui,sans-serif;display:none;">'
        '<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;">'
        '<span id="ts-title" style="font-weight:600;font-size:12px;max-width:280px;'
        'overflow:hidden;text-overflow:ellipsis;white-space:nowrap;"></span>'
        '<button id="ts-close" style="background:none;border:none;cursor:pointer;'
        'font-size:18px;line-height:1;color:#555;flex-shrink:0;">×</button>'
        '</div>'
        '<div id="ts-no-data" style="display:none;color:#888;font-size:11px;font-style:italic;">'
        'No time-series data available (re-run without --map-only).</div>'
        '<div style="position:relative;height:200px;">'
        '<canvas id="ts-chart"></canvas>'
        '</div>'
        '</div>'
    )
    fmap.get_root().html.add_child(folium.Element(_ts_panel_html))

    _ts_script = f"""
(function waitForChart() {{
  if (typeof Chart === 'undefined') {{ setTimeout(waitForChart, 50); return; }}

  var TS_DATA    = {_ts_records_json};
  var LAYER_VARS = {_ts_layer_vars_json};
  var COLORS     = {_ts_colors_json};
  var REF_DATE   = {_ts_ref_date};
  var REF_LABEL  = {_ts_ref_label};
  var INIT_LAYERS = {_ts_init_json};

  var panel   = document.getElementById('ts-panel');
  var noData  = document.getElementById('ts-no-data');
  var chartWrap = panel.querySelector('div[style*="height:200px"]');
  var tsChart = null;

  document.getElementById('ts-close').addEventListener('click', function() {{
    panel.style.display = 'none';
  }});

  function showPanel(layerName) {{
    var vars = LAYER_VARS[layerName];
    if (!vars || vars.length === 0) return;

    document.getElementById('ts-title').textContent = layerName;
    panel.style.display = 'block';

    if (TS_DATA.length === 0) {{
      noData.style.display = 'block';
      if (chartWrap) chartWrap.style.display = 'none';
      return;
    }}
    noData.style.display = 'none';
    if (chartWrap) chartWrap.style.display = 'block';

    var datasets = vars.map(function(v) {{
      var pts = TS_DATA
        .filter(function(r) {{ return r[v] !== null && r[v] !== undefined; }})
        .map(function(r)   {{ return {{ x: r.date, y: r[v] }}; }});
      return {{
        label: v,
        data: pts,
        borderColor: COLORS[v] || '#4a9',
        backgroundColor: (COLORS[v] || '#4a9') + '22',
        borderWidth: 1.5,
        pointRadius: 2.5,
        pointHoverRadius: 4,
        tension: 0.2,
        fill: false,
      }};
    }});

    var annotations = {{}};
    if (REF_DATE) {{
      annotations['ref'] = {{
        type: 'line',
        xMin: REF_DATE, xMax: REF_DATE,
        borderColor: 'rgba(210,40,40,0.8)',
        borderWidth: 1.5,
        borderDash: [6, 4],
        label: {{
          display: true,
          content: REF_LABEL,
          position: 'start',
          backgroundColor: 'rgba(210,40,40,0.1)',
          color: '#c00',
          font: {{ size: 9 }},
          padding: 3,
        }},
      }};
    }}

    if (tsChart) tsChart.destroy();
    var ctx = document.getElementById('ts-chart').getContext('2d');
    tsChart = new Chart(ctx, {{
      type: 'line',
      data: {{ datasets: datasets }},
      options: {{
        responsive: true,
        maintainAspectRatio: false,
        interaction: {{ mode: 'index', intersect: false }},
        scales: {{
          x: {{
            type: 'time',
            time: {{ unit: 'year', tooltipFormat: 'yyyy-MM-dd' }},
            ticks: {{ maxTicksLimit: 6, font: {{ size: 9 }} }},
            grid: {{ color: 'rgba(0,0,0,0.06)' }},
          }},
          y: {{
            ticks: {{ maxTicksLimit: 5, font: {{ size: 9 }} }},
            grid: {{ color: 'rgba(0,0,0,0.06)' }},
          }},
        }},
        plugins: {{
          legend: {{
            display: vars.length > 1,
            labels: {{ font: {{ size: 9 }}, boxWidth: 12 }},
          }},
          tooltip: {{
            callbacks: {{
              label: function(ctx) {{
                var v = ctx.parsed.y;
                return ctx.dataset.label + ': ' + (v !== null ? v.toFixed(3) : 'N/A');
              }},
            }},
          }},
          annotation: {{ annotations: annotations }},
        }},
      }},
    }});
  }}

  (function pollMap() {{
    var m = window['{_map_var}'];
    if (!m) {{ setTimeout(pollMap, 50); return; }}

    m.on('overlayadd', function(e) {{
      showPanel(e.name);
    }});
    m.on('overlayremove', function(e) {{
      var title = document.getElementById('ts-title');
      if (title && title.textContent === e.name) {{
        panel.style.display = 'none';
      }}
    }});

    for (var i = 0; i < INIT_LAYERS.length; i++) {{
      if (LAYER_VARS[INIT_LAYERS[i]]) {{
        showPanel(INIT_LAYERS[i]);
        break;
      }}
    }}
  }})();
}})();
"""
    fmap.get_root().script.add_child(folium.Element(_ts_script))

    folium.LayerControl(collapsed=False).add_to(fmap)
    out_map = SRC_DIR / f"map_{region_name}.html"
    fmap.save(str(out_map))
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
    out3 = SRC_DIR / f"raster_{region_name}.png"
    fig3.savefig(out3, dpi=150)
    plt.show()
    print(f"Saved → {out3}")

    # Change-band thumbnails (vision M2)
    fig4, axes4 = plt.subplots(2, 3, figsize=(12, 8))
    fig4.suptitle(
        f"{region_name} — Change layers "
        f"(recent {RECENT_YEARS[0]}–{RECENT_YEARS[1]} − "
        f"early {EARLY_YEARS[0]}–{EARLY_YEARS[1]}, months "
        f"{CMP_MONTHS[0]}–{CMP_MONTHS[1]-1})",
        fontsize=11,
    )
    for ax, (band, params) in zip(axes4.flat, VIS_CHANGE.items()):
        url = change_img.select(band).getThumbURL(
            {**params, "region": bounds, "dimensions": 512, "format": "png"}
        )
        with urllib.request.urlopen(url) as resp:
            img = Image.open(io.BytesIO(resp.read()))
        ax.imshow(img)
        ax.set_title(band, fontsize=9)
        ax.axis("off")
    plt.tight_layout()
    out4 = SRC_DIR / f"change_{region_name}.png"
    fig4.savefig(out4, dpi=150)
    plt.show()
    print(f"Saved → {out4}")
