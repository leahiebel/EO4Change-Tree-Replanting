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

with open(CONFIG_PATH, encoding="utf-8") as f:
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
PREFIRE_YEARS  = tuple(cmp_cfg.get("prefire_years",  [2017, 2017]))
PREFIRE_MONTHS = tuple(cmp_cfg.get("prefire_months", [4, 6]))   # Apr–May: before June 17 fire
POSTFIRE_YEARS = tuple(cmp_cfg.get("postfire_years",
                        cmp_cfg.get("early_years",   [2017, 2017])))  # backward compat
RECENT_YEARS   = tuple(cmp_cfg.get("recent_years",  [2022, 2025]))
CMP_MONTHS     = tuple(cmp_cfg.get("months",        [5, 9]))

# Forest mask — Dynamic World pre-fire tree probability.
# Pixels outside this mask are excluded from recovery classification so that
# cropland, urban areas, and water bodies do not pollute the statistics.
fm_cfg  = cfg.get("forest_mask", {}) or {}
FM_ENABLED   = bool(fm_cfg.get("enabled", False))
FM_YEAR      = int(fm_cfg.get("year", 2016))
FM_MONTHS    = tuple(fm_cfg.get("months", [5, 9]))
FM_THRESHOLD = float(fm_cfg.get("tree_prob_threshold", 0.40))

# Recovery class and guard-flag thresholds.
# The composite index is the rule-based recovery_class map; each threshold
# has a direct ecological interpretation and is validated by OAT sensitivity
# analysis in calibrate_thresholds.py.
cls_cfg = cfg.get("classification", {}) or {}
# RRI thresholds (primary recovery metric)
T_RRI_GOOD      = float(cls_cfg.get("rri_good_threshold",        0.80))
T_RRI_LOW       = float(cls_cfg.get("rri_low_threshold",         0.30))
T_RRI_MIN_DENOM = float(cls_cfg.get("rri_min_denominator",       0.10))
T_BURN          = float(cls_cfg.get("burn_severity_threshold",   0.10))
# Secondary canopy-quality check (Class 1) and guard-flag thresholds
T_NDRE_GOOD  = float(cls_cfg.get("ndre_good_threshold",          0.18))
T_DNDVI_GOOD = float(cls_cfg.get("ndvi_change_good_threshold",   0.09))
T_BSI_HIGH   = float(cls_cfg.get("bsi_high_threshold",           0.09))
T_WATER      = float(cls_cfg.get("ndwi_water_threshold",         0.15))
T_DISTURB    = float(cls_cfg.get("disturbance_threshold",        0.09))

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
# To replace the proxies above with the full neural-network model, uncomment:
#
#   import gee_biophys
#   results = gee_biophys.retrieve(cfg)
#
# cfg["variables"]["variable"] is already a list: ["laie", "fapar", "fcover"]
# gee_biophys.retrieve() reads this and submits one export task per composite.

# ── Composite builder ─────────────────────────────────────────────────────────
ALL_BANDS = SPECTRAL_BANDS + BIOPHYS_VARS

# Only keep the bands we actually need in final composites.
# This avoids Earth Engine errors caused by inconsistent Sentinel-2 QA/mask bands.
RGB_BANDS = ["B4", "B3", "B2"]
COMPOSITE_BANDS = RGB_BANDS + ALL_BANDS

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
        we = f"{y}-{month_end:02d}-01"

        cols.append(build_s2_col(ws, we)
                      .map(add_spectral_indices)
                      .map(add_requested_biophys)
                      .map(lambda img: img.select(COMPOSITE_BANDS)))

    merged = cols[0]
    for c in cols[1:]:
        merged = merged.merge(c)

    return merged.median().select(COMPOSITE_BANDS).clip(aoi)

def build_forest_mask() -> ee.Image | None:
    """Return a binary mask (1 = was forest pre-fire) from Dynamic World tree
    probability, or None when forest masking is disabled in the config.
    Uses a summer median of the DW 'trees' band for FM_YEAR.
    """
    if not FM_ENABLED:
        return None
    ws = f"{FM_YEAR}-{FM_MONTHS[0]:02d}-01"
    we = f"{FM_YEAR}-{FM_MONTHS[1]:02d}-01"
    dw = (ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
          .filterBounds(aoi)
          .filterDate(ws, we)
          .select("trees")
          .median()
          .clip(aoi))
    return dw.gte(FM_THRESHOLD).rename("forest_mask")


def _recovery_valid_mask(recent: ee.Image, change: ee.Image) -> ee.Image:
    return (recent.select("NDVI").mask()
            .multiply(recent.select("NDRE").mask())
            .multiply(recent.select("BSI").mask())
            .multiply(recent.select("NDWI").mask())
            .multiply(change.select("NDVI_change").mask())
            .multiply(change.select("NBR_change").mask()))


def compute_rri(recent: ee.Image, postfire: ee.Image,
                prefire: ee.Image) -> ee.Image:
    """Relative Recovery Indicator (NDVI-based):
        RRI = (recent_NDVI − postfire_NDVI) / (prefire_NDVI − postfire_NDVI)
    RRI ≈ 1.0 → fully back to pre-fire level.
    RRI ≈ 0   → still at post-fire (damage) level.
    RRI < 0   → worse than post-fire (new disturbance or severe stress).
    Masked where (prefire − postfire) < T_RRI_MIN_DENOM (fire signal too
    small to compute a meaningful ratio).
    """
    ndvi_recent   = recent.select("NDVI")
    ndvi_postfire = postfire.select("NDVI")
    ndvi_prefire  = prefire.select("NDVI")
    denom = ndvi_prefire.subtract(ndvi_postfire)
    rri = (ndvi_recent.subtract(ndvi_postfire)
                      .divide(denom.max(ee.Image.constant(T_RRI_MIN_DENOM)))
                      .clamp(-1, 2)
                      .rename("RRI")
                      .updateMask(denom.gte(T_RRI_MIN_DENOM)))
    return rri.clip(aoi)


# Recovery class values:
#   0 = Uncertain / water / insufficient observations
#   1 = Recovering well   (burned + RRI ≥ T_RRI_GOOD + NDRE ≥ T_NDRE_GOOD)
#   2 = Recovering, but weak  (burned + T_RRI_LOW ≤ RRI < T_RRI_GOOD)
#   3 = Not recovering / failed  (burned + RRI < T_RRI_LOW)
#   4 = Stable / unburned forest  (dNBR < T_BURN; forest not significantly fire-affected)
def build_recovery_class(recent: ee.Image, postfire: ee.Image,
                         prefire: ee.Image, change: ee.Image) -> ee.Image:
    rri = compute_rri(recent, postfire, prefire)

    # Fire-affected pixels: dNBR ≥ T_BURN (USGS low-severity floor).
    # Pixels below this were not meaningfully burned → class 4.
    dnbr  = prefire.select("NBR").subtract(postfire.select("NBR"))
    burned = dnbr.gte(T_BURN)

    ndre_good = recent.select("NDRE").gte(T_NDRE_GOOD)
    water     = recent.select("NDWI").gt(T_WATER)
    valid     = _recovery_valid_mask(recent, change)

    good   = rri.gte(T_RRI_GOOD).And(ndre_good).And(water.Not()).And(burned)
    failed = rri.lt(T_RRI_LOW).And(water.Not()).And(burned)

    cls = ee.Image.constant(2).rename("recovery_class")  # default: burned but undecided
    cls = cls.where(good, 1)
    cls = cls.where(failed, 3)
    # Class 4: valid data but pixel not significantly fire-affected
    cls = cls.where(valid.eq(1).And(water.Not()).And(burned.Not()), 4)
    # Class 0: water or insufficient observations (overrides all)
    cls = cls.where(valid.eq(0).Or(water), 0)
    return cls.clip(aoi).toUint8().rename("recovery_class")


# Guard flags make the limitations explicit instead of burying them inside the
# class. Class values:
#   0 = Clear
#   1 = Possible water / non-vegetation pixel
#   2 = Possible disturbance
#   3 = Mixed signal (green-up without strong canopy/soil-closure evidence)
#   4 = Insufficient observations
def build_guard_flags(recent: ee.Image, postfire: ee.Image,
                      prefire: ee.Image, change: ee.Image) -> ee.Image:
    """
    Uncertainty flags for the burned-area recovery assessment.

    0 = Clear
    1 = Water / invalid observations
    2 = Possible later disturbance inside burned area
    3 = Mixed recovery signal inside burned area
    4 = Outside burn recovery assessment
    """

    # Same burn mask as recovery_class
    dnbr = prefire.select("NBR").subtract(postfire.select("NBR"))
    burned = dnbr.gte(T_BURN)

    water = recent.select("NDWI").gt(T_WATER)
    valid = _recovery_valid_mask(recent, change)

    # Only meaningful inside the burned assessment area
    disturbance = (
        change.select("NBR_change")
        .lt(-T_DISTURB)
        .And(burned)
        .And(water.Not())
        .And(valid.eq(1))
    )

    mixed_signal = (
        change.select("NDVI_change").gte(T_DNDVI_GOOD)
        .And(
            recent.select("NDRE").lt(T_NDRE_GOOD)
            .Or(recent.select("BSI").gt(T_BSI_HIGH))
        )
        .And(burned)
        .And(water.Not())
        .And(valid.eq(1))
    )

    out = ee.Image.constant(0).rename("guard_flags")

    # 4 = outside burned recovery assessment
    out = out.where(valid.eq(1).And(water.Not()).And(burned.Not()), 4)

    # 2 and 3 only inside burned area
    out = out.where(disturbance, 2)
    out = out.where(mixed_signal.And(disturbance.Not()), 3)

    # 1 overrides everything: invalid / water
    out = out.where(valid.eq(0).Or(water), 1)

    return out.clip(aoi).toUint8().rename("guard_flags")

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
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

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
# ── Pre-fire, post-fire, and recent seasonal composites ──────────────────────
print(f"Building pre-fire composite (months {PREFIRE_MONTHS[0]}-{PREFIRE_MONTHS[1]-1}, "
      f"year {PREFIRE_YEARS[0]})\u2026")
prefire_composite  = make_seasonal_composite(PREFIRE_YEARS[0],  PREFIRE_YEARS[1],
                                             PREFIRE_MONTHS[0],  PREFIRE_MONTHS[1])
print(f"Building post-fire composite (months {CMP_MONTHS[0]}-{CMP_MONTHS[1]-1}, "
      f"years {POSTFIRE_YEARS[0]}-{POSTFIRE_YEARS[1]})…")
postfire_composite = make_seasonal_composite(POSTFIRE_YEARS[0], POSTFIRE_YEARS[1],
                                             CMP_MONTHS[0],      CMP_MONTHS[1])
print(f"Building recent composite (months {CMP_MONTHS[0]}-{CMP_MONTHS[1]-1}, "
      f"years {RECENT_YEARS[0]}-{RECENT_YEARS[1]})…")
recent_composite = make_seasonal_composite(RECENT_YEARS[0], RECENT_YEARS[1],
                                           CMP_MONTHS[0],    CMP_MONTHS[1])

# ── Change image: recent vs. post-fire damage baseline (vision M2) ────────────
# Using 2017 (fire year only) as damage baseline excludes 2018 replanting
# contamination. The pre-fire composite is available separately for context.
print("Computing change image (recent − post-fire)…")
CHANGE_BANDS = [f"{b}_change" for b in SPECTRAL_BANDS]
change_img = (recent_composite.select(SPECTRAL_BANDS)
                              .subtract(postfire_composite.select(SPECTRAL_BANDS))
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

# ── Forest mask (Dynamic World pre-fire tree probability) ───────────────
if FM_ENABLED:
    print(f"Building forest mask (Dynamic World {FM_YEAR}, tree prob >= {FM_THRESHOLD})…")
forest_mask_img = build_forest_mask()

# ── Simplified recovery products ─────────────────────────────────────────
print("Computing RRI…")
rri_img = compute_rri(recent_composite, postfire_composite, prefire_composite)

print("Building recovery class…")
recovery_class = build_recovery_class(recent_composite, postfire_composite,
                                      prefire_composite, change_img)

print("Building guard flags…")
guard_flags = build_guard_flags(
    recent_composite,
    postfire_composite,
    prefire_composite,
    change_img
)
# Apply pre-fire forest mask: pixels outside the mask are excluded.
if forest_mask_img is not None:
    recovery_class = recovery_class.updateMask(forest_mask_img)
    guard_flags    = guard_flags.updateMask(forest_mask_img)

# ── Sensitivity analysis for key thresholds ──────────────────────────────────
# One-at-a-time sensitivity analysis:
#   1. vary burn_severity_threshold while keeping RRI thresholds fixed
#   2. vary rri_good_threshold while keeping burn + low RRI fixed
#   3. vary rri_low_threshold while keeping burn + good RRI fixed
#
# Output:
#   - CSV table with class counts and percentages
#   - one plot per tested threshold

sens_cfg = cfg.get("sensitivity_analysis", {}) or {}
SENS_ENABLED = bool(sens_cfg.get("enabled", False))

SENS_BURN_VALUES = sens_cfg.get(
    "burn_severity_threshold_values",
    [0.05, 0.075, T_BURN, 0.125, 0.15],
)

SENS_RRI_GOOD_VALUES = sens_cfg.get(
    "rri_good_threshold_values",
    [0.70, 0.75, T_RRI_GOOD, 0.85, 0.90],
)

SENS_RRI_LOW_VALUES = sens_cfg.get(
    "rri_low_threshold_values",
    [0.20, 0.25, T_RRI_LOW, 0.35, 0.40],
)


def build_recovery_class_sensitivity(
    recent: ee.Image,
    postfire: ee.Image,
    prefire: ee.Image,
    change: ee.Image,
    t_burn: float,
    t_rri_good: float,
    t_rri_low: float,
) -> ee.Image:
    """
    Same logic as build_recovery_class(), but with threshold values passed
    explicitly. This avoids changing the global config thresholds.
    """
    rri = compute_rri(recent, postfire, prefire)

    dnbr = prefire.select("NBR").subtract(postfire.select("NBR"))
    burned = dnbr.gte(t_burn)

    ndre_good = recent.select("NDRE").gte(T_NDRE_GOOD)
    water = recent.select("NDWI").gt(T_WATER)
    valid = _recovery_valid_mask(recent, change)

    good = (
        rri.gte(t_rri_good)
        .And(ndre_good)
        .And(water.Not())
        .And(burned)
    )

    failed = (
        rri.lt(t_rri_low)
        .And(water.Not())
        .And(burned)
    )

    cls = ee.Image.constant(2).rename("recovery_class")

    # 1 = recovering well
    cls = cls.where(good, 1)

    # 3 = not recovering / failed
    cls = cls.where(failed, 3)

    # 4 = outside burn recovery assessment
    cls = cls.where(
        valid.eq(1)
        .And(water.Not())
        .And(burned.Not()),
        4,
    )

    # 0 = water / invalid / insufficient observations
    cls = cls.where(valid.eq(0).Or(water), 0)

    return cls.clip(aoi).toUint8().rename("recovery_class")


def summarize_recovery_class(
    cls: ee.Image,
    parameter_name: str,
    parameter_value: float,
    t_burn: float,
    t_rri_good: float,
    t_rri_low: float,
) -> dict:
    """
    Count pixels in each recovery class and compute percentages.

    We report percentages in two ways:
      - over all valid pixels inside the forest mask
      - over the burned assessment area only: classes 1 + 2 + 3
    """
    if forest_mask_img is not None:
        cls = cls.updateMask(forest_mask_img)

    hist = cls.reduceRegion(
        reducer=ee.Reducer.frequencyHistogram(),
        geometry=aoi,
        scale=int(exp.get("scale", 20)),
        maxPixels=int(exp.get("max_pixels", 1e9)),
        bestEffort=True,
    ).getInfo()

    counts_raw = hist.get("recovery_class") or {}
    counts = {int(float(k)): int(v) for k, v in counts_raw.items()}

    c0 = counts.get(0, 0)
    c1 = counts.get(1, 0)
    c2 = counts.get(2, 0)
    c3 = counts.get(3, 0)
    c4 = counts.get(4, 0)

    total_px = c0 + c1 + c2 + c3 + c4
    burned_px = c1 + c2 + c3

    if total_px == 0:
        total_px = 1

    if burned_px == 0:
        burned_px = 1

    row = {
        "parameter": parameter_name,
        "value": parameter_value,

        "burn_severity_threshold": t_burn,
        "rri_good_threshold": t_rri_good,
        "rri_low_threshold": t_rri_low,

        "class0_uncertain_px": c0,
        "class1_recovering_well_px": c1,
        "class2_recovering_weak_px": c2,
        "class3_failed_px": c3,
        "class4_outside_assessment_px": c4,

        "total_forest_mask_px": total_px,
        "burned_assessment_px": burned_px,

        "class0_uncertain_pct_total": 100.0 * c0 / total_px,
        "class1_recovering_well_pct_total": 100.0 * c1 / total_px,
        "class2_recovering_weak_pct_total": 100.0 * c2 / total_px,
        "class3_failed_pct_total": 100.0 * c3 / total_px,
        "class4_outside_assessment_pct_total": 100.0 * c4 / total_px,

        "recovering_well_pct_burned": 100.0 * c1 / burned_px,
        "recovering_weak_pct_burned": 100.0 * c2 / burned_px,
        "failed_pct_burned": 100.0 * c3 / burned_px,
    }

    return row


def run_threshold_sensitivity_analysis() -> None:
    import pandas as pd
    import matplotlib.pyplot as plt

    rows = []

    print("Running threshold sensitivity analysis...")

    # Baseline row
    baseline_cls = build_recovery_class_sensitivity(
        recent=recent_composite,
        postfire=postfire_composite,
        prefire=prefire_composite,
        change=change_img,
        t_burn=T_BURN,
        t_rri_good=T_RRI_GOOD,
        t_rri_low=T_RRI_LOW,
    )

    rows.append(
        summarize_recovery_class(
            cls=baseline_cls,
            parameter_name="baseline",
            parameter_value=0.0,
            t_burn=T_BURN,
            t_rri_good=T_RRI_GOOD,
            t_rri_low=T_RRI_LOW,
        )
    )

    # 1. Sensitivity to burn severity threshold
    for value in SENS_BURN_VALUES:
        value = float(value)

        print(f"  testing burn_severity_threshold = {value}")

        cls = build_recovery_class_sensitivity(
            recent=recent_composite,
            postfire=postfire_composite,
            prefire=prefire_composite,
            change=change_img,
            t_burn=value,
            t_rri_good=T_RRI_GOOD,
            t_rri_low=T_RRI_LOW,
        )

        rows.append(
            summarize_recovery_class(
                cls=cls,
                parameter_name="burn_severity_threshold",
                parameter_value=value,
                t_burn=value,
                t_rri_good=T_RRI_GOOD,
                t_rri_low=T_RRI_LOW,
            )
        )

    # 2. Sensitivity to RRI good threshold
    for value in SENS_RRI_GOOD_VALUES:
        value = float(value)

        print(f"  testing rri_good_threshold = {value}")

        cls = build_recovery_class_sensitivity(
            recent=recent_composite,
            postfire=postfire_composite,
            prefire=prefire_composite,
            change=change_img,
            t_burn=T_BURN,
            t_rri_good=value,
            t_rri_low=T_RRI_LOW,
        )

        rows.append(
            summarize_recovery_class(
                cls=cls,
                parameter_name="rri_good_threshold",
                parameter_value=value,
                t_burn=T_BURN,
                t_rri_good=value,
                t_rri_low=T_RRI_LOW,
            )
        )

    # 3. Sensitivity to RRI low / failed threshold
    for value in SENS_RRI_LOW_VALUES:
        value = float(value)

        print(f"  testing rri_low_threshold = {value}")

        cls = build_recovery_class_sensitivity(
            recent=recent_composite,
            postfire=postfire_composite,
            prefire=prefire_composite,
            change=change_img,
            t_burn=T_BURN,
            t_rri_good=T_RRI_GOOD,
            t_rri_low=value,
        )

        rows.append(
            summarize_recovery_class(
                cls=cls,
                parameter_name="rri_low_threshold",
                parameter_value=value,
                t_burn=T_BURN,
                t_rri_good=T_RRI_GOOD,
                t_rri_low=value,
            )
        )

    df = pd.DataFrame(rows)

    # Add deltas relative to the baseline result.
    baseline = df[df["parameter"] == "baseline"].iloc[0]

    for col in [
        "recovering_well_pct_burned",
        "recovering_weak_pct_burned",
        "failed_pct_burned",
        "class4_outside_assessment_pct_total",
        "burned_assessment_px",
    ]:
        df[f"delta_{col}"] = df[col] - baseline[col]

    out_csv = SRC_DIR / f"sensitivity_thresholds_{region_name}.csv"
    df.to_csv(out_csv, index=False)

    print(f"Sensitivity CSV → {out_csv}")

    # Plot one figure per parameter.
    for parameter_name in [
        "burn_severity_threshold",
        "rri_good_threshold",
        "rri_low_threshold",
    ]:
        sub = df[df["parameter"] == parameter_name].copy()

        if sub.empty:
            continue

        sub = sub.sort_values("value")

        fig, ax = plt.subplots(figsize=(8, 5))

        ax.plot(
            sub["value"],
            sub["recovering_well_pct_burned"],
            marker="o",
            label="Recovering well (% of burned assessment)",
        )

        ax.plot(
            sub["value"],
            sub["recovering_weak_pct_burned"],
            marker="o",
            label="Recovering weak (% of burned assessment)",
        )

        ax.plot(
            sub["value"],
            sub["failed_pct_burned"],
            marker="o",
            label="Failed (% of burned assessment)",
        )

        # Class 4 is especially relevant for burn threshold sensitivity.
        if parameter_name == "burn_severity_threshold":
            ax.plot(
                sub["value"],
                sub["class4_outside_assessment_pct_total"],
                marker="o",
                linestyle="--",
                label="Outside assessment (% of total forest mask)",
            )

        ax.axvline(
            {
                "burn_severity_threshold": T_BURN,
                "rri_good_threshold": T_RRI_GOOD,
                "rri_low_threshold": T_RRI_LOW,
            }[parameter_name],
            linestyle=":",
            linewidth=1.5,
            label="Baseline value",
        )

        ax.set_xlabel(parameter_name)
        ax.set_ylabel("Pixel percentage")
        ax.set_title(f"Sensitivity analysis: {parameter_name}")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)

        plt.tight_layout()

        out_png = SRC_DIR / f"sensitivity_{parameter_name}_{region_name}.png"
        fig.savefig(out_png, dpi=150)
        plt.close(fig)

        print(f"Sensitivity plot → {out_png}")

    # Print a compact terminal summary.
    print("\nSensitivity summary:")
    cols_to_print = [
        "parameter",
        "value",
        "recovering_well_pct_burned",
        "recovering_weak_pct_burned",
        "failed_pct_burned",
        "class4_outside_assessment_pct_total",
    ]
    print(df[cols_to_print].to_string(index=False))


if SENS_ENABLED:
    run_threshold_sensitivity_analysis()
else:
    print("Sensitivity analysis disabled.")

# ── Failed-pixel diagnostics ─────────────────────────────────────────────────
# We pre-compute monthly time series for a few pixels classified as:
#   3 = Not recovering / failed
#
# The output is:
#   1. A CSV file with the extracted time series
#   2. Clickable markers on the Folium map with a small HTML table

diag_cfg = cfg.get("failed_diagnostics", {}) or {}
FAILED_DIAG_ENABLED = bool(diag_cfg.get("enabled", True))
FAILED_DIAG_N = int(diag_cfg.get("n_examples", 5))
FAILED_DIAG_RADIUS_M = int(diag_cfg.get("radius_m", 40))

failed_examples = []
failed_diag_rows = []


def _fmt(v, digits=3):
    """Format numbers nicely for HTML tables."""
    if v is None:
        return ""
    try:
        return f"{float(v):.{digits}f}"
    except Exception:
        return str(v)


def sample_failed_locations() -> list[dict]:
    """
    Select diagnostic locations from the largest connected clusters of
    recovery class 3 = Not recovering / failed.

    This is better than random sampling because it guarantees that the main
    visible red patches are represented.
    """
    scale = int(exp.get("scale", 20))

    failed_mask = recovery_class.eq(3).selfMask().rename("failed")

    failed_vectors = failed_mask.reduceToVectors(
        geometry=aoi,
        scale=scale,
        geometryType="polygon",
        eightConnected=True,
        labelProperty="failed",
        reducer=ee.Reducer.countEvery(),
        maxPixels=int(exp.get("max_pixels", 1e9)),
        tileScale=4,
    )

    def _add_cluster_props(f):
        centroid = f.geometry().centroid(1)
        coords = centroid.coordinates()
        area_m2 = f.geometry().area(1)
        return f.set({
            "lon": coords.get(0),
            "lat": coords.get(1),
            "area_m2": area_m2,
            "pixel_count_est": area_m2.divide(scale * scale),
        })

    top_clusters = (
        failed_vectors
        .map(_add_cluster_props)
        .sort("area_m2", False)
        .limit(FAILED_DIAG_N)
    )

    features = top_clusters.getInfo().get("features", [])

    out = []
    for i, f in enumerate(features):
        props = f["properties"]

        lon = float(props["lon"])
        lat = float(props["lat"])
        area_m2 = float(props["area_m2"])
        pixel_count_est = float(props["pixel_count_est"])

        out.append({
            "id": f"failed_cluster_{i+1}",
            "cluster_rank": i + 1,
            "lon": lon,
            "lat": lat,
            "area_m2": area_m2,
            "pixel_count_est": pixel_count_est,
            "geometry": ee.Geometry.Point([lon, lat]).buffer(FAILED_DIAG_RADIUS_M),
        })

    return out

def extract_failed_timeseries(example: dict) -> list[dict]:
    """
    Extract monthly mean spectral indices around one failed location.
    RRI is recomputed for every monthly composite.
    """
    rows = []
    geom = example["geometry"]

    for ws, we, label in windows:
        try:
            current = make_composite(ws, we)

            rri_t = compute_rri(
                recent=current,
                postfire=postfire_composite,
                prefire=prefire_composite,
            )

            img = (
                current
                .select(["NDVI", "NDRE", "NDMI", "BSI", "NDWI", "NBR"])
                .addBands(rri_t)
            )

            stats = img.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=geom,
                scale=int(exp.get("scale", 20)),
                maxPixels=int(exp.get("max_pixels", 1e9)),
                bestEffort=True,
            ).getInfo() or {}

        except ee.EEException as e:
            print(f"  {example['id']} {label} skipped: {e}")
            stats = {}

        rows.append({
            "example_id": example["id"],
            "cluster_rank": example.get("cluster_rank"),
            "lon": example["lon"],
            "lat": example["lat"],
            "radius_m": FAILED_DIAG_RADIUS_M,
            "cluster_area_m2": example.get("area_m2"),
            "cluster_pixel_count_est": example.get("pixel_count_est"),
            "recovery_class": 3,
            "date": ws,
            "label": label,
            "NDVI": stats.get("NDVI"),
            "RRI": stats.get("RRI"),
            "NDRE": stats.get("NDRE"),
            "BSI": stats.get("BSI"),
            "NDMI": stats.get("NDMI"),
            "NDWI": stats.get("NDWI"),
            "NBR": stats.get("NBR"),
})

    return rows


def make_failed_popup_table(example_id: str) -> str:
    """
    Build a compact scrollable HTML table for the Folium popup.
    """
    rows = [r for r in failed_diag_rows if r["example_id"] == example_id]

    table_rows = ""
    for r in rows:
        table_rows += (
            "<tr>"
            f"<td>{r['label']}</td>"
            f"<td>{_fmt(r['NDVI'])}</td>"
            f"<td>{_fmt(r['RRI'])}</td>"
            f"<td>{_fmt(r['NDRE'])}</td>"
            f"<td>{_fmt(r['BSI'])}</td>"
            f"<td>{_fmt(r['NDMI'])}</td>"
            f"<td>{_fmt(r['NBR'])}</td>"
            "</tr>"
        )

    return f"""
    <div style="max-height:300px; overflow-y:auto;">
      <table style="border-collapse:collapse; font-size:11px; width:100%;">
        <thead>
          <tr>
            <th style="border-bottom:1px solid #999;">Date</th>
            <th style="border-bottom:1px solid #999;">NDVI</th>
            <th style="border-bottom:1px solid #999;">RRI</th>
            <th style="border-bottom:1px solid #999;">NDRE</th>
            <th style="border-bottom:1px solid #999;">BSI</th>
            <th style="border-bottom:1px solid #999;">NDMI</th>
            <th style="border-bottom:1px solid #999;">NBR</th>
          </tr>
        </thead>
        <tbody>
          {table_rows}
        </tbody>
      </table>
    </div>
    """


if FAILED_DIAG_ENABLED:
    import csv

    print("Sampling failed recovery pixels for diagnostics...")
    failed_examples = sample_failed_locations()

    print(f"  sampled {len(failed_examples)} failed diagnostic locations")

    for ex in failed_examples:
        print(f"Extracting time series for {ex['id']}...")
        failed_diag_rows.extend(extract_failed_timeseries(ex))

    out_diag_csv = SRC_DIR / f"failed_diagnostics_{region_name}.csv"

    if failed_diag_rows:
        with open(out_diag_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "example_id", "cluster_rank",
                    "lon", "lat", "radius_m",
                    "cluster_area_m2", "cluster_pixel_count_est",
                    "recovery_class",
                    "date", "label",
                    "NDVI", "RRI", "NDRE", "BSI", "NDMI", "NDWI", "NBR",
]
            )
            writer.writeheader()
            writer.writerows(failed_diag_rows)

        print(f"Failed-pixel diagnostics CSV → {out_diag_csv}")
    else:
        print("No failed pixels found for diagnostics.")

RECOVERY_CLASSES = [
    (0, "Uncertain / water / insufficient data", "#888888"),
    (1, "Recovering well",                       "#1f5e1f"),
    (2, "Recovering, but weak",                 "#a8d666"),
    (3, "Not recovering / failed",              "#ff1a1a"),  # glowing/bright red
    (4, "Outside burn recovery assessment",             "#4b2e1f"),  # dark brown
]
RECOVERY_PALETTE = [c for _, _, c in sorted(RECOVERY_CLASSES)]
VIS_RECOVERY_CLASS = {"min": 0, "max": 4, "palette": RECOVERY_PALETTE}

GUARD_CLASSES = [
    (0, "Clear recovery signal",                  "#2f7f3f"),
    (1, "Water / insufficient observations",      "#888888"),
    (2, "Possible later disturbance",             "#7d2222"),
    (3, "Mixed recovery signal",                  "#f4a261"),
    (4, "Outside burn recovery assessment",       "#4a7c59"),
]
GUARD_PALETTE = [c for _, _, c in sorted(GUARD_CLASSES)]
VIS_GUARDS = {"min": 0, "max": 4, "palette": GUARD_PALETTE}

# Diagnostic: per-class pixel count and percentage over the AOI.
try:
    _hist = recovery_class.reduceRegion(
        reducer   = ee.Reducer.frequencyHistogram(),
        geometry  = aoi,
        scale     = int(exp.get("scale", 20)),
        maxPixels = int(exp.get("max_pixels", 1e9)),
        bestEffort= True,
    ).getInfo()
    _counts = _hist.get("recovery_class") or {}
    # frequencyHistogram returns string-keyed counts; coerce + total
    _counts = {int(float(k)): int(v) for k, v in _counts.items()}
    _total  = sum(_counts.values()) or 1
    print(f"Recovery class histogram ({region_name}):")
    for _idx, _name, _ in RECOVERY_CLASSES:
        _n = _counts.get(_idx, 0)
        _pct = 100.0 * _n / _total
        print(f"  {_name:26s}  {_n:>10,d} px  ({_pct:5.1f}%)")
except ee.EEException as _e:
    print(f"  (recovery-class histogram unavailable: {_e})")

try:
    _hist = guard_flags.reduceRegion(
        reducer   = ee.Reducer.frequencyHistogram(),
        geometry  = aoi,
        scale     = int(exp.get("scale", 20)),
        maxPixels = int(exp.get("max_pixels", 1e9)),
        bestEffort= True,
    ).getInfo()
    _counts = _hist.get("guard_flags") or {}
    _counts = {int(float(k)): int(v) for k, v in _counts.items()}
    _total  = sum(_counts.values()) or 1
    print(f"Guard-flag histogram ({region_name}):")
    for _idx, _name, _ in GUARD_CLASSES:
        _n = _counts.get(_idx, 0)
        _pct = 100.0 * _n / _total
        print(f"  {_name:34s}  {_n:>10,d} px  ({_pct:5.1f}%)")
except ee.EEException as _e:
    print(f"  (guard-flag histogram unavailable: {_e})")

VIS = {
    "RGB":    {"bands": ["B4","B3","B2"], "params": {"min":0,"max":2800,"gamma":1.4}},
    "NDVI":   {"bands": ["NDVI"],  "params": {"min":-0.1,"max":0.85,
               "palette":["saddlebrown","khaki","limegreen","darkgreen"]}},
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

# ── Failed diagnostic graphs from CSV ────────────────────────────────────────
# This reads failed_diagnostics_<region>.csv, keeps only failed examples,
# limits them to FAILED_DIAG_N, creates one PNG graph per example, and embeds
# those graphs in the Folium map popups.

failed_diag_plot_paths = {}
failed_diag_points = []


def build_failed_diagnostic_plots_from_csv() -> None:
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    csv_path = SRC_DIR / f"failed_diagnostics_{region_name}.csv"

    if not csv_path.exists():
        print(f"No failed diagnostics CSV found: {csv_path}")
        return

    df = pd.read_csv(csv_path)

    # Keep only failed pixels if the CSV contains the recovery class.
    if "recovery_class" in df.columns:
        df = df[df["recovery_class"].astype(int) == 3]

    # Exclude outside burn recovery assessment if this column exists.
    if "guard_flag" in df.columns:
        df = df[df["guard_flag"].astype(int) != 4]

    # Force only the requested number of examples.
    # If cluster_rank exists, keep the largest clusters first.
    if "cluster_rank" in df.columns:
        keep_ids = (
            df[["example_id", "cluster_rank"]]
            .drop_duplicates()
            .sort_values("cluster_rank")
            .head(FAILED_DIAG_N)["example_id"]
            .tolist()
        )
    else:
        keep_ids = (
            df["example_id"]
            .drop_duplicates()
            .head(FAILED_DIAG_N)
            .tolist()
        )

    df = df[df["example_id"].isin(keep_ids)]

    if df.empty:
        print("No failed diagnostic points left after filtering.")
        return

    df["date"] = pd.to_datetime(df["date"])

    for ex_id, sub in df.groupby("example_id"):
        sub = sub.sort_values("date").copy()

        lon = float(sub["lon"].iloc[0])
        lat = float(sub["lat"].iloc[0])

        if "radius_m" in sub.columns:
            radius_m = int(sub["radius_m"].iloc[0])
        else:
            radius_m = FAILED_DIAG_RADIUS_M

        fig, axes = plt.subplots(3, 1, figsize=(8.5, 7), sharex=True)

        # Panel 1 — recovery / canopy signal
        axes[0].plot(
            sub["date"], sub["NDVI"],
            marker="o", linewidth=1.5, label="NDVI"
        )
        axes[0].plot(
            sub["date"], sub["RRI"],
            marker="s", linestyle="--", linewidth=1.5, label="RRI"
        )
        axes[0].plot(
            sub["date"], sub["NDRE"],
            marker="o", linewidth=1.2, label="NDRE"
        )
        axes[0].axhline(
            T_RRI_LOW,
            linestyle=":",
            linewidth=1,
            label=f"RRI failed threshold ({T_RRI_LOW})",
        )
        axes[0].set_ylabel("Recovery / canopy")
        axes[0].legend(fontsize=8, ncol=2)
        axes[0].grid(alpha=0.25)

        # Panel 2 — soil / moisture signal
        axes[1].plot(
            sub["date"], sub["BSI"],
            marker="o", linewidth=1.5, label="BSI"
        )
        axes[1].plot(
            sub["date"], sub["NDMI"],
            marker="o", linewidth=1.5, label="NDMI"
        )
        axes[1].plot(
            sub["date"], sub["NDWI"],
            marker="o", linewidth=1.2, label="NDWI"
        )
        axes[1].axhline(
            T_BSI_HIGH,
            linestyle=":",
            linewidth=1,
            label=f"High BSI threshold ({T_BSI_HIGH})",
        )
        axes[1].set_ylabel("Soil / moisture")
        axes[1].legend(fontsize=8, ncol=2)
        axes[1].grid(alpha=0.25)

        # Panel 3 — burn / disturbance signal
        axes[2].plot(
            sub["date"], sub["NBR"],
            marker="o", linewidth=1.5, label="NBR"
        )
        axes[2].set_ylabel("NBR")
        axes[2].legend(fontsize=8)
        axes[2].grid(alpha=0.25)

        if PLANTING is not None:
            for ax in axes:
                ax.axvline(PLANTING, linestyle="--", linewidth=1)

        axes[2].xaxis.set_major_locator(mdates.YearLocator())
        axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

        fig.suptitle(
            f"{ex_id} — failed recovery diagnostic\n"
            f"lon={lon:.5f}, lat={lat:.5f}, radius={radius_m} m",
            fontsize=11,
        )

        plt.tight_layout()

        out_png = SRC_DIR / f"{ex_id}_{region_name}_diagnostic.png"
        fig.savefig(out_png, dpi=140)
        plt.close(fig)

        failed_diag_plot_paths[ex_id] = out_png
        failed_diag_points.append({
            "id": ex_id,
            "lon": lon,
            "lat": lat,
            "radius_m": radius_m,
        })

        print(f"Saved diagnostic plot → {out_png}")


if FAILED_DIAG_ENABLED:
    build_failed_diagnostic_plots_from_csv()


# ── Folium map export ────────────────────────────────────────────────────────
try:
    import folium
    import base64

    # Centre the map on the AOI bounds.
    _lons = [c[0] for c in bounds]
    _lats = [c[1] for c in bounds]
    _center = [(min(_lats) + max(_lats)) / 2, (min(_lons) + max(_lons)) / 2]

    fmap = folium.Map(
        location=_center,
        zoom_start=10,
        tiles="OpenStreetMap",
        control_scale=True,
    )

    def _add_ee_layer(
        ee_image: ee.Image,
        params: dict,
        name: str,
        show: bool = True,
    ) -> None:
        map_id = ee_image.getMapId(params)

        folium.TileLayer(
            tiles=map_id["tile_fetcher"].url_format,
            attr="Google Earth Engine",
            name=name,
            overlay=True,
            control=True,
            show=show,
        ).add_to(fmap)

    # Single-period summer raster layers.
    for title, vc in VIS_AVAILABLE.items():
        _add_ee_layer(
            raster.select(vc["bands"]),
            vc["params"],
            title,
            show=True,
        )

    # RRI continuous layer.
    _add_ee_layer(
        rri_img,
        {
            "min": 0,
            "max": 1,
            "palette": [
                "#7d2222",
                "#ff6600",
                "#ffe066",
                "#a8d666",
                "#1f5e1f",
            ],
        },
        f"RRI (good≥{T_RRI_GOOD}, failed<{T_RRI_LOW})",
        show=False,
    )

    # Main recovery class map.
    _add_ee_layer(
        recovery_class,
        VIS_RECOVERY_CLASS,
        "Recovery class",
        show=True,
    )

    # Uncertainty / guard flags.
    _add_ee_layer(
        guard_flags,
        VIS_GUARDS,
        "Uncertainty flags",
        show=False,
    )

    # Forest mask.
    if forest_mask_img is not None:
        _add_ee_layer(
            forest_mask_img,
            {
                "min": 0,
                "max": 1,
                "palette": ["#dddddd", "#1a6b1a"],
            },
            f"Forest mask (DW {FM_YEAR}, p≥{FM_THRESHOLD})",
            show=False,
        )

    # Pre-fire, post-fire, and recent RGB composites.
    _RGB_PARAMS = {
        "min": 0,
        "max": 2800,
        "gamma": 1.4,
    }

    _add_ee_layer(
        prefire_composite.select(["B4", "B3", "B2"]),
        _RGB_PARAMS,
        f"Pre-fire RGB ({PREFIRE_YEARS[0]} Apr–May)",
        show=False,
    )

    _add_ee_layer(
        postfire_composite.select(["B4", "B3", "B2"]),
        _RGB_PARAMS,
        f"Post-fire RGB ({POSTFIRE_YEARS[0]}–{POSTFIRE_YEARS[1]})",
        show=False,
    )

    _add_ee_layer(
        recent_composite.select(["B4", "B3", "B2"]),
        _RGB_PARAMS,
        f"Recent RGB ({RECENT_YEARS[0]}–{RECENT_YEARS[1]})",
        show=False,
    )

    # AOI outline.
    folium.GeoJson(
        aoi.getInfo(),
        name="AOI",
        style_function=lambda _f: {
            "color": "red",
            "weight": 2,
            "fillOpacity": 0,
        },
    ).add_to(fmap)

    # Failed diagnostic graph markers.
    if FAILED_DIAG_ENABLED and failed_diag_points:
        failed_group = folium.FeatureGroup(
            name="Failed-pixel diagnostic graphs",
            show=True,
        )

        for ex in failed_diag_points:
            ex_id = ex["id"]
            plot_path = failed_diag_plot_paths.get(ex_id)

            if plot_path is None or not Path(plot_path).exists():
                continue

            with open(plot_path, "rb") as f:
                img64 = base64.b64encode(f.read()).decode("utf-8")

            popup_html = f"""
            <div style="width:720px;">
              <h4 style="margin:4px 0 6px 0;">{ex_id}</h4>
              <div style="font-size:12px; margin-bottom:8px;">
                <b>Recovery class:</b> 3 — Not recovering / failed<br>
                <b>Location:</b> {ex['lat']:.5f}, {ex['lon']:.5f}<br>
                <b>Buffer radius:</b> {ex['radius_m']} m
              </div>
              <img src="data:image/png;base64,{img64}" width="700">
            </div>
            """

            folium.CircleMarker(
                location=[ex["lat"], ex["lon"]],
                radius=7,
                color="#ff0000",
                fill=True,
                fill_color="#ff0000",
                fill_opacity=0.9,
                popup=folium.Popup(popup_html, max_width=750),
                tooltip=f"{ex_id} — failed recovery graphs",
            ).add_to(failed_group)

        failed_group.add_to(fmap)

    # Legend.
    _legend_rows = "".join(
        f'<div style="display:flex;align-items:center;margin:2px 0;">'
        f'  <span style="display:inline-block;width:14px;height:14px;'
        f'background:{c};border:1px solid #555;margin-right:6px;"></span>'
        f'  <span>{idx} — {name}</span>'
        f'</div>'
        for idx, name, c in RECOVERY_CLASSES
    )

    _stress_rows = "".join(
        f'<div style="display:flex;align-items:center;margin:2px 0;">'
        f'  <span style="display:inline-block;width:14px;height:14px;'
        f'background:{c};border:1px solid #555;margin-right:6px;"></span>'
        f'  <span>{idx} — {name}</span>'
        f'</div>'
        for idx, name, c in GUARD_CLASSES
    )

    _legend_html = (
        '<div style="position: fixed; bottom: 30px; right: 30px; z-index: 9999;'
        '  background: rgba(255,255,255,0.95); padding: 10px 12px;'
        '  border: 1px solid #888; border-radius: 4px;'
        '  font: 12px/1.3 system-ui, sans-serif; max-width: 360px;'
        '  box-shadow: 0 2px 6px rgba(0,0,0,0.2);">'
        '<div style="font-weight:600;margin-bottom:6px;">Recovery class</div>'
        f'{_legend_rows}'
        '<div style="font-weight:600;margin:10px 0 6px 0;">Uncertainty flags</div>'
        f'{_stress_rows}'
        '<div style="font-size:10px;color:#666;margin-top:6px;">'
        'Recovery class expresses a spectral trajectory consistent with woody '
        'vegetation establishment, not individual tree health.'
        '</div>'
        '</div>'
    )

    fmap.get_root().html.add_child(folium.Element(_legend_html))

    folium.LayerControl(collapsed=False).add_to(fmap)

    out_map = SRC_DIR / f"map_{region_name}.html"
    fmap.save(str(out_map))

    print(f"Interactive map → {out_map}")


except ImportError:
    import urllib.request
    import io
    import matplotlib.pyplot as plt
    from PIL import Image

    n_panels = len(VIS_AVAILABLE)

    fig3, axes3 = plt.subplots(1, n_panels, figsize=(4 * n_panels, 4))

    if n_panels == 1:
        axes3 = [axes3]

    fig3.suptitle(
        f"{region_name} — {RASTER_START} to {RASTER_END}",
        fontsize=11,
    )

    for ax, (title, vc) in zip(axes3, VIS_AVAILABLE.items()):
        url = raster.select(vc["bands"]).getThumbURL(
            {
                **vc["params"],
                "region": bounds,
                "dimensions": 512,
                "format": "png",
            }
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

    # Change-band thumbnails.
    fig4, axes4 = plt.subplots(2, 3, figsize=(12, 8))

    fig4.suptitle(
        f"{region_name} — Change layers "
        f"(recent {RECENT_YEARS[0]}–{RECENT_YEARS[1]} − "
        f"post-fire {POSTFIRE_YEARS[0]}–{POSTFIRE_YEARS[1]}, months "
        f"{CMP_MONTHS[0]}–{CMP_MONTHS[1] - 1})",
        fontsize=11,
    )

    for ax, (band, params) in zip(axes4.flat, VIS_CHANGE.items()):
        url = change_img.select(band).getThumbURL(
            {
                **params,
                "region": bounds,
                "dimensions": 512,
                "format": "png",
            }
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