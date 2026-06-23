"""
EO4Change Group 4 — configuration, constants, GEE init, AOI, date windows.
All pipeline parameters driven by a YAML config (default: data/config_DK.yaml).
Import this module first; it executes GEE initialisation as a side effect.
"""
import argparse
import ee
import yaml
import json
from pathlib import Path
from datetime import datetime
from dateutil.relativedelta import relativedelta

# ── Paths & CLI ───────────────────────────────────────────────────────────────
SRC_DIR     = Path(__file__).resolve().parent
DATA_DIR    = SRC_DIR.parent / "data"
FIGURES_DIR = SRC_DIR / "figures"
HTML_DIR    = SRC_DIR / "html"
CSV_DIR     = SRC_DIR / "csv"

for _d in (FIGURES_DIR, HTML_DIR, CSV_DIR):
    _d.mkdir(exist_ok=True)

_parser = argparse.ArgumentParser(add_help=True,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
_parser.add_argument("--config", default="config_DK.yaml",
                     help="config filename inside data/ (or absolute path). Default: config_DK.yaml")
_parser.add_argument("--map-only", action="store_true",
                     help="Skip the slow per-window time-series extraction and plots; only build composites + change image + map.")
_parser.add_argument("--force-ts", action="store_true",
                     help="Force re-extraction of time-series from GEE even if a cached CSV already exists.")
_args, _ = _parser.parse_known_args()
MAP_ONLY = _args.map_only
FORCE_TS = _args.force_ts

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

SPECTRAL_BANDS = ["NDVI", "NDRE", "NDMI", "BSI", "NDWI", "NBR"]
RGB_BANDS = ["B4", "B3", "B2"]
ALL_BANDS = SPECTRAL_BANDS + BIOPHYS_VARS
COMPOSITE_BANDS = RGB_BANDS + ALL_BANDS

# Comparison periods
cmp_cfg = cfg.get("comparison", {}) or {}
PREFIRE_YEARS  = tuple(cmp_cfg.get("prefire_years",  [2017, 2017]))
PREFIRE_MONTHS = tuple(cmp_cfg.get("prefire_months", [4, 6]))
POSTFIRE_YEARS = tuple(cmp_cfg.get("postfire_years",
                        cmp_cfg.get("early_years",   [2017, 2017])))
RECENT_YEARS   = tuple(cmp_cfg.get("recent_years",  [2022, 2025]))
CMP_MONTHS     = tuple(cmp_cfg.get("months",        [5, 9]))

# Forest mask
fm_cfg = cfg.get("forest_mask", {}) or {}
FM_ENABLED   = bool(fm_cfg.get("enabled", False))
FM_YEAR      = int(fm_cfg.get("year", 2016))
FM_MONTHS    = tuple(fm_cfg.get("months", [5, 9]))
FM_THRESHOLD = float(fm_cfg.get("tree_prob_threshold", 0.40))

# Recovery class and guard-flag thresholds
cls_cfg = cfg.get("classification", {}) or {}
T_RRI_GOOD      = float(cls_cfg.get("rri_good_threshold",        0.80))
T_RRI_LOW       = float(cls_cfg.get("rri_low_threshold",         0.30))
T_RRI_MIN_DENOM = float(cls_cfg.get("rri_min_denominator",       0.10))
T_BURN          = float(cls_cfg.get("burn_severity_threshold",   0.10))
T_NDRE_GOOD  = float(cls_cfg.get("ndre_good_threshold",          0.18))
T_DNDVI_GOOD = float(cls_cfg.get("ndvi_change_good_threshold",   0.09))
T_BSI_HIGH   = float(cls_cfg.get("bsi_high_threshold",           0.09))
T_WATER      = float(cls_cfg.get("ndwi_water_threshold",         0.15))
T_DISTURB    = float(cls_cfg.get("disturbance_threshold",        0.09))

# Sensitivity analysis
sens_cfg = cfg.get("sensitivity_analysis", {}) or {}
SENS_ENABLED = bool(sens_cfg.get("enabled", False))
SENS_BURN_VALUES = sens_cfg.get("burn_severity_threshold_values", [0.05, 0.075, T_BURN, 0.125, 0.15])
SENS_RRI_GOOD_VALUES = sens_cfg.get("rri_good_threshold_values", [0.70, 0.75, T_RRI_GOOD, 0.85, 0.90])
SENS_RRI_LOW_VALUES = sens_cfg.get("rri_low_threshold_values", [0.20, 0.25, T_RRI_LOW, 0.35, 0.40])

# Failed diagnostics
diag_cfg = cfg.get("failed_diagnostics", {}) or {}
FAILED_DIAG_ENABLED = bool(diag_cfg.get("enabled", True))
FAILED_DIAG_N = int(diag_cfg.get("n_examples", 5))
FAILED_DIAG_RADIUS_M = int(diag_cfg.get("radius_m", 40))

# Reference event
REF_DATE  = tmp.get("reference_date",  "2022-04-01")
REF_LABEL = tmp.get("reference_label", "Reference T₀")
PLANTING_DATE = datetime.strptime(REF_DATE, "%Y-%m-%d") if REF_DATE else None

region_name = sp.get("region_name", "AOI")

# Cloud masking
CS_BAND   = opt.get("csplus_band", "cs")
CS_THRESH = opt.get("cs_plus_threshold", 0.65)
MAX_CLOUD = opt.get("max_cloud_cover", 50)

# Plot colors
SPEC_COLORS  = {
    "NDVI": "#2d8a2d", "NDRE": "#7baf27", "NDMI": "#1e6eb5", "BSI": "#c47a1e",
    "NDWI": "#1ab1c4", "NBR":  "#a83232",
}
BIOPH_COLORS = {"laie": "#8b4513", "fapar": "#ff6600", "fcover": "#009900"}
BIOPH_LABELS = {"laie": "LAI-e [m²/m²]", "fapar": "FAPAR [0-1]", "fcover": "FCOVER [0-1]"}

# ── GEE init ──────────────────────────────────────────────────────────────────
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
