"""
EO4Change Group 4 — Sentinel-2 reforestation / vegetation monitor.
All parameters driven by a YAML config (default: data/config_DK.yaml).

Run: python src/gee.py [--config config_pt.yaml] [--map-only]
"""
import csv as _csv
import ee
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from gee_config import (
    cfg, sp, tmp, exp, MAP_ONLY, FORCE_TS, region_name,
    aoi, windows,
    PREFIRE_YEARS, PREFIRE_MONTHS, POSTFIRE_YEARS, RECENT_YEARS, CMP_MONTHS,
    FM_ENABLED, FM_THRESHOLD, SENS_ENABLED, FAILED_DIAG_ENABLED, FAILED_DIAG_N,
    T_RRI_GOOD, T_RRI_LOW, REF_DATE, REF_LABEL, SRC_DIR,
    FIGURES_DIR, HTML_DIR, CSV_DIR,
    BIOPHYS_VARS, ALL_BANDS, SPECTRAL_BANDS, PLANTING_DATE,
    SPEC_COLORS, BIOPH_COLORS, BIOPH_LABELS,
)
from gee_s2 import make_composite, make_seasonal_composite
from gee_recovery import (
    build_forest_mask, compute_rri,
    build_recovery_class, build_guard_flags,
    RECOVERY_CLASSES, GUARD_CLASSES,
)
from gee_sensitivity import run_threshold_sensitivity_analysis
from gee_diagnostics import (
    sample_failed_locations, extract_failed_timeseries,
    build_failed_diagnostic_plots_from_csv,
)
from gee_map import build_folium_map

# ── Output hygiene: migrate old CSVs, wipe stale figures ─────────────────────
import shutil as _shutil

# Move any CSVs sitting at the old src/ root into csv/
for _legacy_csv in [
    SRC_DIR / f"timeseries_{region_name}.csv",
    SRC_DIR / f"failed_diagnostics_{region_name}.csv",
    SRC_DIR / f"sensitivity_thresholds_{region_name}.csv",
    SRC_DIR / f"sensitivity_{region_name}.csv",
]:
    _dest = CSV_DIR / _legacy_csv.name
    if _legacy_csv.exists() and not _dest.exists():
        _shutil.move(str(_legacy_csv), str(_dest))
        print(f"Migrated → csv/{_legacy_csv.name}")

# Delete all region-specific PNGs (figures/ and legacy src/ root) so only
# fresh images remain after this run.
for _stale in (
    list(FIGURES_DIR.glob(f"*{region_name}*.png"))
    + list(SRC_DIR.glob(f"*{region_name}*.png"))
):
    _stale.unlink()
    print(f"Removed stale figure: {_stale.name}")

# ── Shared plot style ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "#f8f8f8",
    "axes.grid": True,
    "grid.alpha": 0.4,
    "grid.linestyle": "--",
    "grid.color": "#cccccc",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.family": "DejaVu Sans",
    "font.size": 10,
})
_RRI_COLOR = "#9b59b6"

# ── Build prefire/postfire composites (needed before TS loop for RRI) ─────────
print(f"Building pre-fire composite (months {PREFIRE_MONTHS[0]}–{PREFIRE_MONTHS[1]-1}, "
      f"year {PREFIRE_YEARS[0]})…")
prefire_composite  = make_seasonal_composite(PREFIRE_YEARS[0],  PREFIRE_YEARS[1],
                                             PREFIRE_MONTHS[0],  PREFIRE_MONTHS[1])
print(f"Building post-fire composite (months {CMP_MONTHS[0]}–{CMP_MONTHS[1]-1}, "
      f"years {POSTFIRE_YEARS[0]}–{POSTFIRE_YEARS[1]})…")
postfire_composite = make_seasonal_composite(POSTFIRE_YEARS[0], POSTFIRE_YEARS[1],
                                             CMP_MONTHS[0],      CMP_MONTHS[1])

# ── Time-series extraction (or load from cache) ───────────────────────────────
_ts_csv = CSV_DIR / f"timeseries_{region_name}.csv"
_ts_fields = ["date", "label"] + ALL_BANDS + ["RRI"]

_extracted_fresh = False
if _ts_csv.exists() and not FORCE_TS:
    print(f"Loading cached time-series from {_ts_csv.name}  (--force-ts to re-extract)…")
    with open(_ts_csv, encoding="utf-8") as _f:
        _reader = _csv.DictReader(_f)
        records = []
        for _row in _reader:
            for _k in ALL_BANDS + ["RRI"]:
                _raw = _row.get(_k, "")
                _row[_k] = float(_raw) if _raw not in ("", "None", None) else None
            records.append(dict(_row))
    print(f"  Loaded {len(records)} records.")
    if records and all(r.get("RRI") is None for r in records):
        print("  WARNING: RRI column is empty in the cached CSV (CSV predates RRI extraction).")
        print("  Run with --force-ts to re-extract and populate RRI.")
elif MAP_ONLY:
    print("Extracting time series — skipped (--map-only, no cache found).")
    records: list[dict] = []
else:
    print("Extracting time series — please wait…")
    records = []
    for ws, we, label in windows:
        try:
            current = make_composite(ws, we)
            rri_t   = compute_rri(current, postfire_composite, prefire_composite)
            img     = current.select(ALL_BANDS).addBands(rri_t)
            row = img.reduceRegion(
                reducer   = ee.Reducer.mean(),
                geometry  = aoi,
                scale     = int(exp.get("scale", 20)),
                maxPixels = int(exp.get("max_pixels", 1e9)),
            ).getInfo() or {}
        except ee.EEException as e:
            print(f"  {label:12s}  skipped ({e})")
            row = {}
        row["date"]  = ws
        row["label"] = label
        records.append(row)
        v = row.get("NDVI")
        print(f"  {label:12s}  NDVI={v:.3f}" if v is not None else f"  {label:12s}  (no data)")
    _extracted_fresh = True

# ── CSV export (only when freshly extracted) ──────────────────────────────────
if _extracted_fresh and records:
    with open(_ts_csv, "w", newline="", encoding="utf-8") as _f:
        _w = _csv.DictWriter(_f, fieldnames=_ts_fields, extrasaction="ignore")
        _w.writeheader()
        _w.writerows(records)
    print(f"Saved → {_ts_csv}")

# ── Individual time-series plots ──────────────────────────────────────────────
if records:
    PLANTING = PLANTING_DATE
    dates = [datetime.strptime(r["date"], "%Y-%m-%d") for r in records]

    def _save_ts_plot(band: str, color: str, title: str, ylabel: str) -> None:
        vals = [r.get(band) for r in records]
        x = [d for d, v in zip(dates, vals) if v is not None]
        y = [v for v in vals if v is not None]
        if not x:
            print(f"  (no data for {band}, skipping plot)")
            return
        fig, ax = plt.subplots(figsize=(13, 4))
        ax.plot(x, y, marker="o", ms=4, lw=1.8, color=color, zorder=3,
                markerfacecolor="white", markeredgecolor=color, markeredgewidth=1.3)
        if PLANTING is not None:
            ax.axvline(PLANTING, color="#cc0000", ls="--", lw=1.5,
                       label=REF_LABEL or "Reference event", zorder=4)
            ax.legend(fontsize=9, framealpha=0.8, loc="upper left")
        ax.set_title(f"{region_name}  —  {title}", fontsize=13, fontweight="bold", pad=10)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[4, 7, 10]))
        ax.tick_params(axis="both", labelsize=9)
        fig.autofmt_xdate()
        plt.tight_layout()
        out = FIGURES_DIR / f"timeseries_{band.lower()}_{region_name}.png"
        fig.savefig(out, dpi=180, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved → {out}")

    _save_ts_plot("RRI",  _RRI_COLOR,          "Relative Recovery Indicator (RRI)",
                  "RRI  (0 = post-fire level,  1 = pre-fire level)")
    _save_ts_plot("NDVI", SPEC_COLORS["NDVI"], "NDVI — Vegetation greenness",          "NDVI  (−1 … 1)")
    _save_ts_plot("NDRE", SPEC_COLORS["NDRE"], "NDRE — Chlorophyll / canopy quality",  "NDRE  (−1 … 1)")
    _save_ts_plot("NBR",  SPEC_COLORS["NBR"],  "NBR — Burn severity / fire recovery",  "NBR  (−1 … 1)")

# ── Build raster + recent composite ──────────────────────────────────────────
RASTER_START, RASTER_END = "2023-06-01", "2023-08-31"
raster = make_composite(RASTER_START, RASTER_END)
bounds = aoi.bounds().getInfo()["coordinates"][0]

print(f"Building recent composite (months {CMP_MONTHS[0]}–{CMP_MONTHS[1]-1}, "
      f"years {RECENT_YEARS[0]}–{RECENT_YEARS[1]})…")
recent_composite = make_seasonal_composite(RECENT_YEARS[0], RECENT_YEARS[1],
                                           CMP_MONTHS[0],    CMP_MONTHS[1])

# ── Change image ──────────────────────────────────────────────────────────────
print("Computing change image (recent − post-fire)…")
CHANGE_BANDS = [f"{b}_change" for b in SPECTRAL_BANDS]
change_img = (recent_composite.select(SPECTRAL_BANDS)
                              .subtract(postfire_composite.select(SPECTRAL_BANDS))
                              .rename(CHANGE_BANDS))

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

# ── Forest mask ───────────────────────────────────────────────────────────────
if FM_ENABLED:
    print(f"Building forest mask (Dynamic World, tree prob >= {FM_THRESHOLD})…")
forest_mask_img = build_forest_mask()

# ── Recovery products ─────────────────────────────────────────────────────────
print("Computing RRI…")
rri_img = compute_rri(recent_composite, postfire_composite, prefire_composite)

print("Building recovery class…")
recovery_class = build_recovery_class(recent_composite, postfire_composite,
                                      prefire_composite, change_img)

print("Building guard flags…")
guard_flags = build_guard_flags(recent_composite, postfire_composite,
                                prefire_composite, change_img)

if forest_mask_img is not None:
    recovery_class = recovery_class.updateMask(forest_mask_img)
    guard_flags    = guard_flags.updateMask(forest_mask_img)

# ── Recovery class histogram ──────────────────────────────────────────────────
_hist_counts: dict = {}
try:
    _hist = recovery_class.reduceRegion(
        reducer   = ee.Reducer.frequencyHistogram(),
        geometry  = aoi,
        scale     = int(exp.get("scale", 20)),
        maxPixels = int(exp.get("max_pixels", 1e9)),
        bestEffort= True,
    ).getInfo()
    _counts = _hist.get("recovery_class") or {}
    _counts = {int(float(k)): int(v) for k, v in _counts.items()}
    _hist_counts = _counts
    _total  = sum(_counts.values()) or 1
    print(f"Recovery class histogram ({region_name}):")
    for _idx, _name, _ in RECOVERY_CLASSES:
        _n = _counts.get(_idx, 0)
        _pct = 100.0 * _n / _total
        print(f"  {_name:40s}  {_n:>10,d} px  ({_pct:5.1f}%)")
except ee.EEException as _e:
    print(f"  (recovery-class histogram unavailable: {_e})")

# Plot recovery class histogram
if _hist_counts:
    _cls_idxs   = [idx for idx, _, _  in RECOVERY_CLASSES]
    _cls_colors = [col for _, _, col  in RECOVERY_CLASSES]
    _cls_vals   = [_hist_counts.get(i, 0) for i in _cls_idxs]
    _cls_total  = sum(_cls_vals) or 1

    _short_labels = [
        "Uncertain\n(no data / water)",
        "Good\nrecovery",
        "Partial\nrecovery",
        "Not\nrecovering",
        "Outside\nassessment",
    ]

    fig_h, ax_h = plt.subplots(figsize=(10, 5))
    _bars = ax_h.bar(range(len(_cls_idxs)), _cls_vals, color=_cls_colors,
                     edgecolor="white", linewidth=0.8, width=0.6)
    ax_h.set_xticks(range(len(_cls_idxs)))
    ax_h.set_xticklabels(_short_labels, fontsize=9, ha="center")
    ax_h.set_ylabel("Number of pixels", fontsize=11)
    ax_h.set_title(f"{region_name}  —  Recovery class pixel distribution",
                   fontsize=13, fontweight="bold", pad=12)
    _y_pad = max(_cls_vals) * 0.015
    for _bar, _val in zip(_bars, _cls_vals):
        _pct = 100.0 * _val / _cls_total
        ax_h.text(
            _bar.get_x() + _bar.get_width() / 2,
            _bar.get_height() + _y_pad,
            f"{_val:,}\n({_pct:.1f}%)",
            ha="center", va="bottom", fontsize=8.5, color="#222",
        )
    ax_h.set_xlim(-0.5, len(_cls_idxs) - 0.5)
    ax_h.set_ylim(0, max(_cls_vals) * 1.28)
    plt.tight_layout()
    _out_hist = FIGURES_DIR / f"histogram_recovery_class_{region_name}.png"
    fig_h.savefig(_out_hist, dpi=180, bbox_inches="tight")
    plt.close(fig_h)
    print(f"Saved → {_out_hist}")

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
        print(f"  {_name:40s}  {_n:>10,d} px  ({_pct:5.1f}%)")
except ee.EEException as _e:
    print(f"  (guard-flag histogram unavailable: {_e})")

# ── Sensitivity analysis ──────────────────────────────────────────────────────
if SENS_ENABLED:
    run_threshold_sensitivity_analysis(recent_composite, postfire_composite,
                                       prefire_composite, change_img, forest_mask_img)
else:
    print("Sensitivity analysis disabled.")

# ── Failed-pixel diagnostics ──────────────────────────────────────────────────
failed_examples  = []
failed_diag_rows = []

if FAILED_DIAG_ENABLED and not MAP_ONLY:
    import csv as _csv_diag
    print("Sampling failed recovery pixels for diagnostics...")
    failed_examples = sample_failed_locations(recovery_class)
    print(f"  sampled {len(failed_examples)} failed diagnostic locations")

    for ex in failed_examples:
        print(f"Extracting time series for {ex['id']}...")
        failed_diag_rows.extend(
            extract_failed_timeseries(ex, postfire_composite, prefire_composite)
        )

    out_diag_csv = CSV_DIR / f"failed_diagnostics_{region_name}.csv"
    if failed_diag_rows:
        with open(out_diag_csv, "w", newline="", encoding="utf-8") as _f:
            _writer = _csv_diag.DictWriter(_f, fieldnames=[
                "example_id", "cluster_rank",
                "lon", "lat", "radius_m",
                "cluster_area_m2", "cluster_pixel_count_est",
                "recovery_class",
                "date", "label",
                "NDVI", "RRI", "NDRE", "BSI", "NDMI", "NDWI", "NBR",
            ])
            _writer.writeheader()
            _writer.writerows(failed_diag_rows)
        print(f"Failed-pixel diagnostics CSV → {out_diag_csv}")
    else:
        print("No failed pixels found for diagnostics.")

# Always load cluster markers + rows from CSV if one exists from a previous run.
# This ensures CLUSTER_DATA in the JS panel is populated even when
# failed_diagnostics.enabled is false and failed_diag_rows was never extracted.
_csv_paths, _csv_points, _csv_rows = build_failed_diagnostic_plots_from_csv()
failed_diag_plot_paths = _csv_paths
failed_diag_points     = _csv_points
if not failed_diag_rows:
    failed_diag_rows = _csv_rows

# ── LST layer ─────────────────────────────────────────────────────────────────
lst_products = None
_lst_cfg = cfg.get("landsat")
if _lst_cfg:
    from landsat_lst import build_lst_layer
    print("Building Landsat LST layers…")
    try:
        lst_products = build_lst_layer(aoi, _lst_cfg)
    except Exception as _e:
        print(f"  (LST build failed: {_e})")

# ── VIS_CHANGE dict ───────────────────────────────────────────────────────────
_DIVERGE_RG    = ["#7d2222", "white", "#1f5e1f"]
_DIVERGE_GR    = ["#1f5e1f", "white", "#7d2222"]
_DIVERGE_BROWN = ["saddlebrown", "white", "steelblue"]
VIS_CHANGE = {
    "NDVI_change": {"min": -0.3, "max": 0.3, "palette": _DIVERGE_RG},
    "NDRE_change": {"min": -0.3, "max": 0.3, "palette": _DIVERGE_RG},
    "NDMI_change": {"min": -0.3, "max": 0.3, "palette": _DIVERGE_RG},
    "BSI_change":  {"min": -0.3, "max": 0.3, "palette": _DIVERGE_GR},
    "NDWI_change": {"min": -0.3, "max": 0.3, "palette": _DIVERGE_BROWN},
    "NBR_change":  {"min": -0.3, "max": 0.3, "palette": _DIVERGE_RG},
}

# ── Build Folium map ──────────────────────────────────────────────────────────
try:
    out_map = build_folium_map(
        aoi_bounds            = bounds,
        raster                = raster,
        rri_img               = rri_img,
        recovery_class        = recovery_class,
        guard_flags           = guard_flags,
        forest_mask_img       = forest_mask_img,
        prefire_composite     = prefire_composite,
        postfire_composite    = postfire_composite,
        recent_composite      = recent_composite,
        change_img            = change_img,
        VIS_CHANGE            = VIS_CHANGE,
        lst_products          = lst_products,
        failed_diag_points    = failed_diag_points,
        failed_diag_plot_paths= failed_diag_plot_paths,
        failed_diag_rows      = failed_diag_rows,
        records               = records,
        lst_cfg               = _lst_cfg,
    )
    print(f"Interactive map → {out_map}")

except ImportError:
    print("folium not available — cannot build HTML map. Install it with: pip install folium")
