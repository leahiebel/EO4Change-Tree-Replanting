"""
EO4Change Group 4 — failed-pixel diagnostics: sampling, time-series extraction, plot generation.
"""
import csv
from pathlib import Path
import ee
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from gee_config import (SRC_DIR, FIGURES_DIR, CSV_DIR,
                         region_name, exp, aoi, windows,
                         FAILED_DIAG_N, FAILED_DIAG_RADIUS_M,
                         T_RRI_LOW, T_BSI_HIGH, PLANTING_DATE)
from gee_s2 import make_composite
from gee_recovery import compute_rri


def _fmt(v, digits=3):
    """Format numbers nicely for HTML tables."""
    if v is None:
        return ""
    try:
        return f"{float(v):.{digits}f}"
    except Exception:
        return str(v)


def sample_failed_locations(recovery_class: ee.Image) -> list[dict]:
    """Select diagnostic locations from the largest connected clusters of recovery class 3."""
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


def extract_failed_timeseries(example: dict, postfire_composite: ee.Image,
                               prefire_composite: ee.Image) -> list[dict]:
    """Extract monthly mean spectral indices around one failed location. RRI recomputed per window."""
    rows = []
    geom = example["geometry"]
    for ws, we, label in windows:
        try:
            current = make_composite(ws, we)
            rri_t = compute_rri(recent=current, postfire=postfire_composite, prefire=prefire_composite)
            img = (current
                   .select(["NDVI", "NDRE", "NDMI", "BSI", "NDWI", "NBR"])
                   .addBands(rri_t))
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


def make_failed_popup_table(example_id: str, failed_diag_rows: list[dict]) -> str:
    """Build a compact scrollable HTML table for the Folium popup."""
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


def build_failed_diagnostic_plots_from_csv() -> tuple[dict, list, list]:
    """Read failed_diagnostics_<region>.csv, generate per-cluster PNGs.

    Returns (plot_paths, points, rows) where rows is the full time-series
    list needed to populate CLUSTER_DATA in the Folium JS panel.
    """
    csv_path = CSV_DIR / f"failed_diagnostics_{region_name}.csv"
    failed_diag_plot_paths: dict = {}
    failed_diag_points: list = []

    if not csv_path.exists():
        print(f"No failed diagnostics CSV found: {csv_path}")
        return failed_diag_plot_paths, failed_diag_points, []

    df = pd.read_csv(csv_path)
    if "recovery_class" in df.columns:
        df = df[df["recovery_class"].astype(int) == 3]
    if "guard_flag" in df.columns:
        df = df[df["guard_flag"].astype(int) != 4]

    if "cluster_rank" in df.columns:
        keep_ids = (
            df[["example_id", "cluster_rank"]]
            .drop_duplicates()
            .sort_values("cluster_rank")
            .head(FAILED_DIAG_N)["example_id"]
            .tolist()
        )
    else:
        keep_ids = df["example_id"].drop_duplicates().head(FAILED_DIAG_N).tolist()

    df = df[df["example_id"].isin(keep_ids)]
    if df.empty:
        print("No failed diagnostic points left after filtering.")
        return failed_diag_plot_paths, failed_diag_points, []

    df["date"] = pd.to_datetime(df["date"])

    for ex_id, sub in df.groupby("example_id"):
        sub = sub.sort_values("date").copy()
        lon = float(sub["lon"].iloc[0])
        lat = float(sub["lat"].iloc[0])
        radius_m = int(sub["radius_m"].iloc[0]) if "radius_m" in sub.columns else FAILED_DIAG_RADIUS_M

        fig, axes = plt.subplots(3, 1, figsize=(8.5, 7), sharex=True)

        axes[0].plot(sub["date"], sub["NDVI"], marker="o", linewidth=1.5, label="NDVI")
        axes[0].plot(sub["date"], sub["RRI"], marker="s", linestyle="--", linewidth=1.5, label="RRI")
        axes[0].plot(sub["date"], sub["NDRE"], marker="o", linewidth=1.2, label="NDRE")
        axes[0].axhline(T_RRI_LOW, linestyle=":", linewidth=1, label=f"RRI failed threshold ({T_RRI_LOW})")
        axes[0].set_ylabel("Recovery / canopy")
        axes[0].legend(fontsize=8, ncol=2)
        axes[0].grid(alpha=0.25)

        axes[1].plot(sub["date"], sub["BSI"], marker="o", linewidth=1.5, label="BSI")
        axes[1].plot(sub["date"], sub["NDMI"], marker="o", linewidth=1.5, label="NDMI")
        axes[1].plot(sub["date"], sub["NDWI"], marker="o", linewidth=1.2, label="NDWI")
        axes[1].axhline(T_BSI_HIGH, linestyle=":", linewidth=1, label=f"High BSI threshold ({T_BSI_HIGH})")
        axes[1].set_ylabel("Soil / moisture")
        axes[1].legend(fontsize=8, ncol=2)
        axes[1].grid(alpha=0.25)

        axes[2].plot(sub["date"], sub["NBR"], marker="o", linewidth=1.5, label="NBR")
        axes[2].set_ylabel("NBR")
        axes[2].legend(fontsize=8)
        axes[2].grid(alpha=0.25)

        if PLANTING_DATE is not None:
            for ax in axes:
                ax.axvline(PLANTING_DATE, linestyle="--", linewidth=1)

        axes[2].xaxis.set_major_locator(mdates.YearLocator())
        axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        fig.suptitle(
            f"{ex_id} — failed recovery diagnostic\n"
            f"lon={lon:.5f}, lat={lat:.5f}, radius={radius_m} m",
            fontsize=11,
        )
        plt.tight_layout()
        out_png = FIGURES_DIR / f"{ex_id}_{region_name}_diagnostic.png"
        fig.savefig(out_png, dpi=140)
        plt.close(fig)

        failed_diag_plot_paths[ex_id] = out_png
        failed_diag_points.append({"id": ex_id, "lon": lon, "lat": lat, "radius_m": radius_m})
        print(f"Saved diagnostic plot → {out_png}")

    # Convert DataFrame rows → list of dicts for CLUSTER_DATA in the JS panel.
    # Dates are Timestamps after pd.to_datetime; convert back to "YYYY-MM-DD" strings.
    _ts_cols = ["NDVI", "RRI", "NDRE", "BSI", "NDMI", "NDWI", "NBR"]
    failed_diag_rows_out: list = []
    for _, r in df.sort_values(["example_id", "date"]).iterrows():
        row_dict: dict = {
            "example_id": str(r["example_id"]),
            "cluster_rank": int(r["cluster_rank"]) if "cluster_rank" in r.index and pd.notna(r["cluster_rank"]) else None,
            "date": r["date"].strftime("%Y-%m-%d"),
        }
        for _k in _ts_cols:
            _v = r.get(_k)
            row_dict[_k] = float(_v) if (_v is not None and pd.notna(_v)) else None
        failed_diag_rows_out.append(row_dict)

    return failed_diag_plot_paths, failed_diag_points, failed_diag_rows_out
