"""
make_config.py — generate config_<name>.yaml + area_<name>.geojson from a point.

Given a longitude/latitude inside a replanted area, build an AOI polygon and
emit a config file that gee.py can consume.

Two polygon-extraction methods:
  --method bbox   simple square buffer around the seed point (no GEE call)
  --method snic   SNIC superpixel from a Sentinel-2 summer median composite,
                  picking the cluster that contains the seed point.
                  Requires GEE auth (ee.Authenticate() run once).

Examples
--------
  # Rugballegård (matches the existing config_DK.yaml AOI):
  python make_config.py --lon 9.795 --lat 55.866 --name rugballegaard-skov

  # Coarse bbox only, no GEE call:
  python make_config.py --lon 12.34 --lat 56.78 --name new-site --method bbox --buffer 600
"""
import argparse
import json
import math
from copy import deepcopy
from pathlib import Path

import ee
import yaml


# ── Config skeleton (mirrors config_DK.yaml) ─────────────────────────────────
CONFIG_TEMPLATE: dict = {
    "spatial": {
        "type": "geojson",
        "geojson_path": None,            # filled in per-region
        "geojson_clip": True,
        "region_name": None,             # filled in per-region
    },
    "temporal": {
        "start": "2018-01-02",
        "end":   "2024-01-01",
        "cadence": {"type": "fixed", "interval": "monthly"},
    },
    "variables": {
        "model": "s2biophys",
        "variable": ["laie", "fapar", "fcover"],
    },
    "export": {
        "destination": "drive",
        "collection_path": None,         # filled in: "EO4Change/<name>"
        "project_id": "eo4change-ex",
        "crs": "EPSG:32632",
        "scale": 100,
        "max_pixels": 100_000_000_000,
    },
    "options": {
        "max_cloud_cover": 50,
        "csplus_band": "cs",
        "cs_plus_threshold": 0.70,
        "clip_min_max": True,
    },
    "version": "v02",
}


# ── Method A: plain buffered bbox ────────────────────────────────────────────
def polygon_bbox(lon: float, lat: float, buffer_m: float) -> dict:
    """Square AOI of half-side `buffer_m` around (lon, lat). Returns GeoJSON geometry."""
    dlat = buffer_m / 111_320.0
    dlon = buffer_m / (111_320.0 * max(math.cos(math.radians(lat)), 1e-6))
    ring = [
        [lon - dlon, lat - dlat],
        [lon + dlon, lat - dlat],
        [lon + dlon, lat + dlat],
        [lon - dlon, lat + dlat],
        [lon - dlon, lat - dlat],
    ]
    return {"type": "Polygon", "coordinates": [ring]}


# ── Method B: SNIC superpixel containing the seed point ──────────────────────
def polygon_snic(
    lon: float, lat: float, buffer_m: float,
    project_id: str,
    year: str = "2023",
    snic_size: int = 10,
    snic_compactness: float = 0.1,
    simplify_error_m: float = 5.0,
) -> dict:
    """Extract the SNIC superpixel containing (lon, lat) from a S2 summer median.

    The seed lon/lat must fall on a pixel that has at least one cloud-free
    Sentinel-2 scene in the summer window of `year`.
    """
    ee.Initialize(project=project_id)

    seed   = ee.Geometry.Point([lon, lat])
    region = seed.buffer(buffer_m).bounds()

    s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(region)
            .filterDate(f"{year}-06-01", f"{year}-09-01")
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30)))
    composite = (s2.median()
                   .select(["B2", "B3", "B4", "B8", "B11"])
                   .clip(region))

    snic = ee.Algorithms.Image.Segmentation.SNIC(
        image=composite,
        size=snic_size,
        compactness=snic_compactness,
        connectivity=8,
        neighborhoodSize=64,
    )
    clusters = snic.select("clusters")

    # Round-trip the cluster ID at the seed so we can give a clear error if missing
    seed_id_val = clusters.reduceRegion(
        reducer=ee.Reducer.first(), geometry=seed, scale=10,
    ).get("clusters").getInfo()
    if seed_id_val is None:
        raise RuntimeError(
            f"No cloud-free Sentinel-2 pixel at ({lon}, {lat}) in summer {year}. "
            f"Try a different --year or use --method bbox."
        )

    mask = clusters.eq(ee.Number(seed_id_val)).selfMask()
    fc = mask.reduceToVectors(
        geometry=region, scale=10, eightConnected=True,
        geometryType="polygon", maxPixels=int(1e9),
    )
    # Safety net: pick the polygon that actually contains the seed
    feature = ee.Feature(fc.filterBounds(seed).first())
    geom = feature.geometry().simplify(maxError=simplify_error_m)
    return geom.getInfo()


# ── Output writers ───────────────────────────────────────────────────────────
def write_geojson(path: Path, geometry: dict, seed: tuple[float, float], name: str) -> None:
    fc = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": {"name": name, "seed_lonlat": list(seed)},
            "geometry": geometry,
        }],
    }
    path.write_text(json.dumps(fc))


def write_config(path: Path, name: str, geojson_filename: str, project_id: str) -> None:
    cfg = deepcopy(CONFIG_TEMPLATE)
    cfg["spatial"]["geojson_path"]   = geojson_filename
    cfg["spatial"]["region_name"]    = name
    cfg["export"]["collection_path"] = f"EO4Change/{name}"
    cfg["export"]["project_id"]      = project_id
    header = (f"# Auto-generated by make_config.py for region '{name}'.\n"
              f"# Edit temporal / variables / export blocks as needed.\n\n")
    with open(path, "w") as f:
        f.write(header)
        yaml.safe_dump(cfg, f, sort_keys=False)


# ── CLI ──────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--lon",  type=float, required=True, help="seed longitude (deg)")
    ap.add_argument("--lat",  type=float, required=True, help="seed latitude (deg)")
    ap.add_argument("--name", required=True, help="region name (no spaces; hyphens ok)")
    ap.add_argument("--method", choices=["bbox", "snic"], default="snic",
                    help="polygon extraction method (default: snic)")
    ap.add_argument("--buffer", type=float, default=500.0,
                    help="half-side search/buffer radius around seed, metres (default 500)")
    ap.add_argument("--year", default="2023",
                    help="year of the SNIC summer composite (default 2023)")
    ap.add_argument("--snic-size", type=int, default=10,
                    help="SNIC seed spacing in pixels (default 10)")
    ap.add_argument("--snic-compactness", type=float, default=0.1,
                    help="SNIC compactness (default 0.1)")
    ap.add_argument("--project", default=CONFIG_TEMPLATE["export"]["project_id"],
                    help="GEE project id (default %(default)s)")
    ap.add_argument("--out-dir", default=None,
                    help="output directory (default: same dir as this script)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir) if args.out_dir else Path(__file__).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.method == "bbox":
        geometry = polygon_bbox(args.lon, args.lat, args.buffer)
    else:
        print(f"→ SNIC around ({args.lon}, {args.lat}) on {args.year}-summer S2 composite "
              f"(buffer {args.buffer:.0f} m, size {args.snic_size}, "
              f"compactness {args.snic_compactness})…")
        geometry = polygon_snic(
            args.lon, args.lat, args.buffer, args.project,
            year=args.year,
            snic_size=args.snic_size,
            snic_compactness=args.snic_compactness,
        )

    gj_path  = out_dir / f"area_{args.name}.geojson"
    cfg_path = out_dir / f"config_{args.name}.yaml"
    write_geojson(gj_path, geometry, (args.lon, args.lat), args.name)
    write_config(cfg_path, args.name, gj_path.name, args.project)

    print(f"Wrote {gj_path}")
    print(f"Wrote {cfg_path}")
    print(
        "\nNext: point gee.py at the new config — edit CONFIG_PATH in gee.py, "
        f"or rename {cfg_path.name} to config_DK.yaml."
    )


if __name__ == "__main__":
    main()
