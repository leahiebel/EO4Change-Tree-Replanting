"""Web-layer style, prototype tile and TiTiler manifest helpers."""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path

import ee


def default_style_config() -> dict[str, object]:
    """Return display-only style configuration for LST products."""
    return {
        "lst_target": {
            "min": 15,
            "max": 45,
            "unit": "degC",
            "palette": [
                "#313695",
                "#4575b4",
                "#74add1",
                "#abd9e9",
                "#e0f3f8",
                "#ffffbf",
                "#fee090",
                "#fdae61",
                "#f46d43",
                "#d73027",
                "#a50026",
            ],
        },
        "lst_anomaly": {
            "min": -5,
            "max": 5,
            "unit": "degC",
            "palette": [
                "#313695",
                "#4575b4",
                "#74add1",
                "#e0f3f8",
                "#ffffff",
                "#fee090",
                "#fdae61",
                "#f46d43",
                "#a50026",
            ],
        },
        "thermal_stress_class": {
            "unit": "class",
            "classes": {
                "0": {"label": "No valid data", "color": "transparent"},
                "1": {"label": "Near-normal temperature", "color": "#2ca25f"},
                "2": {"label": "Moderately warm anomaly", "color": "#feb24c"},
                "3": {"label": "Strong warm anomaly", "color": "#de2d26"},
            },
            "min": 0,
            "max": 3,
            "palette": ["#00000000", "#2ca25f", "#feb24c", "#de2d26"],
        },
    }


def write_json(path: str | Path, payload: Mapping[str, object] | list[object]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def get_prototype_tile_layer(
    image: ee.Image,
    layer_name: str,
    vis_params: Mapping[str, object],
) -> dict[str, object]:
    """Create an Earth Engine XYZ tile URL for development/testing."""
    map_id = image.getMapId(dict(vis_params))
    return {
        "layer_name": layer_name,
        "tile_url": map_id["tile_fetcher"].url_format,
        "visualization": dict(vis_params),
        "mode": "prototype",
        "warning": "Development/testing only. Do not expose private credentials in frontend code.",
    }


def build_prototype_layers(products: Mapping[str, ee.Image], styles: Mapping[str, object]) -> list[dict[str, object]]:
    """Return Earth Engine prototype tile definitions for key layers."""
    return [
        get_prototype_tile_layer(
            products["01_lst_median_target_C.tif"],
            "Seasonal median LST",
            styles["lst_target"],
        ),
        get_prototype_tile_layer(
            products["09_lst_anomaly_C.tif"],
            "LST anomaly",
            styles["lst_anomaly"],
        ),
        get_prototype_tile_layer(
            products["11_thermal_stress_class.tif"],
            "Thermal stress anomaly class",
            styles["thermal_stress_class"],
        ),
    ]


def build_cog_manifest(config: Mapping[str, object], styles: Mapping[str, object]) -> list[dict[str, object]]:
    """Build a manifest compatible with COG tile services such as TiTiler."""
    gcs_bucket = config.get("gcs_bucket")
    prefix = str(config.get("gcs_prefix") or "").strip("/")
    if gcs_bucket:
        base = f"https://storage.googleapis.com/{gcs_bucket}/{prefix}" if prefix else f"https://storage.googleapis.com/{gcs_bucket}"
    else:
        base = "replace-with-public-cog-base-url"

    titiler_base = str(config.get("titiler_base_url") or "https://tile-server.example.com/cog/tiles").rstrip("/")
    layers = [
        (
            "lst_median_2025",
            "Median land surface temperature, summer 2025",
            "01_lst_median_target_C.tif",
            "lst_target",
            "degC",
        ),
        ("lst_anomaly_2025", "LST anomaly, summer 2025", "09_lst_anomaly_C.tif", "lst_anomaly", "degC"),
        (
            "thermal_stress_class_2025",
            "Surface-temperature anomaly class, summer 2025",
            "11_thermal_stress_class.tif",
            "thermal_stress_class",
            "class",
        ),
    ]
    manifest = []
    for layer_id, title, filename, style_key, unit in layers:
        cog_url = f"{base}/{filename}"
        manifest.append(
            {
                "id": layer_id,
                "title": title,
                "type": "raster",
                "unit": unit,
                "cog_url": cog_url,
                "tile_url": f"{titiler_base}/{{z}}/{{x}}/{{y}}.png?url={cog_url}",
                "legend": styles[style_key],
                "bounds": None,
                "min_zoom": 5,
                "max_zoom": 14,
                "opacity": 0.75,
                "attribution": "USGS Landsat 8-9 Collection 2 Level-2",
            }
        )
    return manifest

