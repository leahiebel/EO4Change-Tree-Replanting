"""Run the EO4C Landsat LST workflow."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import ee
import yaml

from src.authentication import initialize_earth_engine
from src.composites import (
    baseline_season_collection,
    build_all_raster_products,
    target_season_collection,
)
from src.exports import export_all_rasters, export_table_csv, wait_for_tasks
from src.forest_analysis import (
    difference_in_differences,
    forest_control_statistics,
    yearly_lst_timeseries,
)
from src.landsat_collection import build_merged_landsat_lst_collection
from src.metadata import build_processing_metadata, write_processing_metadata
from src.vector_io import load_vector_inputs, write_web_boundaries
from src.web_layers import (
    build_cog_manifest,
    build_prototype_layers,
    default_style_config,
    write_json,
)

LOGGER = logging.getLogger("eo4c_lst")


def load_config(path: str | Path) -> dict[str, Any]:
    """Load and validate the YAML configuration."""
    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    start_date = str(config["start_date"])
    end_date = str(config["end_date"])
    if start_date >= end_date:
        raise ValueError("start_date must be earlier than end_date.")
    if int(config["baseline_start_year"]) > int(config["baseline_end_year"]):
        raise ValueError("baseline_start_year must be <= baseline_end_year.")
    if not 1 <= int(config["season_start_month"]) <= 12:
        raise ValueError("season_start_month must be between 1 and 12.")
    if not 1 <= int(config["season_end_month"]) <= 12:
        raise ValueError("season_end_month must be between 1 and 12.")

    config.setdefault("minimum_baseline_images", 3)
    config.setdefault("wait_for_exports", False)
    config.setdefault("export_yearly_timeseries", True)
    return config


def _collection_size(collection: ee.ImageCollection, label: str) -> int:
    size = int(collection.size().getInfo())
    LOGGER.info("%s image count: %s", label, size)
    return size


def _assert_valid_pixel_support(valid_image: ee.Image, aoi: ee.Geometry, config: dict[str, Any]) -> None:
    coverage = valid_image.gt(0).reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=aoi,
        scale=int(config["export_scale_m"]),
        crs=str(config["output_crs"]),
        bestEffort=True,
        tileScale=4,
        maxPixels=1_000_000_000,
    ).getInfo()
    value = coverage.get("lst_valid_observations")
    if value is None or float(value) <= 0:
        raise RuntimeError("Zero valid target-season pixels after QA masking.")


def run(config_path: str | Path, authenticate: bool) -> None:
    """Run the complete LST workflow."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    config = load_config(config_path)
    output_dir = Path(config["output_directory"])
    output_dir.mkdir(parents=True, exist_ok=True)

    initialize_earth_engine(str(config["earth_engine_project"]), authenticate=authenticate)

    vectors = load_vector_inputs(config["aoi_path"], config["forest_path"], config.get("control_path"))
    write_web_boundaries(vectors, Path(__file__).parent / "web")

    landsat8_raw, landsat9_raw, collection = build_merged_landsat_lst_collection(
        vectors.aoi_geometry,
        config,
    )
    l8_count = _collection_size(landsat8_raw, "Landsat 8 L2SP")
    l9_count = _collection_size(landsat9_raw, "Landsat 9 L2SP")
    if l8_count + l9_count == 0:
        raise RuntimeError("No Landsat L2SP temperature scenes found for the configured AOI/date range.")

    target_collection = target_season_collection(collection, config)
    baseline_collection = baseline_season_collection(collection, config)
    target_count = _collection_size(target_collection, "Target seasonal")
    baseline_count = _collection_size(baseline_collection, "Baseline seasonal")
    if target_count == 0:
        raise RuntimeError("No target-year seasonal Landsat L2SP acquisitions found.")
    if baseline_count < int(config["minimum_baseline_images"]):
        raise RuntimeError(
            f"Insufficient baseline observations: {baseline_count}; "
            f"minimum is {config['minimum_baseline_images']}."
        )

    products = build_all_raster_products(collection, config)
    _assert_valid_pixel_support(products["06_lst_valid_observations.tif"], vectors.aoi_geometry, config)

    style_config = default_style_config()
    write_json(output_dir / "style_config.json", style_config)
    manifest = build_cog_manifest(config, style_config)
    write_json(output_dir / "layer_manifest.json", manifest)

    if config.get("prototype_tile_layers"):
        prototype_layers = build_prototype_layers(products, style_config)
        write_json(output_dir / "prototype_tile_layers.json", prototype_layers)

    tasks = export_all_rasters(products, vectors.aoi_geometry, config)

    forest_stats = forest_control_statistics(target_collection, vectors.forest_fc, vectors.control_fc, config)
    if vectors.control_fc is not None:
        did = difference_in_differences(collection, vectors.forest_fc, vectors.control_fc, config).map(
            lambda f: ee.Feature(f).set("analysis_type", "difference_in_differences")
        )
        forest_stats = forest_stats.map(lambda f: ee.Feature(f).set("analysis_type", "forest_control_stats")).merge(did)
    else:
        LOGGER.warning("Missing control polygons: exporting forest-only exploratory statistics.")
    tasks.append(export_table_csv(forest_stats, "forest_control_statistics.csv", config))

    if config.get("export_yearly_timeseries"):
        yearly = yearly_lst_timeseries(
            collection,
            vectors.aoi_geometry,
            vectors.forest_fc,
            vectors.control_fc,
            config,
        )
        tasks.append(export_table_csv(yearly, "yearly_lst_timeseries.csv", config))

    metadata = build_processing_metadata(
        config=config,
        aoi_gdf=vectors.aoi_gdf,
        aoi=vectors.aoi_geometry,
        landsat8_raw=landsat8_raw,
        landsat9_raw=landsat9_raw,
        valid_collection=collection,
        valid_observation_image=products["06_lst_valid_observations.tif"],
        export_tasks=[task.__dict__ for task in tasks],
    )
    write_processing_metadata(output_dir / "processing_metadata.json", metadata)

    if config.get("wait_for_exports"):
        wait_for_tasks(tasks)
    else:
        LOGGER.info("Exports started. Re-run with wait_for_exports: true to block until completion.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Landsat LST products for EO4C.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to YAML configuration.")
    parser.add_argument(
        "--no-auth",
        action="store_true",
        help="Skip ee.Authenticate() and only call ee.Initialize(project=...).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.config, authenticate=not args.no_auth)

