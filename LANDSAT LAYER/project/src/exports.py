"""Earth Engine batch export helpers."""

from __future__ import annotations

import logging
import time
from collections.abc import Mapping
from dataclasses import dataclass

import ee

LOGGER = logging.getLogger(__name__)

CONTINUOUS_OUTPUTS = {
    "01_lst_median_target_C.tif",
    "02_lst_mean_target_C.tif",
    "03_lst_p10_target_C.tif",
    "04_lst_p90_target_C.tif",
    "05_lst_std_target_C.tif",
    "07_lst_uncertainty_K.tif",
    "08_lst_baseline_median_C.tif",
    "09_lst_anomaly_C.tif",
    "10_lst_zscore.tif",
}
CATEGORICAL_OUTPUTS = {"06_lst_valid_observations.tif", "11_thermal_stress_class.tif"}


@dataclass(frozen=True)
class StartedTask:
    """Small export task descriptor suitable for metadata JSON."""

    description: str
    task_id: str
    destination: str
    filename: str


def _export_description(filename: str) -> str:
    return filename.replace(".tif", "").replace("-", "_")


def _prepared_image(image: ee.Image, filename: str, nodata: int | float) -> ee.Image:
    image = image.unmask(nodata)
    if filename in CATEGORICAL_OUTPUTS:
        return image.toInt16()
    return image.toFloat()


def export_image_cog(
    image: ee.Image,
    filename: str,
    aoi: ee.Geometry,
    config: Mapping[str, object],
) -> StartedTask:
    """Start a Cloud Optimized GeoTIFF image export to GCS or Drive."""
    description = _export_description(filename)
    nodata = config["nodata_value"]
    export_image = _prepared_image(image, filename, nodata)
    common = {
        "image": export_image,
        "description": description,
        "region": aoi,
        "scale": int(config["export_scale_m"]),
        "crs": str(config["output_crs"]),
        "maxPixels": 1_000_000_000_000,
        "fileFormat": "GeoTIFF",
        "formatOptions": {"cloudOptimized": True, "noData": nodata},
    }

    gcs_bucket = config.get("gcs_bucket")
    if gcs_bucket:
        prefix = str(config.get("gcs_prefix") or "").strip("/")
        file_name_prefix = f"{prefix}/{filename[:-4]}" if prefix else filename[:-4]
        task = ee.batch.Export.image.toCloudStorage(
            **common,
            bucket=str(gcs_bucket),
            fileNamePrefix=file_name_prefix,
        )
        destination = f"gs://{gcs_bucket}/{file_name_prefix}.tif"
    else:
        task = ee.batch.Export.image.toDrive(
            **common,
            folder=str(config["drive_folder"]),
            fileNamePrefix=filename[:-4],
        )
        destination = f"Google Drive/{config['drive_folder']}/{filename}"

    task.start()
    LOGGER.info("Started image export %s -> %s", description, destination)
    return StartedTask(description=description, task_id=task.id, destination=destination, filename=filename)


def export_all_rasters(
    products: Mapping[str, ee.Image],
    aoi: ee.Geometry,
    config: Mapping[str, object],
) -> list[StartedTask]:
    """Start all required raster exports."""
    tasks: list[StartedTask] = []
    for filename, image in products.items():
        tasks.append(export_image_cog(image, filename, aoi, config))
    return tasks


def export_table_csv(
    collection: ee.FeatureCollection,
    filename: str,
    config: Mapping[str, object],
) -> StartedTask:
    """Start a table export as CSV to GCS or Drive."""
    description = filename.replace(".csv", "")
    gcs_bucket = config.get("gcs_bucket")
    if gcs_bucket:
        prefix = str(config.get("gcs_prefix") or "").strip("/")
        file_name_prefix = f"{prefix}/{filename[:-4]}" if prefix else filename[:-4]
        task = ee.batch.Export.table.toCloudStorage(
            collection=collection,
            description=description,
            bucket=str(gcs_bucket),
            fileNamePrefix=file_name_prefix,
            fileFormat="CSV",
        )
        destination = f"gs://{gcs_bucket}/{file_name_prefix}.csv"
    else:
        task = ee.batch.Export.table.toDrive(
            collection=collection,
            description=description,
            folder=str(config["drive_folder"]),
            fileNamePrefix=filename[:-4],
            fileFormat="CSV",
        )
        destination = f"Google Drive/{config['drive_folder']}/{filename}"

    task.start()
    LOGGER.info("Started table export %s -> %s", description, destination)
    return StartedTask(description=description, task_id=task.id, destination=destination, filename=filename)


def wait_for_tasks(tasks: list[StartedTask], poll_seconds: int = 30) -> None:
    """Poll Earth Engine batch tasks and raise if any export fails."""
    remaining = {task.task_id: task for task in tasks}
    while remaining:
        for task_id, descriptor in list(remaining.items()):
            status = ee.data.getTaskStatus(task_id)[0]
            state = status["state"]
            if state == "COMPLETED":
                LOGGER.info("Export completed: %s", descriptor.description)
                remaining.pop(task_id)
            elif state in {"FAILED", "CANCELLED"}:
                raise RuntimeError(f"Earth Engine export failed: {descriptor.description}: {status}")
        if remaining:
            time.sleep(poll_seconds)

