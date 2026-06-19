"""Local vector validation and conversion to Earth Engine objects."""

from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import ee
import geopandas as gpd
from shapely.geometry.base import BaseGeometry

LOGGER = logging.getLogger(__name__)

POLYGON_TYPES = {"Polygon", "MultiPolygon"}
IDENTIFIER_FIELDS = ["project_id", "site_name", "planting_year", "species", "management_type"]


@dataclass(frozen=True)
class VectorInputs:
    """Validated local and Earth Engine vector inputs."""

    aoi_gdf: gpd.GeoDataFrame
    forest_gdf: gpd.GeoDataFrame
    control_gdf: gpd.GeoDataFrame | None
    aoi_geometry: ee.Geometry
    forest_fc: ee.FeatureCollection
    control_fc: ee.FeatureCollection | None


def _repair_geometry(geometry: BaseGeometry) -> BaseGeometry:
    if geometry.is_valid:
        return geometry
    try:
        repaired = geometry.make_valid()
    except AttributeError:
        repaired = geometry.buffer(0)
    if repaired.is_empty or not repaired.is_valid:
        raise ValueError("Geometry could not be repaired safely.")
    return repaired


def _stable_geometry_id(geometry: BaseGeometry, index: int) -> str:
    digest = hashlib.sha1(geometry.wkb).hexdigest()[:12]
    return f"feature_{index:05d}_{digest}"


def _ensure_identifier(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gdf = gdf.copy()
    if "feature_id" not in gdf.columns:
        gdf["feature_id"] = [_stable_geometry_id(geom, idx) for idx, geom in enumerate(gdf.geometry)]
    return gdf


def _select_properties(row: object, columns: Iterable[str]) -> dict[str, object]:
    props: dict[str, object] = {}
    for column in columns:
        if column == "geometry":
            continue
        value = getattr(row, column)
        if value is not None:
            props[column] = value
    return props


def _feature_collection_from_gdf(gdf: gpd.GeoDataFrame) -> ee.FeatureCollection:
    features = []
    property_columns = [c for c in gdf.columns if c != "geometry"]
    for row in gdf.itertuples(index=False):
        geom = getattr(row, "geometry")
        props = _select_properties(row, property_columns)
        features.append(ee.Feature(ee.Geometry(geom.__geo_interface__), props))
    return ee.FeatureCollection(features)


def validate_vector_layer(path: str | Path, layer_name: str) -> gpd.GeoDataFrame:
    """Read and validate a polygon vector layer."""
    path = Path(path).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"{layer_name} file does not exist: {path}")

    gdf = gpd.read_file(path)
    if gdf.empty:
        raise ValueError(f"{layer_name} contains no features.")
    if gdf.geometry.isna().any() or gdf.geometry.is_empty.any():
        raise ValueError(f"{layer_name} contains empty geometries.")

    gdf = gdf.copy()
    gdf["geometry"] = gdf.geometry.apply(_repair_geometry)
    geom_types = set(gdf.geometry.geom_type)
    if not geom_types.issubset(POLYGON_TYPES):
        raise ValueError(f"{layer_name} must contain only polygons or multipolygons, found {geom_types}.")

    if gdf.crs is None:
        raise ValueError(f"{layer_name} has no CRS. Assign a CRS before running the workflow.")

    gdf = gdf.drop_duplicates(subset="geometry").to_crs("EPSG:4326")
    gdf = _ensure_identifier(gdf)

    for field in IDENTIFIER_FIELDS:
        if field not in gdf.columns:
            gdf[field] = None

    return gdf


def _assert_intersects_aoi(layer: gpd.GeoDataFrame, aoi: gpd.GeoDataFrame, layer_name: str) -> None:
    aoi_union = aoi.geometry.union_all()
    intersects = layer.geometry.intersects(aoi_union)
    if not bool(intersects.any()):
        raise ValueError(f"{layer_name} does not intersect the AOI.")


def load_vector_inputs(
    aoi_path: str | Path,
    forest_path: str | Path,
    control_path: str | Path | None,
) -> VectorInputs:
    """Load, validate and convert AOI, forest and optional control vectors."""
    aoi_gdf = validate_vector_layer(aoi_path, "AOI")
    forest_gdf = validate_vector_layer(forest_path, "Forest polygons")
    _assert_intersects_aoi(forest_gdf, aoi_gdf, "Forest polygons")

    control_gdf = None
    control_fc = None
    if control_path:
        control_gdf = validate_vector_layer(control_path, "Control polygons")
        _assert_intersects_aoi(control_gdf, aoi_gdf, "Control polygons")
        control_fc = _feature_collection_from_gdf(control_gdf)
    else:
        LOGGER.warning(
            "No control polygons configured. Forest cooling analysis will be skipped or labelled exploratory."
        )

    aoi_geometry = ee.Geometry(aoi_gdf.geometry.union_all().__geo_interface__)
    return VectorInputs(
        aoi_gdf=aoi_gdf,
        forest_gdf=forest_gdf,
        control_gdf=control_gdf,
        aoi_geometry=aoi_geometry,
        forest_fc=_feature_collection_from_gdf(forest_gdf),
        control_fc=control_fc,
    )


def write_web_boundaries(vectors: VectorInputs, web_directory: str | Path) -> None:
    """Write AOI and forest boundaries for the Leaflet prototype."""
    web_directory = Path(web_directory)
    web_directory.mkdir(parents=True, exist_ok=True)
    vectors.aoi_gdf.to_file(web_directory / "aoi.geojson", driver="GeoJSON")
    vectors.forest_gdf.to_file(web_directory / "forest.geojson", driver="GeoJSON")

