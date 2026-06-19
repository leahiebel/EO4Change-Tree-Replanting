from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import LineString, Polygon

from src.vector_io import validate_vector_layer


def _write_layer(path: Path, geometry) -> None:
    gdf = gpd.GeoDataFrame({"site_name": ["test"]}, geometry=[geometry], crs="EPSG:4326")
    gdf.to_file(path, driver="GeoJSON")


def test_validate_vector_layer_adds_feature_id(tmp_path: Path) -> None:
    path = tmp_path / "polygon.geojson"
    _write_layer(path, Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]))
    result = validate_vector_layer(path, "test")
    assert "feature_id" in result.columns
    assert result.crs.to_string() == "EPSG:4326"


def test_validate_vector_layer_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        validate_vector_layer(tmp_path / "missing.geojson", "missing")


def test_validate_vector_layer_rejects_non_polygon(tmp_path: Path) -> None:
    path = tmp_path / "line.geojson"
    _write_layer(path, LineString([(0, 0), (1, 1)]))
    with pytest.raises(ValueError, match="polygons or multipolygons"):
        validate_vector_layer(path, "line")

