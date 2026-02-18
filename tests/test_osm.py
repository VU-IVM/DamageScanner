"""Tests for damagescanner.osm module."""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import (
    GeometryCollection,
    MultiPolygon,
    Point,
    Polygon,
)

from damagescanner.osm import (
    DICT_CIS_OSM,
    _combine_columns,
    _extract_value,
    _filter_dataframe,
    _remove_contained_assets,
    _remove_contained_points,
    _remove_contained_polys,
    create_point_from_polygon,
    extract_first_geom,
    read_osm_data,
)

from .helpers import data_path


class TestCombineColumns:
    """Tests for _combine_columns function."""

    def test_only_first_value(self) -> None:
        """Test when only first value is present."""
        result = _combine_columns("highway", None)
        assert result == "highway"

    def test_only_second_value(self) -> None:
        """Test when only second value is present."""
        result = _combine_columns(None, "residential")
        assert result == "residential"

    def test_both_values_same(self) -> None:
        """Test when both values are the same."""
        result = _combine_columns("school", "school")
        assert result == "school"

    def test_both_values_different(self) -> None:
        """Test when both values are different."""
        result = _combine_columns("primary", "secondary")
        assert result == "primary"  # First value takes precedence

    def test_first_value_is_yes(self) -> None:
        """Test when first value is 'yes'."""
        result = _combine_columns("yes", "residential")
        assert result == "residential"

    def test_second_value_is_yes(self) -> None:
        """Test when second value is 'yes'."""
        result = _combine_columns("residential", "yes")
        assert result == "residential"

    def test_both_none(self) -> None:
        """Test when both values are None."""
        result = _combine_columns(None, None)
        assert result is None

    def test_with_nan(self) -> None:
        """Test with NaN values."""
        result = _combine_columns(np.nan, "value")
        assert result == "value"


class TestFilterDataframe:
    """Tests for _filter_dataframe function."""

    def test_filter_with_two_columns(self) -> None:
        """Test filtering with two columns."""
        gdf = gpd.GeoDataFrame(
            {
                "col1": ["highway", "building", None],
                "col2": [None, "residential", "school"],
                "geometry": [Point(0, 0), Point(1, 1), Point(2, 2)],
            }
        )

        result = _filter_dataframe(gdf, ["col1", "col2"])

        assert "object_type" in result.columns
        assert "col1" not in result.columns
        assert "col2" not in result.columns

    def test_filter_with_three_columns(self) -> None:
        """Test filtering with three columns."""
        gdf = gpd.GeoDataFrame(
            {
                "col1": ["highway", None],
                "col2": [None, "building"],
                "col3": ["primary", "residential"],
                "geometry": [Point(0, 0), Point(1, 1)],
            }
        )

        result = _filter_dataframe(gdf, ["col1", "col2", "col3"])

        assert "object_type" in result.columns
        assert "col1" not in result.columns
        assert "col2" not in result.columns
        assert "col3" not in result.columns


class TestExtractFirstGeom:
    """Tests for extract_first_geom function."""

    def test_geometry_collection(self) -> None:
        """Test extracting from GeometryCollection."""
        point = Point(0, 0)
        polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        gc = GeometryCollection([point, polygon])

        result = extract_first_geom(gc)

        assert result == point

    def test_empty_geometry_collection(self) -> None:
        """Test with empty GeometryCollection."""
        gc = GeometryCollection([])
        result = extract_first_geom(gc)
        assert result == gc  # Returns unchanged

    def test_regular_geometry(self) -> None:
        """Test with regular geometry (not collection)."""
        point = Point(0, 0)
        result = extract_first_geom(point)
        assert result == point

    def test_polygon(self) -> None:
        """Test with polygon geometry."""
        polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        result = extract_first_geom(polygon)
        assert result == polygon


class TestRemoveContainedPoints:
    """Tests for _remove_contained_points function."""

    def test_point_inside_polygon_removed(self) -> None:
        """Test that points inside polygons are removed."""
        polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        point_inside = Point(5, 5)
        point_outside = Point(20, 20)

        gdf = gpd.GeoDataFrame(
            {
                "id": [1, 2, 3],
                "geometry": [polygon, point_inside, point_outside],
            }
        )

        result = _remove_contained_points(gdf)

        assert len(result) == 2
        assert Point(20, 20) in result.geometry.values
        assert Point(5, 5) not in result.geometry.values

    def test_no_points_inside(self) -> None:
        """Test when no points are inside polygons."""
        polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        point = Point(10, 10)

        gdf = gpd.GeoDataFrame(
            {
                "id": [1, 2],
                "geometry": [polygon, point],
            }
        )

        result = _remove_contained_points(gdf)

        assert len(result) == 2


class TestRemoveContainedPolys:
    """Tests for _remove_contained_polys function."""

    def test_smaller_polygon_inside_larger_removed(self) -> None:
        """Test that smaller polygons inside larger ones are removed."""
        large_polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        small_polygon = Polygon([(2, 2), (4, 2), (4, 4), (2, 4)])
        outside_polygon = Polygon([(20, 20), (25, 20), (25, 25), (20, 25)])

        gdf = gpd.GeoDataFrame(
            {
                "id": [1, 2, 3],
                "geometry": [large_polygon, small_polygon, outside_polygon],
            }
        )

        result = _remove_contained_polys(gdf)

        assert len(result) == 2

    def test_no_contained_polygons(self) -> None:
        """Test when no polygons are contained."""
        polygon1 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        polygon2 = Polygon([(5, 5), (6, 5), (6, 6), (5, 6)])

        gdf = gpd.GeoDataFrame(
            {
                "id": [1, 2],
                "geometry": [polygon1, polygon2],
            }
        )

        result = _remove_contained_polys(gdf)

        assert len(result) == 2


class TestRemoveContainedAssets:
    """Tests for _remove_contained_assets function."""

    def test_removes_points_and_polygons(self) -> None:
        """Test that both contained points and polygons are removed."""
        large_polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        small_polygon = Polygon([(2, 2), (4, 2), (4, 4), (2, 4)])
        point_inside = Point(5, 5)
        point_outside = Point(20, 20)

        gdf = gpd.GeoDataFrame(
            {
                "id": [1, 2, 3, 4],
                "geometry": [large_polygon, small_polygon, point_inside, point_outside],
            }
        )

        result = _remove_contained_assets(gdf)

        # Should keep large_polygon and point_outside
        assert len(result) == 2


class TestCreatePointFromPolygon:
    """Tests for create_point_from_polygon function."""

    def test_polygon_converted_to_centroid(self) -> None:
        """Test that polygons are converted to their centroids."""
        polygon = Polygon([(0, 0), (4, 0), (4, 4), (0, 4)])

        gdf = gpd.GeoDataFrame(
            {
                "id": [1],
                "geometry": [polygon],
            }
        )

        result = create_point_from_polygon(gdf)

        # Centroid of the square should be (2, 2)
        assert result.geometry.iloc[0].geom_type == "Point"
        assert result.geometry.iloc[0].x == 2.0
        assert result.geometry.iloc[0].y == 2.0

    def test_multipolygon_converted(self) -> None:
        """Test that multipolygons are converted to centroids."""
        poly1 = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
        poly2 = Polygon([(4, 4), (6, 4), (6, 6), (4, 6)])
        multipoly = MultiPolygon([poly1, poly2])

        gdf = gpd.GeoDataFrame(
            {
                "id": [1],
                "geometry": [multipoly],
            }
        )

        result = create_point_from_polygon(gdf)

        assert result.geometry.iloc[0].geom_type == "Point"


class TestExtractValue:
    """Tests for _extract_value function."""

    def test_extract_existing_key(self) -> None:
        """Test extracting an existing key."""
        text = '"highway"=>"primary","name"=>"Main Street"'
        result = _extract_value(text, "highway")
        assert result == "primary"

    def test_extract_name(self) -> None:
        """Test extracting name key."""
        text = '"highway"=>"primary","name"=>"Main Street"'
        result = _extract_value(text, "name")
        assert result == "Main Street"

    def test_extract_missing_key(self) -> None:
        """Test extracting a missing key returns None."""
        text = '"highway"=>"primary"'
        result = _extract_value(text, "lanes")
        assert result is None

    def test_extract_from_none(self) -> None:
        """Test extracting from None returns None."""
        result = _extract_value(None, "highway")
        assert result is None

    def test_extract_complex_value(self) -> None:
        """Test extracting complex value with special characters."""
        text = '"voltage"=>"110000","operator"=>"National Grid"'
        result = _extract_value(text, "voltage")
        assert result == "110000"


class TestDictCisOsm:
    """Tests for DICT_CIS_OSM configuration."""

    def test_all_asset_types_have_required_keys(self) -> None:
        """Test that all asset types have osm_keys and osm_query."""
        for asset_type, config in DICT_CIS_OSM.items():
            assert "osm_keys" in config, f"{asset_type} missing osm_keys"
            assert "osm_query" in config, f"{asset_type} missing osm_query"

    def test_osm_keys_are_lists(self) -> None:
        """Test that osm_keys are lists."""
        for asset_type, config in DICT_CIS_OSM.items():
            assert isinstance(config["osm_keys"], list), (
                f"{asset_type} osm_keys not a list"
            )

    def test_osm_query_are_dicts(self) -> None:
        """Test that osm_query are dicts."""
        for asset_type, config in DICT_CIS_OSM.items():
            assert isinstance(config["osm_query"], dict), (
                f"{asset_type} osm_query not a dict"
            )

    def test_expected_asset_types_exist(self) -> None:
        """Test that expected asset types exist."""
        expected_types = [
            "roads",
            "main_roads",
            "rail",
            "air",
            "telecom",
            "water_supply",
            "waste_solid",
            "waste_water",
            "education",
            "healthcare",
            "power",
            "gas",
            "oil",
            "buildings",
        ]
        for asset_type in expected_types:
            assert asset_type in DICT_CIS_OSM, f"{asset_type} not in DICT_CIS_OSM"


# Tests requiring OSM data file
class TestReadOsmData:
    """Tests for read_osm_data function using Jamaica OSM data."""

    @pytest.fixture
    def osm_file(self) -> Path:
        """Fixture for OSM file path.

        Returns:
            Path to Jamaica OSM data file, or skips tests if not available.
        """
        jamaica_path = data_path / "jamaica" / "exposure" / "jamaica-latest.osm.pbf"
        if not jamaica_path.exists():
            pytest.skip("OSM test data not available")
        return jamaica_path

    def test_read_osm_data_returns_geodataframe(self, osm_file: Path) -> None:
        """Test that read_osm_data returns a GeoDataFrame."""
        result = read_osm_data(osm_file, "main_roads")
        assert isinstance(result, gpd.GeoDataFrame)

    def test_read_osm_data_has_required_columns(self, osm_file: Path) -> None:
        """Test that result has required columns."""
        result = read_osm_data(osm_file, "main_roads")

        assert "osm_id" in result.columns
        assert "geometry" in result.columns
        assert "object_type" in result.columns

    def test_read_osm_data_invalid_asset_type(self, osm_file: Path) -> None:
        """Test that invalid asset type raises warning."""
        with pytest.raises(ImportWarning):
            read_osm_data(osm_file, "invalid_asset_type")

    def test_read_osm_data_geometries_valid(self, osm_file: Path) -> None:
        """Test that all geometries are valid."""
        result = read_osm_data(osm_file, "main_roads")

        if len(result) > 0:
            assert result.geometry.is_valid.all()

    def test_read_osm_main_roads(self, osm_file: Path) -> None:
        """Test extraction of main roads."""
        result = read_osm_data(osm_file, "main_roads")

        assert len(result) > 0
        assert all(result.geometry.geom_type.isin(["LineString", "MultiLineString"]))

    def test_read_osm_roads(self, osm_file: Path) -> None:
        """Test extraction of all roads."""
        result = read_osm_data(osm_file, "roads")

        assert len(result) > 0
        # Roads should have more features than main_roads
        main_roads = read_osm_data(osm_file, "main_roads")
        assert len(result) >= len(main_roads)

    def test_read_osm_buildings(self, osm_file: Path) -> None:
        """Test extraction of buildings."""
        result = read_osm_data(osm_file, "buildings")

        assert len(result) > 0
        # Buildings can be points, polygons, or lines (unclosed outlines)
        valid_types = [
            "Point",
            "Polygon",
            "MultiPolygon",
            "LineString",
            "MultiLineString",
        ]
        assert all(result.geometry.geom_type.isin(valid_types))

    def test_read_osm_healthcare(self, osm_file: Path) -> None:
        """Test extraction of healthcare facilities."""
        result = read_osm_data(osm_file, "healthcare")

        # May or may not have healthcare facilities
        assert isinstance(result, gpd.GeoDataFrame)
        if len(result) > 0:
            valid_types = ["Point", "Polygon", "MultiPolygon"]
            assert all(result.geometry.geom_type.isin(valid_types))

    def test_read_osm_education(self, osm_file: Path) -> None:
        """Test extraction of education facilities."""
        result = read_osm_data(osm_file, "education")

        assert isinstance(result, gpd.GeoDataFrame)
        if len(result) > 0:
            valid_types = ["Point", "Polygon", "MultiPolygon"]
            assert all(result.geometry.geom_type.isin(valid_types))

    def test_read_osm_power(self, osm_file: Path) -> None:
        """Test extraction of power infrastructure."""
        result = read_osm_data(osm_file, "power")

        assert isinstance(result, gpd.GeoDataFrame)
        if len(result) > 0:
            # Power can be points, lines, or polygons
            valid_types = [
                "Point",
                "LineString",
                "MultiLineString",
                "Polygon",
                "MultiPolygon",
            ]
            assert all(result.geometry.geom_type.isin(valid_types))

    def test_read_osm_rail(self, osm_file: Path) -> None:
        """Test extraction of rail infrastructure."""
        result = read_osm_data(osm_file, "rail")

        assert isinstance(result, gpd.GeoDataFrame)
        if len(result) > 0:
            assert all(
                result.geometry.geom_type.isin(["LineString", "MultiLineString"])
            )

    def test_read_osm_telecom(self, osm_file: Path) -> None:
        """Test extraction of telecom infrastructure."""
        result = read_osm_data(osm_file, "telecom")

        assert isinstance(result, gpd.GeoDataFrame)
        if len(result) > 0:
            valid_types = ["Point", "Polygon", "MultiPolygon"]
            assert all(result.geometry.geom_type.isin(valid_types))

    def test_read_osm_water_supply(self, osm_file: Path) -> None:
        """Test extraction of water supply infrastructure."""
        result = read_osm_data(osm_file, "water_supply")

        assert isinstance(result, gpd.GeoDataFrame)

    def test_read_osm_no_geometry_collections(self, osm_file: Path) -> None:
        """Test that no GeometryCollections remain after extraction."""
        result = read_osm_data(osm_file, "main_roads")

        if len(result) > 0:
            assert not any(result.geometry.geom_type == "GeometryCollection")

    def test_read_osm_unique_geometries(self, osm_file: Path) -> None:
        """Test that contained geometries are removed."""
        result = read_osm_data(osm_file, "buildings")

        if len(result) > 0:
            # Check that result has been processed for duplicates
            # (index should be reset)
            assert result.index.is_unique

    def test_read_osm_object_types_in_vulnerability_dict(self, osm_file: Path) -> None:
        """Test that object types are filtered to vulnerability dict."""
        from damagescanner.config import DICT_CIS_VULNERABILITY_FLOOD

        result = read_osm_data(osm_file, "main_roads")

        if len(result) > 0:
            valid_types = list(DICT_CIS_VULNERABILITY_FLOOD["main_roads"].keys())
            assert all(result["object_type"].isin(valid_types))


class TestExtract:
    """Tests for extract function using Jamaica OSM data."""

    @pytest.fixture
    def osm_file(self) -> Path:
        """Fixture for OSM file path.

        Returns:
            Path to Jamaica OSM data file, or skips tests if not available.
        """
        jamaica_path = data_path / "jamaica" / "exposure" / "jamaica-latest.osm.pbf"
        if not jamaica_path.exists():
            pytest.skip("OSM test data not available")
        return jamaica_path

    def test_extract_lines(self, osm_file: Path) -> None:
        """Test extracting line features."""
        from damagescanner.osm import extract

        result = extract(
            osm_file,
            "lines",
            DICT_CIS_OSM["main_roads"]["osm_keys"],
            DICT_CIS_OSM["main_roads"]["osm_query"],
        )

        assert isinstance(result, gpd.GeoDataFrame)
        assert len(result) > 0
        assert "object_type" in result.columns
        assert "osm_id" in result.columns

    def test_extract_points(self, osm_file: Path) -> None:
        """Test extracting point features."""
        from damagescanner.osm import extract

        result = extract(
            osm_file,
            "points",
            DICT_CIS_OSM["healthcare"]["osm_keys"],
            DICT_CIS_OSM["healthcare"]["osm_query"],
        )

        assert isinstance(result, gpd.GeoDataFrame)
        assert "object_type" in result.columns

    def test_extract_multipolygons(self, osm_file: Path) -> None:
        """Test extracting multipolygon features."""
        from damagescanner.osm import extract

        result = extract(
            osm_file,
            "multipolygons",
            DICT_CIS_OSM["buildings"]["osm_keys"],
            DICT_CIS_OSM["buildings"]["osm_query"],
        )

        assert isinstance(result, gpd.GeoDataFrame)
        assert "object_type" in result.columns
