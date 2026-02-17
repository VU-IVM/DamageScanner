import pytest
import pandas as pd
import geopandas as gpd
import shapely

from .helpers import data_path

from damagescanner.vector import (
    VectorScanner,
    VectorExposure,
    _create_grid,
    _remove_duplicates,
    _get_cell_area_m2,
    _overlay_vector_vector,
)


# Define test data paths
KAMPEN = data_path / "kampen"
JAMAICA = data_path / "jamaica"


@pytest.fixture
def vector_files():
    """Fixture for vector file paths."""
    return {
        "hazard": KAMPEN / "hazard" / "1in100_inundation_map.tif",
        "exposure": KAMPEN / "exposure" / "landuse.gpkg",
        "curves": KAMPEN / "vulnerability" / "curves_landuse.csv",
        "maxdam": KAMPEN / "vulnerability" / "maxdam_landuse.csv",
    }


@pytest.fixture
def sample_geodataframe():
    """Create a simple GeoDataFrame for testing."""
    return gpd.GeoDataFrame(
        {
            "id": [1, 2, 3],
            "landuse": [111, 112, 120],
            "geometry": [
                shapely.Point(5.9, 52.5),
                shapely.Point(5.91, 52.51),
                shapely.Point(5.92, 52.52),
            ],
        },
        crs="EPSG:4326",
    )


class TestCreateGrid:
    """Tests for _create_grid function."""

    def test_create_grid_returns_polygons(self):
        """Test that _create_grid returns shapely polygons."""
        bbox = shapely.box(0, 0, 10, 10)
        grid = _create_grid(bbox, height=2)

        assert len(grid) > 0
        assert all(shapely.get_type_id(g) == 3 for g in grid)  # 3 = Polygon

    def test_create_grid_correct_number_of_cells(self):
        """Test that grid has expected number of cells."""
        bbox = shapely.box(0, 0, 10, 10)
        grid = _create_grid(bbox, height=2)

        # 10/2 = 5 cells per side, so 5*5 = 25 cells
        assert len(grid) == 25

    def test_create_grid_covers_bbox(self):
        """Test that grid cells cover the bounding box."""
        bbox = shapely.box(0, 0, 10, 10)
        grid = _create_grid(bbox, height=2)

        # Union of all grid cells should cover the bbox
        grid_union = shapely.union_all(grid)
        assert shapely.covers(grid_union, bbox)


class TestRemoveDuplicates:
    """Tests for _remove_duplicates function."""

    def test_remove_duplicates_no_duplicates(self):
        """Test with no duplicate indices."""
        df = pd.DataFrame(
            {
                "coverage": [[0.5, 0.3], [0.2, 0.4]],
                "values": [[1.0, 2.0], [3.0, 4.0]],
            },
            index=[0, 1],
        )

        result = _remove_duplicates(df)

        assert len(result) == 2
        assert list(result.index) == [0, 1]

    def test_remove_duplicates_with_duplicates(self):
        """Test with duplicate indices - should concatenate lists."""
        df = pd.DataFrame(
            {
                "coverage": [[0.5], [0.3], [0.2]],
                "values": [[1.0], [2.0], [3.0]],
            },
            index=[0, 0, 1],
        )

        result = _remove_duplicates(df)

        assert len(result) == 2
        # Index 0 should have concatenated lists
        assert result.loc[0, "coverage"] == [0.5, 0.3]
        assert result.loc[0, "values"] == [1.0, 2.0]


class TestGetCellAreaM2:
    """Tests for _get_cell_area_m2 function."""

    def test_get_cell_area_returns_positive(self, sample_geodataframe):
        """Test that cell area is positive."""
        resolution = 0.001  # degrees
        area = _get_cell_area_m2(sample_geodataframe, resolution)

        assert area > 0

    def test_get_cell_area_scales_with_resolution(self, sample_geodataframe):
        """Test that area scales with resolution squared."""
        area_small = _get_cell_area_m2(sample_geodataframe, 0.001)
        area_large = _get_cell_area_m2(sample_geodataframe, 0.002)

        # Area should scale approximately with resolution squared
        ratio = area_large / area_small
        assert 3.5 < ratio < 4.5  # Should be ~4 (2^2)


class TestVectorExposure:
    """Tests for VectorExposure function."""

    def test_vector_exposure_returns_tuple(self, vector_files):
        """Test that VectorExposure returns correct tuple structure."""
        result = VectorExposure(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
        )

        assert isinstance(result, tuple)
        assert len(result) == 4
        features, object_col, hazard_crs, cell_area_m2 = result

    def test_vector_exposure_returns_geodataframe(self, vector_files):
        """Test that features are returned as GeoDataFrame."""
        features, _, _, _ = VectorExposure(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
        )

        assert isinstance(features, gpd.GeoDataFrame)

    def test_vector_exposure_adds_coverage_and_values(self, vector_files):
        """Test that coverage and values columns are added."""
        features, _, _, _ = VectorExposure(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
        )

        if len(features) > 0:
            assert "coverage" in features.columns
            assert "values" in features.columns

    def test_vector_exposure_cell_area_positive(self, vector_files):
        """Test that cell area is positive."""
        _, _, _, cell_area_m2 = VectorExposure(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
        )

        assert cell_area_m2 > 0


class TestVectorScanner:
    """Tests for VectorScanner function."""

    def test_vector_scanner_returns_geodataframe(self, vector_files):
        """Test that VectorScanner returns a GeoDataFrame."""
        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        assert isinstance(result, (gpd.GeoDataFrame, pd.DataFrame))

    def test_vector_scanner_has_damage_column(self, vector_files):
        """Test that result contains damage column."""
        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        if len(result) > 0:
            assert "damage" in result.columns

    def test_vector_scanner_damages_non_negative(self, vector_files):
        """Test that all damage values are non-negative."""
        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        if len(result) > 0:
            assert (result["damage"] >= 0).all()

    def test_vector_scanner_preserves_geometry(self, vector_files):
        """Test that geometry column is preserved."""
        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        if len(result) > 0:
            assert "geometry" in result.columns

    def test_vector_scanner_with_dataframe_inputs(self, vector_files):
        """Test VectorScanner with DataFrame inputs for curves and maxdam."""
        curves = pd.read_csv(vector_files["curves"], index_col=0)
        maxdam = pd.read_csv(vector_files["maxdam"])

        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=curves,
            maxdam_path=maxdam,
            object_col="landuse",
        )

        assert isinstance(result, (gpd.GeoDataFrame, pd.DataFrame))

    def test_vector_scanner_with_geodataframe_input(self, vector_files):
        """Test VectorScanner with GeoDataFrame input for features."""
        features = gpd.read_file(vector_files["exposure"])

        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=features,
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        assert isinstance(result, (gpd.GeoDataFrame, pd.DataFrame))

    def test_vector_scanner_gridded_vs_non_gridded(self, vector_files):
        """Test that gridded and non-gridded produce similar results."""
        result_gridded = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
            gridded=True,
        )

        result_non_gridded = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
            gridded=False,
        )

        # Both should return similar total damage (within tolerance due to edge effects)
        if len(result_gridded) > 0 and len(result_non_gridded) > 0:
            total_gridded = result_gridded["damage"].sum()
            total_non_gridded = result_non_gridded["damage"].sum()

            assert (
                abs(total_gridded - total_non_gridded)
                / max(total_gridded, total_non_gridded)
                < 0.05
            )

    def test_vector_scanner_total_damage_positive(self, vector_files):
        """Test that total damage is positive."""
        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        if len(result) > 0:
            assert result["damage"].sum() > 0

    def test_vector_scanner_empty_features(self, vector_files):
        """Test VectorScanner handles empty features gracefully."""
        # Create empty GeoDataFrame
        empty_features = gpd.GeoDataFrame(
            {"landuse": [], "geometry": []},
            crs="EPSG:28992",
        )

        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=empty_features,
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        assert len(result) == 0


@pytest.fixture
def vector_hazard_files():
    """Fixture for vector-vector overlay test files."""
    vector_hazard = KAMPEN / "hazard" / "1in100_inundation_map_vector.parquet"
    if not vector_hazard.exists():
        pytest.skip(
            "Vector hazard test data not available - run examples/create_vector_hazard.py first"
        )

    return {
        "hazard": vector_hazard,
        "exposure": KAMPEN / "exposure" / "landuse.gpkg",
        "curves": KAMPEN / "vulnerability" / "curves_landuse.csv",
        "maxdam": KAMPEN / "vulnerability" / "maxdam_landuse.csv",
    }


class TestVectorVectorOverlay:
    """Tests for vector-vector overlay functionality."""

    def test_overlay_vector_vector_returns_geodataframe(self, vector_hazard_files):
        """Test that vector-vector overlay returns a GeoDataFrame."""
        hazard = gpd.read_parquet(vector_hazard_files["hazard"])
        features = gpd.read_file(vector_hazard_files["exposure"])

        result = _overlay_vector_vector(
            hazard=hazard,
            features=features,
            hazard_value_col="band_data",
            disable_progress=False,
        )

        assert isinstance(result, gpd.GeoDataFrame)
        print(f"\nFeatures with exposure: {len(result)}")

    def test_overlay_vector_vector_has_required_columns(self, vector_hazard_files):
        """Test that result has coverage and values columns."""
        hazard = gpd.read_parquet(vector_hazard_files["hazard"])
        features = gpd.read_file(vector_hazard_files["exposure"])

        result = _overlay_vector_vector(
            hazard=hazard,
            features=features,
            hazard_value_col="band_data",
            disable_progress=False,
        )

        assert "coverage" in result.columns
        assert "values" in result.columns

    def test_overlay_vector_vector_values_are_lists(self, vector_hazard_files):
        """Test that values and coverage are lists."""
        hazard = gpd.read_parquet(vector_hazard_files["hazard"])
        features = gpd.read_file(vector_hazard_files["exposure"])

        result = _overlay_vector_vector(
            hazard=hazard,
            features=features,
            hazard_value_col="band_data",
            disable_progress=False,
        )

        if len(result) > 0:
            assert all(isinstance(v, list) for v in result["values"])
            assert all(isinstance(c, list) for c in result["coverage"])

    def test_overlay_vector_vector_coverage_positive(self, vector_hazard_files):
        """Test that coverage values are positive."""
        hazard = gpd.read_parquet(vector_hazard_files["hazard"])
        features = gpd.read_file(vector_hazard_files["exposure"])

        result = _overlay_vector_vector(
            hazard=hazard,
            features=features,
            hazard_value_col="band_data",
            disable_progress=False,
        )

        if len(result) > 0:
            for coverage_list in result["coverage"]:
                assert all(c > 0 for c in coverage_list)

    def test_vector_exposure_with_vector_hazard(self, vector_hazard_files):
        """Test VectorExposure with vector hazard input."""
        features, object_col, hazard_crs, cell_area_m2 = VectorExposure(
            hazard_file=vector_hazard_files["hazard"],
            feature_file=vector_hazard_files["exposure"],
            object_col="landuse",
            hazard_value_col="band_data",
        )

        assert isinstance(features, gpd.GeoDataFrame)
        assert cell_area_m2 == 1  # Vector hazards use coverage in meters directly
        print(f"\nExposed features: {len(features)}")

    def test_vector_scanner_with_vector_hazard(self, vector_hazard_files):
        """Test VectorScanner with vector hazard input."""
        result = VectorScanner(
            hazard_file=vector_hazard_files["hazard"],
            feature_file=vector_hazard_files["exposure"],
            curve_path=vector_hazard_files["curves"],
            maxdam_path=vector_hazard_files["maxdam"],
            object_col="landuse",
            hazard_value_col="band_data",
        )

        assert isinstance(result, gpd.GeoDataFrame)
        print(f"\nFeatures with damage: {len(result)}")

        if len(result) > 0:
            assert "damage" in result.columns
            print(f"Total damage: {result['damage'].sum():,.2f}")

    def test_vector_scanner_with_vector_hazard_damages_non_negative(
        self, vector_hazard_files
    ):
        """Test that damages are non-negative."""
        result = VectorScanner(
            hazard_file=vector_hazard_files["hazard"],
            feature_file=vector_hazard_files["exposure"],
            curve_path=vector_hazard_files["curves"],
            maxdam_path=vector_hazard_files["maxdam"],
            object_col="landuse",
            hazard_value_col="band_data",
        )

        if len(result) > 0:
            assert (result["damage"] >= 0).all()


@pytest.fixture
def osm_files():
    """Fixture for Jamaica OSM test files."""
    osm_path = JAMAICA / "exposure" / "jamaica-latest.osm.pbf"
    hazard_path = JAMAICA / "hazard" / "FD_1in1000.tif"

    if not osm_path.exists() or not hazard_path.exists():
        pytest.skip("Jamaica OSM/hazard test data not available")

    return {
        "hazard": hazard_path,
        "exposure": osm_path,
        "curves": JAMAICA / "vulnerability" / "curves_osm.csv",
        "maxdam": JAMAICA / "vulnerability" / "maxdam_osm.csv",
    }


class TestExtractStrategy:
    """Tests for exactextract strategy comparison."""

    def test_gridded_feature_vs_raster_sequential_jamaica_roads(self, osm_files):
        """Test that feature-sequential and raster-sequential produce same results for gridded."""
        # Run with feature-sequential
        result_feature = VectorScanner(
            hazard_file=osm_files["hazard"],
            feature_file=osm_files["exposure"],
            curve_path=osm_files["curves"],
            maxdam_path=osm_files["maxdam"],
            asset_type="roads",
            gridded=True,
            extract_strategy="feature-sequential",
        )

        # Run with raster-sequential
        result_raster = VectorScanner(
            hazard_file=osm_files["hazard"],
            feature_file=osm_files["exposure"],
            curve_path=osm_files["curves"],
            maxdam_path=osm_files["maxdam"],
            asset_type="roads",
            gridded=True,
            extract_strategy="raster-sequential",
        )

        assert len(result_feature) > 0, "Feature-sequential returned no results"
        assert len(result_raster) > 0, "Raster-sequential returned no results"

        total_feature = result_feature["damage"].sum()
        total_raster = result_raster["damage"].sum()

        print("\n--- Extract Strategy Comparison (Gridded) ---")
        print(f"Feature-sequential total damage: {total_feature:,.2f}")
        print(f"Raster-sequential total damage: {total_raster:,.2f}")
        print(f"Difference: {abs(total_feature - total_raster):,.2f}")

        assert total_feature == pytest.approx(total_raster, rel=0.01), (
            f"Feature-sequential ({total_feature:,.2f}) and raster-sequential ({total_raster:,.2f}) "
            "should produce equal results"
        )
