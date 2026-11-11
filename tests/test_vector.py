"""Tests for vector-based damage assessment functionality in DamageScanner."""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import pytest
import rasterio
import shapely

from damagescanner.vector import (
    VectorExposure,
    VectorScanner,
    _create_grid,
    _get_cell_area_m2,
    _overlay_raster_vector,
    _overlay_vector_vector,
    _remove_duplicates,
)

from .helpers import data_path

# Define test data paths
KAMPEN = data_path / "kampen"
JAMAICA = data_path / "jamaica"


@pytest.fixture
def vector_files() -> dict:
    """Fixture for vector file paths.

    Returns:
        A dictionary with paths to hazard, exposure, curves, and maxdam files for vector tests
    """
    return {
        "hazard": KAMPEN / "hazard" / "1in100_inundation_map.tif",
        "exposure": KAMPEN / "exposure" / "landuse.gpkg",
        "curves": KAMPEN / "vulnerability" / "curves_landuse.csv",
        "maxdam": KAMPEN / "vulnerability" / "maxdam_landuse.csv",
    }


@pytest.fixture
def sample_geodataframe() -> gpd.GeoDataFrame:
    """Create a simple GeoDataFrame for testing.

    Returns:
        A GeoDataFrame with sample points and landuse values.
    """
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

    def test_create_grid_returns_polygons(self) -> None:
        """Test that _create_grid returns shapely polygons."""
        bbox = shapely.box(0, 0, 10, 10)
        grid = _create_grid(bbox, height=2)

        assert len(grid) > 0
        assert all(shapely.get_type_id(g) == 3 for g in grid)  # 3 = Polygon

    def test_create_grid_correct_number_of_cells(self) -> None:
        """Test that grid has expected number of cells."""
        bbox = shapely.box(0, 0, 10, 10)
        grid = _create_grid(bbox, height=2)

        # 10/2 = 5 cells per side, so 5*5 = 25 cells
        assert len(grid) == 25

    def test_create_grid_covers_bbox(self) -> None:
        """Test that grid cells cover the bounding box."""
        bbox = shapely.box(0, 0, 10, 10)
        grid = _create_grid(bbox, height=2)

        # Union of all grid cells should cover the bbox
        grid_union = shapely.union_all(grid)
        assert shapely.covers(grid_union, bbox)


class TestRemoveDuplicates:
    """Tests for _remove_duplicates function."""

    def test_remove_duplicates_no_duplicates(self) -> None:
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

    def test_remove_duplicates_with_duplicates(self) -> None:
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

    def test_get_cell_area_returns_positive(
        self, sample_geodataframe: gpd.GeoDataFrame
    ) -> None:
        """Test that cell area is positive."""
        resolution = 0.001  # degrees
        area = _get_cell_area_m2(
            sample_geodataframe, pyproj.CRS.from_epsg(4326), resolution
        )

        assert area > 0

    def test_get_cell_area_scales_with_resolution(
        self, sample_geodataframe: gpd.GeoDataFrame
    ) -> None:
        """Test that area scales with resolution squared."""
        area_small = _get_cell_area_m2(
            sample_geodataframe, pyproj.CRS.from_epsg(4326), 0.001
        )
        area_large = _get_cell_area_m2(
            sample_geodataframe, pyproj.CRS.from_epsg(4326), 0.002
        )

        # Area should scale approximately with resolution squared
        ratio = area_large / area_small
        assert 3.5 < ratio < 4.5  # Should be ~4 (2^2)


class TestVectorExposure:
    """Tests for VectorExposure function."""

    def test_vector_exposure_returns_tuple(self, vector_files: dict) -> None:
        """Test that VectorExposure returns correct tuple structure."""
        result = VectorExposure(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
        )

        assert isinstance(result, tuple)
        assert len(result) == 4
        features, object_col, hazard_crs, cell_area_m2 = result

    def test_vector_exposure_returns_geodataframe(self, vector_files: dict) -> None:
        """Test that features are returned as GeoDataFrame."""
        features, _, _, _ = VectorExposure(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
        )

        assert isinstance(features, gpd.GeoDataFrame)

    def test_vector_exposure_adds_coverage_and_values(self, vector_files: dict) -> None:
        """Test that coverage and values columns are added."""
        features, _, _, _ = VectorExposure(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
        )

        if len(features) > 0:
            assert "coverage" in features.columns
            assert "values" in features.columns

    def test_vector_exposure_cell_area_positive(self, vector_files: dict) -> None:
        """Test that cell area is positive."""
        _, _, _, cell_area_m2 = VectorExposure(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
        )
        assert cell_area_m2 is not None
        assert cell_area_m2 > 0


class TestVectorScanner:
    """Tests for VectorScanner function."""

    def test_vector_scanner_returns_geodataframe(self, vector_files: dict) -> None:
        """Test that VectorScanner returns a GeoDataFrame."""
        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        assert isinstance(result, (gpd.GeoDataFrame, pd.DataFrame))

    def test_vector_scanner_has_damage_column(self, vector_files: dict) -> None:
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

    def test_vector_scanner_damages_non_negative(self, vector_files: dict) -> None:
        """Test that all damage values are non-negative."""
        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
        )

        if len(result) > 0:
            assert ((result["damage"] >= 0) | result["damage"].isnull()).all()

        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
            return_full=False,
        )

        if len(result) > 0:
            assert (result["damage"] >= 0).all()

    def test_vector_scanner_preserves_geometry(self, vector_files: dict) -> None:
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

    def test_vector_scanner_with_dataframe_inputs(self, vector_files: dict) -> None:
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

    def test_vector_scanner_with_geodataframe_input(self, vector_files: dict) -> None:
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

    def test_vector_scanner_gridded_vs_non_gridded(self, vector_files: dict) -> None:
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

    def test_vector_scanner_total_damage_positive(self, vector_files: dict) -> None:
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

    def test_vector_scanner_empty_features(self, vector_files: dict) -> None:
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

    def test_vector_scanner_return_full_on(self, vector_files: dict) -> None:
        """Test that return_full=True returns all features."""
        # Get total count of features first
        features = gpd.read_file(vector_files["exposure"])

        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
            return_full=True,
        )

        assert len(result) == len(features)

    def test_vector_scanner_return_full_off(self, vector_files: dict) -> None:
        """Test that return_full=False returns only features with hazard intersection."""
        # Get total count of features first
        features_all = gpd.read_file(vector_files["exposure"])

        result = VectorScanner(
            hazard_file=vector_files["hazard"],
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
            return_full=False,
        )

        # Features with no hazard intersection should be removed
        assert len(result) < len(features_all)
        # All returned features should have non-empty values
        assert result["values"].apply(lambda x: len(x) > 0).all()

    def test_vector_scanner_nodata_hazard(
        self, vector_files: dict, tmp_path: Path
    ) -> None:
        """Test VectorScanner with hazard that has no overlapping data."""
        # Create a hazard raster that is far away from the exposure
        with rasterio.open(vector_files["hazard"]) as src:
            data = src.read()
            profile = src.profile.copy()

            # Move the hazard far away
            profile.update(
                {"transform": rasterio.Affine(1, 0, 1000000, 0, -1, 1000000)}
            )

            nodata_hazard = tmp_path / "nodata_hazard.tif"
            with rasterio.open(nodata_hazard, "w", **profile) as dst:
                dst.write(np.zeros_like(data))

        # Test with return_full=True
        result_full = VectorScanner(
            hazard_file=nodata_hazard,
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
            return_full=True,
        )
        assert len(result_full) > 0
        assert (result_full["damage"] == 0).all()

        # Test with return_full=False
        result_none = VectorScanner(
            hazard_file=nodata_hazard,
            feature_file=vector_files["exposure"],
            curve_path=vector_files["curves"],
            maxdam_path=vector_files["maxdam"],
            object_col="landuse",
            return_full=False,
        )
        assert len(result_none) == 0


@pytest.fixture
def vector_hazard_files() -> dict:
    """Fixture for vector-vector overlay test files.

    Returns:
        A dictionary with paths to hazard, exposure, curves, and maxdam files for vector hazard tests.
    """
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

    def test_overlay_vector_vector_returns_geodataframe(
        self, vector_hazard_files: dict
    ) -> None:
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

    def test_overlay_vector_vector_has_required_columns(
        self, vector_hazard_files: dict
    ) -> None:
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

    def test_overlay_vector_vector_values_are_lists(
        self, vector_hazard_files: dict
    ) -> None:
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

    def test_overlay_vector_vector_coverage_positive(
        self, vector_hazard_files: dict
    ) -> None:
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

    def test_vector_exposure_with_vector_hazard(
        self, vector_hazard_files: dict
    ) -> None:
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

    def test_vector_scanner_with_vector_hazard(self, vector_hazard_files: dict) -> None:
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
        self, vector_hazard_files: dict
    ) -> None:
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
def osm_files() -> dict:
    """Fixture for Jamaica OSM test files.

    Returns:
        A dictionary with paths to hazard, exposure, curves, and maxdam files for Jamaica.
    """
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

    def test_gridded_feature_vs_raster_sequential_jamaica_roads(
        self, osm_files: dict
    ) -> None:
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

class TestExposedAreaInvariant:
    """Exposed polygon area must equal the polygon's own area (issues #82, #83)."""
 
    @staticmethod
    def _raster(x0, y0, res, n, crs):
        import rioxarray  # noqa: F401
        import xarray as xr
 
        xs = x0 + res / 2 + np.arange(n) * res
        ys = y0 + n * res - res / 2 - np.arange(n) * res
        return xr.DataArray(
            np.ones((n, n), dtype="float32"), coords={"y": ys, "x": xs}, dims=("y", "x")
        ).rio.write_crs(crs)
 
    @pytest.mark.parametrize("gridded", [True, False])
    def test_projected_tiles_not_aligned_with_cells(self, gridded):
        """90 m cells do not divide the 50 km tiles; boundary cells must count once."""
        res, x0, y0 = 90.0, 3_000_000.0, 2_000_000.0
        hazard = self._raster(x0, y0, res, 1300, "EPSG:3035")
        polys = [
            shapely.box(x0 + 45_000, y0 + 10_000, x0 + 55_000, y0 + 20_000),
            shapely.box(x0 + 49_900, y0 + 1_000, x0 + 50_100, y0 + 110_000),
        ]
        feats = gpd.GeoDataFrame(geometry=polys, crs="EPSG:3035")
        out = _overlay_raster_vector(
            hazard, feats, pyproj.CRS(3035), gridded=gridded, disable_progress=True
        )
        for cov, poly in zip(out["coverage"], polys):
            assert np.sum(cov) * res * res / poly.area == pytest.approx(1, abs=0.002)
 
    @pytest.mark.parametrize("lat0", [-60, -40, 0, 40, 60])
    def test_geographic_polygon_area_in_m2(self, lat0):
        """Polygon coverage on a geographic raster is converted per cell to m².
 
        Both hemispheres: the area depends on distance from the equator only.
        """
        from pyproj import Geod
 
        res = 1 / 120
        hazard = self._raster(10.0, lat0, res, 400, "EPSG:4326")
        poly = shapely.box(10.5, lat0 + 0.5, 12.6, lat0 + 2.9)
        feats = gpd.GeoDataFrame(geometry=[poly], crs="EPSG:4326")
        out = _overlay_raster_vector(
            hazard, feats, pyproj.CRS(4326), gridded=True, disable_progress=True
        )
        true_area = abs(Geod(ellps="WGS84").geometry_area_perimeter(poly)[0])
        assert np.sum(out["coverage"].iloc[0]) / true_area == pytest.approx(1, abs=0.002)

def _multi_curve_case(hazard_value: float = 2.0):
    """A 200 x 200 m power plant and two power tower points on a uniform hazard.
 
    Values are chosen so every damage can be checked by hand.
    """
    import rioxarray  # noqa: F401
    import xarray as xr
 
    res, n, x0, y0 = 100.0, 50, 500_000.0, 5_600_000.0
    xs = x0 + res / 2 + np.arange(n) * res
    ys = (y0 + res / 2 + np.arange(n) * res)[::-1]
    hazard = xr.DataArray(
        np.full((n, n), hazard_value, dtype="float32"),
        coords={"y": ys, "x": xs},
        dims=("y", "x"),
    ).rio.write_crs("EPSG:32632")
    features = gpd.GeoDataFrame(
        {"osm_id": [1, 2, 3], "object_type": ["plant", "tower", "tower"]},
        geometry=[
            shapely.box(x0 + 1000, y0 + 1000, x0 + 1200, y0 + 1200),
            shapely.Point(x0 + 2050, y0 + 2050),
            shapely.Point(x0 + 3050, y0 + 3050),
        ],
        crs="EPSG:32632",
    )
    depth = [0, 1, 2, 3]
    flood = pd.DataFrame({"plant": [0, 0.2, 0.4, 0.6], "tower": [0, 0.1, 0.2, 0.3]}, index=depth)
    collapse = pd.DataFrame({"tower": [0, 0.5, 1.0, 1.0]}, index=depth)  # towers only
    maxdam = {"plant": 1000, "tower": 100}
    return hazard, features, flood, collapse, maxdam
 
 
class TestMultiCurves:
    """A curve set only needs curves for the object types it applies to (#58)."""
 
    def test_partial_curve_sets(self) -> None:
        hazard, features, flood, collapse, maxdam = _multi_curve_case()
        out = VectorScanner(
            hazard, features, flood, maxdam,
            multi_curves={"flood": flood, "collapse": collapse},
            disable_progress=True,
        )
        # plant: 0.4 damage fraction x 40,000 m2 x 1000 per m2
        assert out.loc[0, "flood"] == pytest.approx(16_000_000)
        # towers: 0.2 x 100 for flood, 1.0 x 100 for collapse
        assert list(out.loc[[1, 2], "flood"]) == pytest.approx([20, 20])
        assert list(out.loc[[1, 2], "collapse"]) == pytest.approx([100, 100])
        # collapse does not apply to the plant: NaN, not 0
        assert np.isnan(out.loc[0, "collapse"])
 
    def test_object_type_in_no_curve_set_warns(self) -> None:
        hazard, features, flood, collapse, maxdam = _multi_curve_case()
        with pytest.warns(UserWarning, match="plant"):
            out = VectorScanner(
                hazard, features, flood, maxdam,
                multi_curves={"collapse": collapse},
                disable_progress=True,
            )
        assert np.isnan(out.loc[0, "collapse"])
 
    def test_single_curve_missing_object_type_raises(self) -> None:
        hazard, features, _, collapse, maxdam = _multi_curve_case()
        with pytest.raises(KeyError, match="plant"):
            VectorScanner(hazard, features, collapse, maxdam, disable_progress=True)