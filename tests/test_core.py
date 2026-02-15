import pytest
import pandas as pd
import geopandas as gpd
import pathlib

from damagescanner import DamageScanner

# Use the same path pattern as conftest.py
data_path = pathlib.Path(__file__).parent.parent / "data"
output_folder = pathlib.Path(__file__).parent / "output"
output_folder.mkdir(exist_ok=True)

KAMPEN = data_path / "kampen"
JAMAICA = data_path / "jamaica"
HAZARDS = data_path / "hazard"


@pytest.fixture
def raster_inputs():
    """Fixture for raster-based damage calculation inputs."""
    return {
        "hazard": KAMPEN / "hazard" / "1in100_inundation_map.tif",
        "exposure": KAMPEN / "exposure" / "landuse_map.tif",
        "curves": KAMPEN / "vulnerability" / "curves.csv",
        "maxdam": KAMPEN / "vulnerability" / "maxdam.csv",
    }


@pytest.fixture
def vector_inputs():
    """Fixture for vector-based damage calculation inputs."""
    return {
        "hazard": KAMPEN / "hazard" / "1in100_inundation_map.tif",
        "exposure": KAMPEN / "exposure" / "landuse.gpkg",
        "curves": KAMPEN / "vulnerability" / "curves_landuse.csv",
        "maxdam": KAMPEN / "vulnerability" / "maxdam_landuse.csv",
    }


@pytest.fixture
def hazard_dict():
    """Fixture for risk calculation with multiple return periods."""
    return {
        10: KAMPEN / "hazard" / "1in10_inundation_map.tif",
        50: KAMPEN / "hazard" / "1in50_inundation_map.tif",
        100: KAMPEN / "hazard" / "1in100_inundation_map.tif",
        500: KAMPEN / "hazard" / "1in500_inundation_map.tif",
        1000: KAMPEN / "hazard" / "1in1000_inundation_map.tif",
    }


@pytest.fixture
def osm_inputs():
    """Fixture for OSM-based damage calculation inputs."""
    osm_file = JAMAICA / "exposure" / "jamaica-latest.osm.pbf"
    if not osm_file.exists():
        pytest.skip("Jamaica OSM test data not available")

    return {
        "hazard": JAMAICA / "hazard" / "FD_1in1000.tif",
        "exposure": osm_file,
        "curves": JAMAICA / "vulnerability" / "curves_osm.csv",
        "maxdam": JAMAICA / "vulnerability" / "maxdam_osm.csv",
    }


class TestDamageScannerInit:
    """Tests for DamageScanner initialization."""

    def test_init_raster_from_path(self, raster_inputs):
        """Test initialization with raster exposure data."""
        ds = DamageScanner(
            hazard_data=raster_inputs["hazard"],
            feature_data=raster_inputs["exposure"],
            curves=raster_inputs["curves"],
            maxdam=raster_inputs["maxdam"],
        )
        assert ds.assessment_type == "raster"

    def test_init_vector_from_path(self, vector_inputs):
        """Test initialization with vector exposure data."""
        ds = DamageScanner(
            hazard_data=vector_inputs["hazard"],
            feature_data=vector_inputs["exposure"],
            curves=vector_inputs["curves"],
            maxdam=vector_inputs["maxdam"],
        )
        assert ds.assessment_type == "vector"

    def test_init_osm_from_path(self, osm_inputs):
        """Test initialization with OSM .pbf exposure data."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        assert ds.assessment_type == "vector"
        assert ds.osm


class TestRasterCalculation:
    """Tests for raster-based damage calculations."""

    def test_calculate_returns_correct_outputs(self, raster_inputs):
        """Test that raster calculation returns damage df, damage map, landuse, and hazard."""
        ds = DamageScanner(
            hazard_data=raster_inputs["hazard"],
            feature_data=raster_inputs["exposure"],
            curves=raster_inputs["curves"],
            maxdam=raster_inputs["maxdam"],
        )
        result = ds.calculate()

        assert isinstance(result, tuple)
        assert len(result) == 4

        damage_df, damage_map, landuse, hazard = result
        assert isinstance(damage_df, pd.DataFrame)
        assert damage_map is not None
        assert landuse is not None
        assert hazard is not None
        assert len(damage_df) > 0

    def test_damages_are_non_negative(self, raster_inputs):
        """Test that calculated damages are non-negative."""
        ds = DamageScanner(
            hazard_data=raster_inputs["hazard"],
            feature_data=raster_inputs["exposure"],
            curves=raster_inputs["curves"],
            maxdam=raster_inputs["maxdam"],
        )
        damage_df, damage_map, _, _ = ds.calculate()

        assert (damage_df["damage"] >= 0).all()
        assert (damage_map >= 0).all()


class TestVectorCalculation:
    """Tests for vector-based damage calculations."""

    def test_calculate_returns_dataframe(self, vector_inputs):
        """Test that vector calculation returns a dataframe."""
        ds = DamageScanner(
            hazard_data=vector_inputs["hazard"],
            feature_data=vector_inputs["exposure"],
            curves=vector_inputs["curves"],
            maxdam=vector_inputs["maxdam"],
        )
        result = ds.calculate(object_col="landuse")

        assert isinstance(result, pd.DataFrame)

    def test_damage_column_exists(self, vector_inputs):
        """Test that result contains damage column."""
        ds = DamageScanner(
            hazard_data=vector_inputs["hazard"],
            feature_data=vector_inputs["exposure"],
            curves=vector_inputs["curves"],
            maxdam=vector_inputs["maxdam"],
        )
        result = ds.calculate(object_col="landuse")

        if len(result) > 0:
            assert "damage" in result.columns


class TestRiskCalculation:
    """Tests for risk assessment across multiple return periods."""

    def test_risk_raster(self, raster_inputs, hazard_dict):
        """Test risk calculation with raster data."""
        ds = DamageScanner(
            hazard_data=raster_inputs["hazard"],
            feature_data=raster_inputs["exposure"],
            curves=raster_inputs["curves"],
            maxdam=raster_inputs["maxdam"],
        )
        result = ds.risk(hazard_dict)

        # Result could be None if no damages, or a DataFrame
        assert result is None or isinstance(result, pd.DataFrame)


class TestOSMExposure:
    """Tests for OSM-based exposure analysis."""

    def test_exposure_returns_geodataframe(self, osm_inputs):
        """Test that exposure() returns a GeoDataFrame for OSM data."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        result = ds.exposure(asset_type="main_roads")

        assert isinstance(result, gpd.GeoDataFrame)

    def test_exposure_has_required_columns(self, osm_inputs):
        """Test that exposure result has required columns."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        result = ds.exposure(asset_type="main_roads")

        if len(result) > 0:
            assert "object_type" in result.columns
            assert "coverage" in result.columns
            assert "values" in result.columns
            assert "geometry" in result.columns

    @pytest.mark.parametrize(
        "asset_type",
        [
            "main_roads",
            "buildings",
            "education",
        ],
    )
    def test_exposure_various_asset_types(self, osm_inputs, asset_type):
        """Test exposure analysis for various asset types."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        result = ds.exposure(asset_type=asset_type)

        assert isinstance(result, gpd.GeoDataFrame)


class TestOSMCalculation:
    """Tests for OSM-based damage calculations."""

    def test_calculate_returns_dataframe(self, osm_inputs):
        """Test that calculate() returns a DataFrame for OSM data."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        result = ds.calculate(asset_type="main_roads")

        assert isinstance(result, (pd.DataFrame, gpd.GeoDataFrame))

    def test_calculate_has_damage_column(self, osm_inputs):
        """Test that OSM calculation result has damage column."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        result = ds.calculate(asset_type="main_roads")

        if len(result) > 0:
            assert "damage" in result.columns

    def test_calculate_damages_non_negative(self, osm_inputs):
        """Test that OSM damages are non-negative."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        result = ds.calculate(asset_type="main_roads")

        if len(result) > 0:
            assert (result["damage"] >= 0).all()

    def test_calculate_total_damage_positive(self, osm_inputs):
        """Test that total damage is positive for main_roads."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        result = ds.calculate(asset_type="main_roads")

        if len(result) > 0:
            total_damage = result["damage"].sum()
            assert total_damage > 0

    @pytest.mark.parametrize(
        "asset_type",
        [
            "main_roads",
            #        "buildings",
        ],
    )
    def test_calculate_various_asset_types(self, osm_inputs, asset_type):
        """Test damage calculation for various asset types."""
        ds = DamageScanner(
            hazard_data=osm_inputs["hazard"],
            feature_data=osm_inputs["exposure"],
            curves=osm_inputs["curves"],
            maxdam=osm_inputs["maxdam"],
        )
        result = ds.calculate(asset_type=asset_type)

        assert isinstance(result, (pd.DataFrame, gpd.GeoDataFrame))
        if len(result) > 0:
            assert "damage" in result.columns


@pytest.fixture
def netcdf_inputs():
    """Fixture for NetCDF windstorm hazard inputs."""
    nc_file = KAMPEN / "hazard" / "windstorm.nc"
    if not nc_file.exists():
        pytest.skip("NetCDF windstorm test data not available")

    # Kampen OSM data
    osm_file = KAMPEN / "exposure" / "kampen.osm.pbf"
    if not osm_file.exists():
        pytest.skip("Kampen OSM test data not available")

    return {
        "hazard": nc_file,
        "exposure": osm_file,
        "curves": KAMPEN / "vulnerability" / "curves_osm.csv",
        "maxdam": KAMPEN / "vulnerability" / "maxdam_osm.csv",
    }


class TestNetCDFExposure:
    """Tests for NetCDF hazard data handling."""

    def test_init_with_netcdf_hazard(self, netcdf_inputs):
        """Test initialization with NetCDF hazard file."""
        ds = DamageScanner(
            hazard_data=netcdf_inputs["hazard"],
            feature_data=netcdf_inputs["exposure"],
            curves=netcdf_inputs["curves"],
            maxdam=netcdf_inputs["maxdam"],
        )
        assert ds.assessment_type == "vector"
        assert ds.osm

    def test_exposure_with_netcdf_buildings(self, netcdf_inputs):
        """Test exposure calculation with NetCDF hazard and OSM buildings."""
        ds = DamageScanner(
            hazard_data=netcdf_inputs["hazard"],
            feature_data=netcdf_inputs["exposure"],
            curves=netcdf_inputs["curves"],
            maxdam=netcdf_inputs["maxdam"],
        )
        result = ds.exposure(asset_type="buildings")

        assert isinstance(result, gpd.GeoDataFrame)
        print("--- NetCDF + Buildings Exposure ---")
        print(f"Features exposed: {len(result)}")

        if len(result) > 0:
            assert "coverage" in result.columns
            assert "values" in result.columns
            assert "object_type" in result.columns
            print(f"Building types: {result['object_type'].value_counts().to_dict()}")

    def test_netcdf_hazard_loads_as_xarray(self, netcdf_inputs):
        """Test that NetCDF hazard is loaded correctly as xarray."""
        import xarray as xr

        hazard = xr.open_dataset(netcdf_inputs["hazard"])

        print("--- NetCDF Structure ---")
        print(f"Variables: {list(hazard.data_vars)}")
        print(f"Dimensions: {dict(hazard.sizes)}")

        assert hazard is not None
        hazard.close()


# class TestOSMAllAssetTypes:
#     """Tests for all OSM asset types (slower, comprehensive tests)."""

#     @pytest.mark.slow
#     @pytest.mark.parametrize("asset_type", [
#         "main_roads",
#         "rail",
#         "air",
#         "telecom",
#         "water_supply",
#         "waste_solid",
#         "waste_water",
#         "education",
#         "healthcare",
#         "power",
#         "gas",
#         "oil",
#         "buildings",
#     ])
#     def test_exposure_all_asset_types(self, osm_inputs, asset_type):
#         """Test exposure analysis for all supported asset types."""
#         ds = DamageScanner(
#             hazard_data=osm_inputs["hazard"],
#             feature_data=osm_inputs["exposure"],
#             curves=osm_inputs["curves"],
#             maxdam=osm_inputs["maxdam"],
#         )
#         result = ds.exposure(asset_type=asset_type)

#         assert isinstance(result, gpd.GeoDataFrame)
#         # Some asset types might have no features in the test area
#         if len(result) > 0:
#             assert "object_type" in result.columns
#             assert "coverage" in result.columns
#             assert "values" in result.columns
