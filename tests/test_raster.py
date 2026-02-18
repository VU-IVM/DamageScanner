"""Tests for raster-based damage assessment functionality in DamageScanner."""

import numpy as np
import pandas as pd
import pytest
import rasterio

from damagescanner.raster import RasterScanner, match_and_load_rasters

from .helpers import data_path, tmp_folder

# Define test data paths
KAMPEN = data_path / "kampen"


@pytest.fixture
def raster_files() -> dict:
    """Fixture for raster file paths.

    Returns:
        A dictionary with paths to landuse, hazard, curves, and maxdam files for raster tests
    """
    return {
        "landuse": KAMPEN / "exposure" / "landuse_map.tif",
        "hazard": KAMPEN / "hazard" / "1in100_inundation_map.tif",
        "curves": KAMPEN / "vulnerability" / "curves.csv",
        "maxdam": KAMPEN / "vulnerability" / "maxdam.csv",
    }


class TestMatchAndLoadRasters:
    """Tests for match_and_load_rasters function."""

    def test_match_and_load_rasters_extended(self, raster_files: dict) -> None:
        """Test matching rasters when one is extended beyond the other."""
        landuse_file = raster_files["landuse"]
        hazard_file = raster_files["hazard"]

        with rasterio.open(landuse_file) as src:
            landuse_data = src.read(1)
            landuse_transform = src.transform

            # Extend the landuse data in all directions
            extended_landuse_data = np.pad(
                landuse_data, 1, mode="constant", constant_values=-1
            )

            # Calculate the new transform
            new_transform = rasterio.Affine(
                landuse_transform.a,
                landuse_transform.b,
                landuse_transform.c - landuse_transform.a,
                landuse_transform.d,
                landuse_transform.e,
                landuse_transform.f - landuse_transform.e,
            )

            # Write the extended landuse data to a new file
            extended_landuse_file = tmp_folder / "extended_landuse_map.tif"
            with rasterio.open(
                extended_landuse_file,
                "w",
                driver="GTiff",
                height=extended_landuse_data.shape[0],
                width=extended_landuse_data.shape[1],
                count=1,
                dtype=extended_landuse_data.dtype,
                crs=src.crs,
                transform=new_transform,
            ) as dst:
                dst.write(extended_landuse_data, 1)

            # Test with extended file as first argument
            land_use_cropped, _, transform = match_and_load_rasters(
                extended_landuse_file, hazard_file
            )
            assert np.array_equal(land_use_cropped, landuse_data)
            assert transform == landuse_transform

            # Test with extended file as second argument
            _, land_use_cropped, transform = match_and_load_rasters(
                hazard_file, extended_landuse_file
            )
            assert np.array_equal(land_use_cropped, landuse_data)
            assert transform == landuse_transform

    def test_match_and_load_rasters_same_extent(self, raster_files: dict) -> None:
        """Test matching rasters with identical extents."""
        landuse_file = raster_files["landuse"]
        hazard_file = raster_files["hazard"]

        data1, data2, transform = match_and_load_rasters(landuse_file, hazard_file)

        assert data1.shape == data2.shape
        assert transform is not None


class TestRasterScanner:
    """Tests for RasterScanner function."""

    def test_raster_scanner_returns_correct_types(self, raster_files: dict) -> None:
        """Test that RasterScanner returns correct output types."""
        damage_df, damagemap, landuse, hazard = RasterScanner(
            exposure_file=raster_files["landuse"],
            hazard_file=raster_files["hazard"],
            curve_path=raster_files["curves"],
            maxdam_path=raster_files["maxdam"],
        )

        assert isinstance(damage_df, pd.DataFrame)
        assert isinstance(damagemap, np.ndarray)
        assert isinstance(landuse, np.ndarray)
        assert isinstance(hazard, np.ndarray)

    def test_raster_scanner_damage_df_structure(self, raster_files: dict) -> None:
        """Test that damage DataFrame has correct structure."""
        damage_df, _, _, _ = RasterScanner(
            exposure_file=raster_files["landuse"],
            hazard_file=raster_files["hazard"],
            curve_path=raster_files["curves"],
            maxdam_path=raster_files["maxdam"],
        )

        assert "damage" in damage_df.columns
        assert len(damage_df) > 0

    def test_raster_scanner_damages_non_negative(self, raster_files: dict) -> None:
        """Test that all damage values are non-negative."""
        damage_df, damagemap, _, _ = RasterScanner(
            exposure_file=raster_files["landuse"],
            hazard_file=raster_files["hazard"],
            curve_path=raster_files["curves"],
            maxdam_path=raster_files["maxdam"],
        )

        assert (damage_df["damage"] >= 0).all()
        assert (damagemap >= 0).all()

    def test_raster_scanner_shapes_match(self, raster_files: dict) -> None:
        """Test that output arrays have matching shapes."""
        _, damagemap, landuse, hazard = RasterScanner(
            exposure_file=raster_files["landuse"],
            hazard_file=raster_files["hazard"],
            curve_path=raster_files["curves"],
            maxdam_path=raster_files["maxdam"],
        )

        assert damagemap.shape == landuse.shape
        assert damagemap.shape == hazard.shape

    def test_raster_scanner_with_dataframe_inputs(self, raster_files: dict) -> None:
        """Test RasterScanner with DataFrame inputs for curves and maxdam."""
        curves = pd.read_csv(raster_files["curves"])
        maxdam = pd.read_csv(raster_files["maxdam"])

        damage_df, damagemap, _, _ = RasterScanner(
            exposure_file=raster_files["landuse"],
            hazard_file=raster_files["hazard"],
            curve_path=curves,
            maxdam_path=maxdam,
        )

        assert isinstance(damage_df, pd.DataFrame)
        assert isinstance(damagemap, np.ndarray)

    def test_raster_scanner_with_save(self, raster_files: dict) -> None:
        """Test RasterScanner with save option."""
        damage_df, damagemap, _, _ = RasterScanner(
            exposure_file=raster_files["landuse"],
            hazard_file=raster_files["hazard"],
            curve_path=raster_files["curves"],
            maxdam_path=raster_files["maxdam"],
            save=True,
            output_path=str(tmp_folder),
            scenario_name="test_output",
        )

        # Check that output files were created
        assert (tmp_folder / "test_output_damages.csv").exists()
        assert (tmp_folder / "test_output_damagemap.tif").exists()

    def test_raster_scanner_total_damage_reasonable(self, raster_files: dict) -> None:
        """Test that total damage is a reasonable positive value."""
        damage_df, _, _, _ = RasterScanner(
            exposure_file=raster_files["landuse"],
            hazard_file=raster_files["hazard"],
            curve_path=raster_files["curves"],
            maxdam_path=raster_files["maxdam"],
        )

        total_damage = damage_df["damage"].sum()
        assert total_damage > 0
