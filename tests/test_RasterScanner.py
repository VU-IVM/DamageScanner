"""Test core objects/concepts"""

import numpy as np
import pandas as pd
import rasterio
from .setup import tmp_folder, data_path
from damagescanner.core import RasterScanner


def test_raster_scanner():
    # Define paths to example files
    landuse_file = data_path / "landuse" / "landuse_map.tif"
    hazard_file = data_path / "hazard" / "inundation_map.tif"
    curve_path = data_path / "curves" / "curves.csv"
    maxdam_path = data_path / "curves" / "maxdam.csv"

    # Call the RasterScanner function
    damage_df, damagemap, landuse_in, hazard = RasterScanner(
        landuse_file=landuse_file,
        hazard_file=hazard_file,
        curve_path=curve_path,
        maxdam_path=maxdam_path,
        lu_crs=28992,
        haz_crs=4326,
        dtype=np.int32,
        save=False,
        # Add any other required parameters or kwargs here
    )

    # Perform assertions to check the correctness of the results
    assert isinstance(damage_df, pd.DataFrame)
    assert isinstance(damagemap, np.ndarray)
    assert isinstance(landuse_in, np.ndarray)
    assert isinstance(hazard, np.ndarray)

    # Add more specific assertions based on your requirements
    landuse_file = data_path / "landuse" / "landuse_map.tif"

    with open(landuse_file, "r") as src:
        # Open landuse file
        with rasterio.open(landuse_file) as src:
            # Read the data and transform from the landuse file
            landuse_data = src.read(1)
            landuse_transform = src.transform

            # Extend the landuse data in all directions
            extended_landuse_data = np.pad(landuse_data, 1, mode="constant")

            # Calculate the new transform
            new_transform = rasterio.Affine(
                landuse_transform.a,
                landuse_transform.b,
                landuse_transform.c - landuse_transform.a,
                landuse_transform.d,
                landuse_transform.e,
                landuse_transform.f - landuse_transform.e,
            )

            # Write the extended landuse data with the new transform to a new file
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

    hazard_file = data_path / "hazard" / "inundation_map.tif"
    curve_path = data_path / "curves" / "curves.csv"
    maxdam_path = data_path / "curves" / "maxdam.csv"

    # Call the RasterScanner function
    damage_df, damagemap, landuse_in, hazard = RasterScanner(
        landuse_file=extended_landuse_file,
        hazard_file=hazard_file,
        curve_path=curve_path,
        maxdam_path=maxdam_path,
        lu_crs=28992,
        haz_crs=4326,
        dtype=np.int32,
        save=False,
        # Add any other required parameters or kwargs here
    )

    # Perform assertions to check the correctness of the results
    assert isinstance(damage_df, pd.DataFrame)
    assert isinstance(damagemap, np.ndarray)
    assert isinstance(landuse_in, np.ndarray)
    assert isinstance(hazard, np.ndarray)

    # Add more specific assertions based on your requirements
