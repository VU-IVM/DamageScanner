import rasterio
import numpy as np

from .setup import tmp_folder, data_path

from damagescanner.raster import match_and_load_rasters


def test_match_and_load_rasters():
    landuse_file = data_path / "landuse" / "landuse_map.tif"
    hazard_file = data_path / "hazard" / "inundation_map.tif"

    with rasterio.open(landuse_file) as src:
        # Read the data and transform from the landuse file
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

        land_use_cropped, _, transform = match_and_load_rasters(
            extended_landuse_file, hazard_file
        )
        assert np.array_equal(land_use_cropped, landuse_data)
        assert transform == landuse_transform

        _, land_use_cropped, transform = match_and_load_rasters(
            hazard_file, extended_landuse_file
        )
        assert np.array_equal(land_use_cropped, landuse_data)
        assert transform == landuse_transform
