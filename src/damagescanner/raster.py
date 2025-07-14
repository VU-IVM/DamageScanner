"""DamageScanner - a directe damage assessment toolkit

Raster specific functions
"""

import rasterio
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
from affine import Affine
import warnings
from pathlib import PurePath
import shapely
import rasterio.transform
from rasterio.windows import Window

from damagescanner.utils import _check_output_path, _check_scenario_name


def match_and_load_rasters(raster_in1, raster_in2):
    """
    Match and clip two raster files to their common spatial extent and resolution.

    Code adapted from http://sciience.tumblr.com/post/101722591382/finding-the-georeferenced-intersection-between-two


    Args:
        raster_in1 (str or Path): Path to the first raster file.
        raster_in2 (str or Path): Path to the second raster file.

    Returns:
        tuple:
            - np.ndarray: Clipped raster array from the first file.
            - np.ndarray: Clipped raster array from the second file.
            - Affine: Affine transform of the intersecting region.
    """
    with rasterio.open(raster_in1) as src1, rasterio.open(raster_in2) as src2:
        if src1.crs != src2.crs:
            raise ValueError("Different CRS: CRS must be the same.")
        if src1.res != src2.res:
            raise ValueError("Different resolution: Cell sizes must be the same.")

        top_delta = round((src2.bounds.top - src1.bounds.top) / src1.transform.e)
        bottom_delta = round(
            (src2.bounds.bottom - src1.bounds.bottom) / src1.transform.e
        )
        left_delta = round((src2.bounds.left - src1.bounds.left) / src1.transform.a)
        right_delta = round((src2.bounds.right - src1.bounds.right) / src1.transform.a)

        data1 = src1.read(
            1,
            window=Window(
                col_off=left_delta,
                row_off=top_delta,
                width=src1.width - left_delta + right_delta,
                height=src1.height - top_delta + bottom_delta,
            ),
        )
        data2 = src2.read(
            1,
            window=Window(
                col_off=abs(min(left_delta, 0)),
                row_off=abs(min(top_delta, 0)),
                width=max(src1.width, src2.width) - abs(left_delta) - abs(right_delta),
                height=max(src1.height, src2.height)
                - abs(top_delta)
                - abs(bottom_delta),
            ),
        )
        transform = rasterio.Affine(
            src1.transform.a,
            src1.transform.b,
            src1.transform.c + src1.transform.a * max(left_delta, 0),
            src1.transform.d,
            src1.transform.e,
            src1.transform.f + src1.transform.e * max(top_delta, 0),
        )

    return data1, data2, transform


def _match_raster_to_vector(hazard, landuse, lu_crs, haz_crs, resolution, hazard_col):
    """
    Align hazard and land-use rasters by CRS, extent, and resolution.

    Args:
        hazard (xr.Dataset): Hazard dataset.
        landuse (xr.Dataset): Land-use dataset.
        lu_crs (int): EPSG code of the land-use CRS.
        haz_crs (int): EPSG code of the hazard CRS.
        resolution (float): Target resolution in projected units.
        hazard_col (str): Column name containing hazard intensity.

    Returns:
        tuple:
            - np.ndarray: Reprojected hazard raster as array.
            - np.ndarray: Reprojected land-use raster as array.
    """
    # Set the crs of the hazard variable to haz_crs
    hazard.rio.write_crs(haz_crs, inplace=True)

    # Rename the latitude and longitude variables to 'y' and 'x' respectively
    hazard = hazard.rename({"Latitude": "y", "Longitude": "x"})

    # Set the x and y dimensions in the hazard variable to 'x' and 'y' respectively
    hazard.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)

    # Set the crs of the landuse variable to lu_crs
    landuse.rio.write_crs(lu_crs, inplace=True)

    # Reproject the landuse variable from EPSG:4326 to EPSG:3857
    landuse = landuse.rio.reproject("EPSG:3857", resolution=resolution)

    # Get the minimum longitude and latitude values in the landuse variable
    min_lon = landuse.x.min().to_dict()["data"]
    min_lat = landuse.y.min().to_dict()["data"]

    # Get the maximum longitude and latitude values in the landuse variable
    max_lon = landuse.x.max().to_dict()["data"]
    max_lat = landuse.y.max().to_dict()["data"]

    # Create a bounding box using the minimum and maximum latitude and longitude values
    area = gpd.GeoDataFrame(
        [shapely.box(min_lon, min_lat, max_lon, max_lat)], columns=["geometry"]
    )

    # Set the crs of the bounding box to EPSG:3857
    area.crs = "epsg:3857"

    # Convert the crs of the bounding box to EPSG:4326
    area = area.to_crs("epsg:4326")

    # Clip the hazard variable to the extent of the bounding box
    hazard = hazard.rio.clip(area.geometry.values, area.crs)

    # Reproject the hazard variable to EPSG:3857 with the desired resolution
    hazard = hazard.rio.reproject("EPSG:3857", resolution=resolution)

    # Clip the hazard variable again to the extent of the bounding box
    hazard = hazard.rio.clip(area.geometry.values, area.crs)

    # If the hazard variable has fewer columns and rows than the landuse variable, reproject
    # the landuse variable to match the hazard variable
    if (len(hazard.x) < len(landuse.x)) & (len(hazard.y) < len(landuse.y)):
        landuse = landuse.rio.reproject_match(hazard)

    # If the hazard variable has more columns and rows than the landuse variable,
    # reproject the hazard variable to match the landuse variable

    elif (len(hazard.x) > len(landuse.x)) & (len(hazard.y) > len(landuse.y)):
        hazard = hazard.rio.reproject_match(landuse)

    # Convert the hazard and landuse variable to a numpy array
    landuse = landuse["band_data"].to_numpy()[0, :, :]
    hazard = hazard[hazard_col].to_numpy()[0, :, :]

    return hazard, landuse


def RasterScanner(
    exposure_file,
    hazard_file,
    curve_path,
    maxdam_path,
    lu_crs=28992,
    haz_crs=4326,
    hazard_col="FX",
    dtype=np.int32,
    save=False,
    **kwargs,
):
    """
    Run a raster-based direct damage assessment using hazard and exposure layers.

    Args:
        exposure_file (Path | np.ndarray): Path to land-use GeoTIFF or numpy array.
        hazard_file (Path | np.ndarray | xr.Dataset): Path to hazard raster or dataset.
        curve_path (Path | pd.DataFrame | np.ndarray): Vulnerability curve(s).
        maxdam_path (Path | pd.DataFrame | np.ndarray): Maximum damage values.
        lu_crs (int): CRS of the land-use file (default EPSG:28992).
        haz_crs (int): CRS of the hazard file (default EPSG:4326).
        hazard_col (str): Column containing hazard intensity (default "FX").
        dtype (type): Output dtype for damage raster.
        save (bool): If True, saves damage results to file.

    Keyword Args:
        nan_value (float): Replace this value in the hazard raster with 0.
        cellsize (float): Cell size (m²) if exposure and hazard are arrays.
        resolution (float): Resolution in target projection (used for reprojection).
        output_path (str or Path): Output directory for saving results.
        scenario_name (str): Scenario name used for filenames.
        in_millions (bool): Convert results to millions.
        crs (CRS): CRS for saving output raster (optional).
        transform (Affine): Affine transform for saving raster (optional).

    Raises:
        ValueError: If cell size is not provided when required.
        ValueError: If vulnerability or max damage file has invalid structure.

    Returns:
        tuple:
            - pd.DataFrame: Damage per land-use category.
            - np.ndarray: Damage map (grid with estimated damages).
            - np.ndarray: Reprojected land-use map.
            - np.ndarray: Reprojected hazard map.
    """
    
    # load land-use map
    if isinstance(exposure_file, PurePath):
        with rasterio.open(exposure_file) as src:
            landuse = src.read()[0, :, :]
            transform = src.transform
            resolution = src.res[0]
            cellsize = src.res[0] * src.res[1]
    else:
        landuse = exposure_file.copy()

    landuse_in = landuse.copy()

    # Load hazard map
    if isinstance(hazard_file, PurePath):
        if hazard_file.parts[-1].endswith(".tif") | hazard_file.parts[-1].endswith(
            ".tiff"
        ):
            with rasterio.open(hazard_file) as src:
                hazard = src.read()[0, :, :]
                transform = src.transform

        elif hazard_file.parts[-1].endswith(".nc"):
            # Open the hazard netcdf file and store it in the hazard variable
            hazard = xr.open_dataset(hazard_file)

            # Open the landuse geotiff file and store it in the landuse variable
            landuse = xr.open_dataset(exposure_file, engine="rasterio")

            # Match raster to vector
            hazard, landuse = _match_raster_to_vector(
                hazard, landuse, lu_crs, haz_crs, resolution, hazard_col
            )

    elif isinstance(hazard_file, xr.Dataset):
        # Open the landuse geotiff file and store it in the landuse variable
        landuse = xr.open_dataset(exposure_file, engine="rasterio")

        # Match raster to vector
        hazard, landuse = _match_raster_to_vector(
            hazard_file, landuse, lu_crs, haz_crs, resolution, hazard_col
        )

    else:
        hazard = hazard_file.copy()

    # check if land-use and hazard map have the same shape.
    if landuse.shape != hazard.shape:
        warnings.warn(
            "WARNING: landuse and hazard maps are not the same shape. Let's fix this first!"
        )

        landuse, hazard, intersection = match_and_load_rasters(
            exposure_file, hazard_file
        )

        # create the right affine for saving the output
        transform = Affine(
            transform[0],
            transform[1],
            intersection[0],
            transform[3],
            transform[4],
            intersection[1],
        )

    # set cellsize:
    if isinstance(exposure_file, PurePath) | isinstance(hazard_file, PurePath):
        cellsize = src.res[0] * src.res[1]
    else:
        try:
            cellsize = kwargs["cellsize"]
        except KeyError:
            raise ValueError("Required `cellsize` not given.")

    # Load curves
    if isinstance(curve_path, pd.DataFrame):
        curves = curve_path.values
    elif isinstance(curve_path, np.ndarray):
        curves = curve_path
    elif curve_path.parts[-1].endswith(".csv"):
        curves = pd.read_csv(curve_path).values

    if ((curves > 1).all()) or ((curves < 0).all()):
        raise ValueError("Stage-damage curve values must be between 0 and 1")

    # Load maximum damages
    if isinstance(maxdam_path, pd.DataFrame):
        maxdam = maxdam_path.values
    elif isinstance(maxdam_path, np.ndarray):
        maxdam = maxdam_path
    elif maxdam_path.parts[-1].endswith(".csv"):
        maxdam = pd.read_csv(maxdam_path).values

    if maxdam.shape[0] != (curves.shape[1] - 1):
        raise ValueError(
            "Dimensions between maximum damages and the number of depth-damage curve do not agree"
        )

    # Speed up calculation by only considering feasible points
    if kwargs.get("nan_value"):
        nan_value = kwargs.get("nan_value")
        hazard[hazard == nan_value] = 0

    haz = hazard * (hazard >= 0) + 0
    haz[haz >= curves[:, 0].max()] = curves[:, 0].max()
    area = haz > 0
    haz_intensity = haz[haz > 0]
    landuse = landuse[haz > 0]

    # Calculate damage per land-use class for structures
    numberofclasses = len(maxdam)
    alldamage = np.zeros(landuse.shape[0])
    damagebin = np.zeros(
        (
            numberofclasses,
            2,
        )
    )
    for i in range(0, numberofclasses):
        n = maxdam[i, 0]
        damagebin[i, 0] = n
        wd = haz_intensity[landuse == n]
        alpha = np.interp(wd, (curves[:, 0]), curves[:, i + 1])
        damage = alpha * (maxdam[i, 1] * cellsize)
        damagebin[i, 1] = sum(damage)
        alldamage[landuse == n] = damage

    # create the damagemap
    damagemap = np.zeros((area.shape[0], area.shape[1]), dtype=dtype)
    damagemap[area] = alldamage

    # create pandas dataframe with output
    damage_df = (
        pd.DataFrame(damagebin.astype(dtype), columns=["landuse", "damages"])
        .groupby("landuse")
        .sum()
    )

    if save:
        crs = kwargs.get("crs", src.crs)
        transform = kwargs.get("transform", transform)

        # requires adding output_path and scenario_name to function call
        # If output path is not defined, will place file in current directory
        output_path = _check_output_path(kwargs)
        scenario_name = _check_scenario_name(kwargs)
        path_prefix = PurePath(output_path, scenario_name)

        damage_fn = "{}_damages.csv".format(path_prefix)
        damage_df.to_csv(damage_fn)

        dmap_fn = "{}_damagemap.tif".format(path_prefix)
        rst_opts = {
            "driver": "GTiff",
            "height": damagemap.shape[0],
            "width": damagemap.shape[1],
            "count": 1,
            "dtype": dtype,
            "crs": crs,
            "transform": transform,
            "compress": "LZW",
        }
        with rasterio.open(dmap_fn, "w", **rst_opts) as dst:
            dst.write(damagemap, 1)

    if "in_millions" in kwargs:
        damage_df = damage_df / 1e6

    # return output
    return damage_df, damagemap, landuse_in, hazard
