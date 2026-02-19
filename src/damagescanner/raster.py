"""DamageScanner - a directe damage assessment toolkit.

Raster specific functions
"""

import warnings
from pathlib import Path, PurePath
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import shapely
import xarray as xr
from affine import Affine
from rasterio.windows import Window

from damagescanner.utils import _check_output_path, _check_scenario_name


def match_and_load_rasters(
    raster_in1: str | Path, raster_in2: str | Path
) -> tuple[np.ndarray, np.ndarray, Affine]:
    """
    Match and clip two raster files to their common spatial extent and resolution.

    Code adapted from http://sciience.tumblr.com/post/101722591382/finding-the-georeferenced-intersection-between-two


    Args:
        raster_in1: Path to the first raster file.
        raster_in2: Path to the second raster file.

    Returns:
        tuple: (data1, data2, transform)
            - Clipped raster array from the first file.
            - Clipped raster array from the second file.
            - Affine transform of the intersecting region.

    Raises:
        ValueError: If the rasters have different CRS or resolution.
    """
    with rasterio.open(raster_in1) as src1, rasterio.open(raster_in2) as src2:
        if src1.crs != src2.crs:
            raise ValueError("Different CRS: CRS must be the same.")
        if src1.res != src2.res:
            raise ValueError("Different resolution: Cell sizes must be the same.")

        top_delta: int = round((src2.bounds.top - src1.bounds.top) / src1.transform.e)
        bottom_delta: int = round(
            (src2.bounds.bottom - src1.bounds.bottom) / src1.transform.e
        )
        left_delta: int = round(
            (src2.bounds.left - src1.bounds.left) / src1.transform.a
        )
        right_delta: int = round(
            (src2.bounds.right - src1.bounds.right) / src1.transform.a
        )

        data1 = src1.read(
            1,
            window=Window(
                col_off=left_delta,  # ty:ignore[unknown-argument]
                row_off=top_delta,  # ty:ignore[unknown-argument]
                width=src1.width - left_delta + right_delta,  # ty:ignore[unknown-argument]
                height=src1.height - top_delta + bottom_delta,  # ty:ignore[unknown-argument]
            ),
        )
        data2 = src2.read(
            1,
            window=Window(
                col_off=abs(min(left_delta, 0)),  # ty:ignore[unknown-argument]
                row_off=abs(min(top_delta, 0)),  # ty:ignore[unknown-argument]
                width=max(src1.width, src2.width) - abs(left_delta) - abs(right_delta),  # ty:ignore[unknown-argument]
                height=max(src1.height, src2.height)  # ty:ignore[unknown-argument]
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


def _match_raster_to_vector(
    hazard: xr.Dataset,
    landuse: xr.Dataset,
    lu_crs: int,
    haz_crs: int,
    resolution: float,
    hazard_col: str,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Align hazard and land-use rasters by CRS, extent, and resolution.

    Args:
        hazard: Hazard dataset.
        landuse: Land-use dataset.
        lu_crs: EPSG code of the land-use CRS.
        haz_crs: EPSG code of the hazard CRS.
        resolution: Target resolution in projected units.
        hazard_col: Column name containing hazard intensity.

    Returns:
        tuple: (hazard_array, landuse_array)
            - Reprojected hazard raster as array.
            - Reprojected land-use raster as array.
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
    exposure_file: Path | np.ndarray | str,
    hazard_file: Path | np.ndarray | xr.Dataset | str,
    curve_path: Path | pd.DataFrame | np.ndarray | str,
    maxdam_path: Path | pd.DataFrame | np.ndarray | str,
    lu_crs: int = 28992,
    haz_crs: int = 4326,
    hazard_col: str = "FX",
    dtype: type = np.int32,
    save: bool = False,
    **kwargs: Any,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, xr.Dataset | np.ndarray]:
    """
    Run a raster-based direct damage assessment using hazard and exposure layers.

    Args:
        exposure_file: Path to land-use GeoTIFF or numpy array.
        hazard_file: Path to hazard raster or dataset.
        curve_path: Vulnerability curve(s).
        maxdam_path: Maximum damage values.
        lu_crs: CRS of the land-use file (default EPSG:28992).
        haz_crs: CRS of the hazard file (default EPSG:4326).
        hazard_col: Column containing hazard intensity (default "FX").
        dtype: Output dtype for damage raster.
        save: If True, saves damage results to file.
        kwargs: Additional keyword arguments for saving output and other options.

    Keyword Args:
        nan_value: Replace this value in the hazard raster with 0.
        cellsize: Cell size (m²) if exposure and hazard are arrays.
        resolution: Resolution in target projection (used for reprojection).
        output_path: Output directory for saving results.
        scenario_name: Scenario name used for filenames.
        in_millions: Convert results to millions.
        crs: CRS for saving output raster (optional).
        transform: Affine transform for saving raster (optional).

    Raises:
        ValueError: If cell size is not provided when required.
        ValueError: If vulnerability or max damage file has invalid structure.

    Returns:
        tuple: (damage_df, damagemap, landuse_in, hazard)
            - Damage per land-use category.
            - Damage map (grid with estimated damages).
            - Reprojected land-use map.
            - Reprojected hazard map.
    """
    if isinstance(curve_path, str):
        curve_path = Path(curve_path)
    if isinstance(maxdam_path, str):
        maxdam_path = Path(maxdam_path)
    if isinstance(exposure_file, str):
        exposure_file = Path(exposure_file)
    if isinstance(hazard_file, str):
        hazard_file = Path(hazard_file)

    # Initialize metadata
    has_raster_metadata = (
        isinstance(exposure_file, Path)
        or (
            isinstance(hazard_file, Path)
            and hazard_file.suffix.lower() in [".tif", ".tiff"]
        )
        or isinstance(hazard_file, xr.Dataset)
    )

    if has_raster_metadata:
        if "cellsize" in kwargs or "transform" in kwargs:
            raise ValueError(
                "cellsize and transform are loaded from the raster files or datasets, please do not set them as keyword arguments."
            )
        cellsize = None
        transform = None
    else:
        cellsize = kwargs.get("cellsize")
        transform = kwargs.get("transform")
        if cellsize is None:
            raise ValueError(
                "When using array inputs for exposure and hazard, you must provide the cell size (in m²) as a keyword argument."
            )
        if transform is None:
            raise ValueError(
                "When using array inputs for exposure and hazard, you must provide an affine transform as a keyword argument."
            )

    crs = kwargs.get("crs")
    resolution = kwargs.get("resolution")

    # load land-use map
    if isinstance(exposure_file, Path):
        with rasterio.open(exposure_file) as src:
            landuse = src.read()[0, :, :]
            if transform is None:
                transform = src.transform
            if resolution is None:
                resolution = src.res[0]
            if cellsize is None:
                cellsize = abs(src.res[0] * src.res[1])
            if crs is None:
                crs = src.crs
    else:
        landuse = exposure_file.copy()

    landuse_in = landuse.copy()

    # Load hazard map
    if isinstance(hazard_file, Path):
        if hazard_file.suffix.lower() in [".tif", ".tiff"]:
            with rasterio.open(hazard_file) as src_haz:
                hazard = src_haz.read()[0, :, :]
                if transform is None:
                    transform = src_haz.transform
                if crs is None:
                    crs = src_haz.crs
                if cellsize is None:
                    cellsize = abs(src_haz.res[0] * src_haz.res[1])

        elif hazard_file.suffix.lower() == ".nc":
            hazard = xr.open_dataset(hazard_file)
    else:
        hazard = hazard_file.copy()

    # Extract metadata from hazard if it is an xarray dataset
    if isinstance(hazard, xr.Dataset):
        try:
            if resolution is None:
                resolution = hazard.rio.resolution()[0]
            if cellsize is None:
                res = hazard.rio.resolution()
                cellsize = abs(res[0] * res[1])
            if transform is None:
                transform = hazard.rio.transform()
            if crs is None:
                crs = hazard.rio.crs
        except Exception:
            # If rioxarray fails to get metadata (e.g. non-spatial NC), continue
            pass

    # Align hazard and land-use if hazard is an xarray dataset
    if isinstance(hazard, xr.Dataset):
        # Open land-use as dataset if it isn't one already for aligning
        if isinstance(exposure_file, Path):
            landuse_ds = xr.open_dataset(exposure_file, engine="rasterio")
        else:
            landuse_ds = landuse

        # Match hazard and land-use
        hazard, landuse = _match_raster_to_vector(
            hazard, landuse_ds, lu_crs, haz_crs, resolution, hazard_col
        )

    # check if land-use and hazard map have the same shape.
    if landuse.shape != hazard.shape:
        warnings.warn(
            "WARNING: landuse and hazard maps are not the same shape. Let's fix this first!"
        )

        landuse, hazard, intersection = match_and_load_rasters(
            exposure_file, hazard_file
        )

        # create the right affine for saving the output
        if transform is not None:
            transform = Affine(
                transform[0],
                transform[1],
                intersection[0],
                transform[3],
                transform[4],
                intersection[1],
            )

    # set cellsize:
    if cellsize is None:
        raise ValueError("Required `cellsize` not given.")

    # Load curves
    curves: np.ndarray
    if isinstance(curve_path, pd.DataFrame):
        curves: np.ndarray = curve_path.values
    elif isinstance(curve_path, np.ndarray):
        curves: np.ndarray = curve_path
    elif curve_path.parts[-1].endswith(".csv"):
        curves: np.ndarray = pd.read_csv(curve_path).values
    else:
        raise ValueError(
            "Invalid curve file format. Must be DataFrame, ndarray, or CSV."
        )

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
        pd.DataFrame(damagebin.astype(dtype), columns=["landuse", "damage"])
        .groupby("landuse")
        .sum()
    )

    if save:
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
