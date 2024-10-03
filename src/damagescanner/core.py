"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2019 Elco Koks. All versions released under the MIT license.
"""

# Get all the needed modules
import rasterio
import rasterio.sample
import xarray as xr
import numpy as np
import pandas as pd
from affine import Affine
import pyproj
import warnings
from pathlib import PurePath
from exactextract import exact_extract

from damagescanner.vector import (
    match_raster_to_vector,
)
from damagescanner.raster import match_and_load_rasters
from damagescanner.utils import check_output_path, check_scenario_name


def RasterScanner(
    landuse_file,
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
    Raster-based implementation of a direct damage assessment.

    Arguments:
        *landuse_file* : GeoTiff with land-use information per grid cell. Make sure
        the land-use categories correspond with the curves and maximum damages
        (see below). Furthermore, the resolution and extend of the land-use map
        has to be exactly the same as the inundation map.

        *hazard_file* : GeoTiff or netCDF4 with hazard intensity per grid cell. Make sure
        that the unit of the hazard map corresponds with the unit of the
        first column of the curves file.

        *curve_path* : File with the stage-damage curves of the different
        land-use classes. Values should be given as ratios, i.e. between 0 and 1.
        Can also be a pandas DataFrame or numpy Array.

        *maxdam_path* : File with the maximum damages per land-use class
        (in euro/m2). Can also be a pandas DataFrame or numpy Array.

        *dtype*: Set the dtype to the requires precision. This will affect the output damage raster as well

    Optional Arguments:
        *save* : Set to True if you would like to save the output. Requires
        several **kwargs**

    kwargs:
        *nan_value* : if nan_value is provided, will mask the inundation file.
        This option can significantly fasten computations

        *cell_size* : If both the landuse and hazard map are numpy arrays,
        manually set the cell size.

        *resolution* : If landuse is a numpy array, but the hazard map
        is a netcdf, you need to specify the resolution of the landuse map.

        *output_path* : Specify where files should be saved.

        *scenario_name*: Give a unique name for the files that are going to be saved.

        *in_millions*: Set to True if all values should be set in millions.

        *crs*: Specify crs if you only read in two numpy array

        *transform*: Specify transform if you only read in numpy arrays in order to save the result raster

    Raises:
        *ValueError* : on missing kwarg options

    Returns:
     *damagebin* : Table with the land-use class numbers (1st column) and the
     damage for that land-use class (2nd column).

     *damagemap* : Map displaying the damage per grid cell of the area.

    """
    # load land-use map
    if isinstance(landuse_file, PurePath):
        with rasterio.open(landuse_file) as src:
            landuse = src.read()[0, :, :]
            transform = src.transform
            resolution = src.res[0]
            cellsize = src.res[0] * src.res[1]
    else:
        landuse = landuse_file.copy()

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
            landuse = xr.open_dataset(landuse_file, engine="rasterio")

            # Match raster to vector
            hazard, landuse = match_raster_to_vector(
                hazard, landuse, lu_crs, haz_crs, resolution, hazard_col
            )

    elif isinstance(hazard_file, xr.Dataset):
        # Open the landuse geotiff file and store it in the landuse variable
        landuse = xr.open_dataset(landuse_file, engine="rasterio")

        # Match raster to vector
        hazard, landuse = match_raster_to_vector(
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
            landuse_file, hazard_file
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
    if isinstance(landuse_file, PurePath) | isinstance(hazard_file, PurePath):
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
        output_path = check_output_path(kwargs)
        scenario_name = check_scenario_name(kwargs)
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


def _get_damage_per_object(asset, curves, cell_area_m2):
    """
    Calculate damage for a given asset based on hazard information.
    Arguments:
        *asset*: Tuple containing information about the asset. It includes:
            - Index or identifier of the asset (asset[0]).
            - Asset-specific information, including hazard points (asset[1]['hazard_point']).
        *maxdam_dict*: Maximum damage value.
    Returns:
        *tuple*: A tuple containing the asset index or identifier and the calculated damage.
    """

    if asset.geometry.geom_type in ("Polygon", "MultiPolygon"):
        coverage = asset["coverage"] * cell_area_m2
    elif asset.geometry.geom_type in ("LineString", "MultiLineString"):
        coverage = asset["coverage"]
    elif asset.geometry.geom_type in ("Point"):
        coverage = 1
    else:
        raise ValueError(f"Geometry type {asset.geometry.geom_type} not supported")

    return (
        np.sum(
            np.interp(
                asset["values"], curves.index, curves[asset["object_type"]].values
            )
            * coverage
        )
        * asset["maximum_damage"]
    )


def object_scanner(objects, hazard, curves):
    """
    Function to calculate the damage per object.

    Arguments
    ----------
    objects : GeoDataFrame
        GeoDataFrame with the objects for which to calculate the damage. It should contain the following columns: 'object_type' and 'maximum_damage'.
        The maximum damage is the total damage for points, the damage per meter for lines and the damage per square meter for polygons.
    hazard : rasterio.io.DatasetReader or xr.DataArray
        The hazard raster.
    curves : pandas.DataFrame
        The curves to use for the damage calculation. The index should be the values of the hazard raster and the columns should be the object types

    Returns
    -------
    pandas.Series
        Series with the damage per object
    """

    if isinstance(hazard, rasterio.io.DatasetReader):
        hazard_crs = hazard.crs
        cell_area_m2 = hazard.res[0] * hazard.res[1]
    elif isinstance(hazard, (xr.Dataset, xr.DataArray)):
        hazard_crs = hazard.rio.crs
        cell_area_m2 = abs(hazard.rio.resolution()[0]) * abs(hazard.rio.resolution()[1])
    else:
        raise ValueError(f"Hazard should be a raster object, {type(hazard)} given")

    # make sure crs are identical
    assert hazard_crs == objects.crs
    # make sure crs is in meters
    assert pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name == "metre"

    area_and_line_objects = objects.geom_type.isin(
        ["Polygon", "MultiPolygon", "LineString", "MultiLineString"]
    )
    point_objects = objects.geom_type == "Point"

    assert area_and_line_objects.sum() + point_objects.sum() == len(objects)

    if area_and_line_objects.sum() > 0:
        values_and_coverage_per_area_and_line_object = exact_extract(
            hazard,
            objects[area_and_line_objects],
            ["coverage", "values"],
            output="pandas",
        )

        objects.loc[area_and_line_objects, "coverage"] = (
            values_and_coverage_per_area_and_line_object["coverage"].values
        )
        objects.loc[area_and_line_objects, "values"] = (
            values_and_coverage_per_area_and_line_object["values"].values
        )

    if point_objects.sum() > 0:
        if isinstance(hazard, rasterio.io.DatasetReader):
            values = np.array(
                [
                    value[0]
                    for value in rasterio.sample.sample_gen(
                        hazard,
                        [
                            (point.x, point.y)
                            for point in objects[point_objects].geometry
                        ],
                    )
                ]
            )
        else:
            values = hazard.sel(
                {
                    hazard.rio.x_dim: xr.DataArray(objects[point_objects].geometry.x),
                    hazard.rio.y_dim: xr.DataArray(objects[point_objects].geometry.y),
                },
                method="nearest",
            ).values[0]
        objects.loc[point_objects, "values"] = values

    damage = objects.apply(
        lambda _object: _get_damage_per_object(_object, curves, cell_area_m2),
        axis=1,
    )
    return damage
