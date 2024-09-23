"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2019 Elco Koks. All versions released under the MIT license.
"""

# Get all the needed modules
import rasterio
import xarray as xr
import numpy as np
import shapely
import pandas as pd
import geopandas as gpd
from affine import Affine
import pyproj
from tqdm import tqdm
import warnings
from pathlib import PurePath
from exactextract import exact_extract

from damagescanner.vector import (
    reproject,
    get_damage_per_object,
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


def calculate_damage_per_object(objects, hazard, curves, maximum_damage, column):
    """
    Function to calculate the damage per object.
    
    Arguments:
        *objects* : GeoDataFrame with the objects for which we want to calculate the damage.
        *hazard* : RasterObject (e.g., xarray dataset or array with crs or rasterio dataset with crs) with the hazard intensity per grid cell.
        *curves* : DataFrame with the vulnerability curves.
        *maximum_damage* : DataFrame with the maximum damage per object type.
        *column* : The column name that specifies the object type. This is used to get the
        associated maximum damage and curve for the object type.

    Returns:
        *damage* : The damage per object.
        
    """
    if isinstance(hazard, rasterio.io.DatasetReader):
        hazard_crs = hazard.crs
    elif isinstance(hazard, (xr.Dataset, xr.DataArray)):
        hazard_crs = hazard.rio.crs
    else:
        raise ValueError(f"Hazard should be a raster object, {type(hazard)} given")
    
    # make sure crs are identical
    assert hazard_crs == objects.crs
    # make sure crs is in meters
    assert pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name == "metre"
    
    # reset the index to make sure we can use the index to assign the damage in the end
    objects = objects.reset_index(drop=True)

    # get the mean severity of the hazard per object
    objects['severity'] = exact_extract(
        hazard, objects, ["mean"], output="pandas"
    )['mean']

    # create an array to store the damage
    damage = np.full_like(objects['severity'], fill_value=np.nan, dtype=np.float64)
    
    # iterate over the unique object types and calculate the damage for each object
    for group_label, group_data in objects.groupby(column):
        # get the curves for the specific object type
        group_curves = curves[group_label]
        # the the maximum damage for the specific object type
        group_maxdam = maximum_damage.loc[group_label, 'damage']
        # calculate the damage ratio. On the left of the curve we assume no damage, on the right we assume the maximum damage
        group_damage_ratio = np.interp(
            group_data['severity'], group_curves.index, group_curves.values, left=0, right=group_curves.values[-1]
        ) 
        # calculate and assign damage
        group_damage = group_damage_ratio * group_maxdam * group_data.area
        damage[group_data.index] = group_damage
    return damage



def VectorScanner(
    exposure_file,
    hazard_file,
    curve_path,
    maxdam_path,
    cell_size=5,
    exp_crs=4326,
    haz_crs=4326,
    object_col="landuse",
    hazard_col="inun_val",
    centimeters=False,
    save=False,
    **kwargs,
):
    """
    Vector based implementation of a direct damage assessment

    Arguments:
        *exposure_file* : Shapefile, Pandas DataFrame or Geopandas GeoDataFrame
        with land-use information of the area.

        *hazard_file* : GeoTiff with inundation depth per grid cell. Make sure
        that the unit of the inundation map corresponds with the unit of the
        first column of the curves file.

        *curve_path* : File with the stage-damage curves of the different
        land-use classes. Can also be a pandas DataFrame (but not a numpy Array).

        *maxdam_path* : File with the maximum damages per land-use class
        (in euro/m2). Can also be a pandas DataFrame (but not a numpy Array).

    Optional Arguments:
        *cell_size* : Specify the cell size of the hazard map.

        *exp_crs* : Specify the coordinate reference system of the exposure. Default is
        set to **4326**. A preferred CRS system is in meters. Please note that the function
        only accepts the CRS system in the form of an integer EPSG code.

        *haz_crs* : Specify the cooordinate reference system of the hazard. Default is
        set to **4326**.  A preferred CRS system is in meters.Please note that the function
        only accepts the CRS system in the form of an integer EPSG code.

        *centimeters* : Set to True if the inundation map and curves are in
        centimeters

        *object_col* : Specify the column name of the unique object id's.
        Default is set to **landuse**.

        *hazard_col* : Specify the column name of the hazard intensity
        Default is set to **inun_val**.

        *save* : Set to True if you would like to save the output. Requires
        several **kwargs**

    kwargs:
        *output_path* : Specify where files should be saved.

        *scenario_name*: Give a unique name for the files that are going to be saved.

    Raises:
        *ValueError* : on missing kwargs

    Returns:
     *damagebin* : Table with the land-use class names (1st column) and the
     damage for that land-use class (2nd column).

    """
    # load exposure data
    if isinstance(exposure_file, PurePath):
        exposure = gpd.read_file(exposure_file)

    elif isinstance(exposure_file, gpd.GeoDataFrame) | isinstance(
        exposure_file, pd.DataFrame
    ):
        exposure = gpd.GeoDataFrame(exposure_file.copy())

    else:
        print(
            "ERROR: exposure data should either be a shapefile, GeoPackage, a GeoDataFrame or a pandas Dataframe with a geometry column"
        )

    # load hazard file
    if isinstance(hazard_file, PurePath):
        if hazard_file.parts[-1].endswith(".tif") | hazard_file.parts[-1].endswith(
            ".tiff"
        ):
            # load dataset
            hazard_map = xr.open_dataset(hazard_file, engine="rasterio")

            # specify hazard_col
            hazard_col = "band_data"

            # convert to dataframe
            hazard = hazard_map["band_data"].to_dataframe().reset_index()

            # drop all non values and zeros to reduce size
            hazard = hazard.loc[
                ~(hazard["band_data"].isna() | hazard["band_data"] <= 0)
            ].reset_index(drop=True)

            # create geometry values and drop lat lon columns
            hazard["geometry"] = shapely.points(
                np.array(list(zip(hazard["x"], hazard["y"])))
            )
            hazard = hazard.drop(["x", "y", "band", "spatial_ref"], axis=1)

            # and turn them into squares again:
            hazard.geometry = shapely.buffer(
                hazard.geometry, distance=cell_size / 2, cap_style="square"
            ).values

        elif hazard_file.parts[-1].endswith(".nc"):
            # load dataset
            hazard_map = xr.open_dataset(hazard_file)

            # convert to dataframe
            hazard = hazard_map[hazard_col].to_dataframe().reset_index()

            # drop all non values and below zeros to reduce size. Might cause issues for cold waves
            hazard = hazard.loc[
                ~(hazard[hazard_col].isna() | hazard[hazard_col] <= 0)
            ].reset_index(drop=True)

            # create geometry values and drop lat lon columns
            hazard["geometry"] = shapely.points(
                np.array(list(zip(hazard["x"], hazard["y"])))
            )
            hazard = hazard.drop(["x", "y", "band", "spatial_ref"], axis=1)

            # and turn them into squares again:
            hazard.geometry = shapely.buffer(
                hazard.geometry, distance=cell_size / 2, cap_style="square"
            ).values

        elif hazard_file.parts[-1].endswith(".shp") | hazard_file.parts[-1].endswith(
            ".gpkg"
        ):
            hazard = gpd.read_file(hazard_file)
        else:
            print(
                "ERROR: hazard data should either be a GeoTIFF, a netCDF4, a shapefile or a GeoPackage."
            )

    elif isinstance(hazard_file, gpd.GeoDataFrame) | isinstance(
        hazard_file, pd.DataFrame
    ):
        hazard = gpd.GeoDataFrame(hazard_file.copy())
    else:
        raise ValueError(
            "ERROR: hazard data should be a GeoTiff, a netCDF4, a shapefile, a GeoDataFrame \
              or any other georeferenced format that can be read by xarray or geopandas"
        )

    print("PROGRESS: Exposure and hazard data loaded")

    # check if exposure and hazard data are in the same CRS and in a CRS in meters
    CRS_exposure = pyproj.CRS.from_epsg(exp_crs)
    CRS_hazard = pyproj.CRS.from_epsg(haz_crs)

    # get the unit name of the CRS
    exp_crs_unit_name = CRS_exposure.axis_info[0].unit_name
    haz_crs_unit_name = CRS_hazard.axis_info[0].unit_name

    # reproject exposure and hazard data to the same CRS
    if exp_crs_unit_name == "degree" and haz_crs_unit_name == "metre":
        exposure.geometry = reproject(
            exposure, current_crs=f"epsg:{exp_crs}", approximate_crs=f"epsg:{haz_crs}"
        )
    elif exp_crs_unit_name == "degree" and haz_crs_unit_name == "degree":
        exposure.geometry = reproject(
            exposure, current_crs=f"epsg:{exp_crs}", approximate_crs="epsg:3857"
        )

    if haz_crs_unit_name == "degree" and exp_crs_unit_name == "metre":
        hazard.geometry = reproject(
            hazard, current_crs=haz_crs, approximate_crs=f"epsg:{exp_crs}"
        )
    elif haz_crs_unit_name == "degree" and exp_crs_unit_name == "degree":
        hazard.geometry = reproject(
            hazard, current_crs=haz_crs, approximate_crs="epsg:3857"
        )

    # rename inundation colum, set values to centimeters if required and to integers:
    hazard = hazard.rename(columns={hazard_col: "haz_val"})
    if not centimeters:
        hazard["haz_val"] = hazard.haz_val * 100
        hazard["haz_val"] = hazard.haz_val.astype(int)

    # rename object col to make sure everything is consistent:
    exposure = exposure.rename({object_col: "obj_type"}, axis=1)

    # Load curves
    if isinstance(curve_path, pd.DataFrame):
        curves = curve_path.copy()
    elif isinstance(curve_path, np.ndarray):
        raise ValueError(
            "ERROR: for the vector-based approach we use a pandas DataFrame, not a Numpy Array"
        )
    elif curve_path.parts[-1].endswith(".csv"):
        curves = pd.read_csv(curve_path, index_col=[0])

    # Load maximum damages
    if isinstance(maxdam_path, PurePath) and maxdam_path.parts[-1].endswith(".csv"):
        maxdam_path = pd.read_csv(maxdam_path)
    elif isinstance(maxdam_path, pd.DataFrame):
        maxdam = dict(zip(maxdam_path[object_col], maxdam_path["damage"]))
    elif isinstance(maxdam_path, np.ndarray):
        maxdam = dict(zip(maxdam_path[:, 0], maxdam_path[:, 1]))
    elif isinstance(maxdam_path, dict):
        maxdam = maxdam_path

    # overlay hazard and exposure data
    hazard_tree = shapely.STRtree(hazard.geometry.values)

    if (shapely.get_type_id(exposure.iloc[0].geometry) == 3) | (
        shapely.get_type_id(exposure.iloc[0].geometry) == 6
    ):
        overlay = hazard_tree.query(exposure.geometry, predicate="intersects")
    else:
        overlay = hazard_tree.query(
            shapely.buffer(exposure.geometry.values, distance=cell_size * 2),
            predicate="intersects",
        )

    overlay_exp_haz = pd.DataFrame(overlay.T, columns=["obj_type", "hazard_point"])

    # Perform calculation
    collect_output = []

    for obj in tqdm(
        overlay_exp_haz.groupby("obj_type"),
        total=len(overlay_exp_haz.obj_type.unique()),
        desc="damage calculation",
    ):
        collect_output.append(
            get_damage_per_object(obj, hazard, exposure, curves, maxdam)
        )

    # Merge results
    damaged_objects = exposure.merge(
        pd.DataFrame(collect_output, columns=["index", "damage"]),
        left_index=True,
        right_on="index",
    )[["obj_type", "geometry", "damage"]]

    # Save output
    if save is True:
        # requires adding output_path and scenario_name to function call
        # If output path is not defined, will place file in current directory
        output_path = check_output_path(kwargs)
        scenario_name = check_scenario_name(kwargs)
        path_prefix = PurePath(output_path, scenario_name)

        damage_fn = f"{path_prefix}_damages.csv"
        damaged_objects.to_csv(damage_fn)
        return damaged_objects

    else:
        return damaged_objects
