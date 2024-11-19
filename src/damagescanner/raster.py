"""DamageScanner - a directe damage assessment toolkit

Raster specific functions
"""

# Get all the needed modules
from osgeo import gdal
import rasterio
import xarray as xr
import numpy as np
import pandas as pd
from affine import Affine
import warnings
from pathlib import PurePath
import shapely
from damagescanner.utils import _check_output_path, _check_scenario_name


def _match_rasters(raster_in1, raster_in2):
    """
    In case of a mismatch between two rasters, return only the intersecting parts.

    Code adapted from http://sciience.tumblr.com/post/101722591382/finding-the-georeferenced-intersection-between-two

    Arguments:
        *raster_in1* : One of the two rasters to be clipped to the overlapping extent.

        *raster_in2* : One of the two rasters to be clipped to the overlapping extent.

    Returns:
        *array1* : Numpy Array of raster1

        *array2* : Numpy Array of raster2

        *intersection* : Bounding box of overlapping part

    """
    # Read rasterdata
    raster1 = gdal.Open(raster_in1)
    raster2 = gdal.Open(raster_in2)

    # load data
    band1 = raster1.GetRasterBand(1)
    band2 = raster2.GetRasterBand(1)
    gt1 = raster1.GetGeoTransform()
    gt2 = raster2.GetGeoTransform()

    # find each image's bounding box
    # r1 has left, top, right, bottom of dataset's bounds in geospatial coordinates.
    r1 = [
        gt1[0],
        gt1[3],
        gt1[0] + (gt1[1] * raster1.RasterXSize),
        gt1[3] + (gt1[5] * raster1.RasterYSize),
    ]
    r2 = [
        gt2[0],
        gt2[3],
        gt2[0] + (gt2[1] * raster2.RasterXSize),
        gt2[3] + (gt2[5] * raster2.RasterYSize),
    ]
    print("\t1 bounding box: %s" % str(r1))
    print("\t2 bounding box: %s" % str(r2))

    # find intersection between bounding boxes
    intersection = [
        max(r1[0], r2[0]),
        min(r1[1], r2[1]),
        min(r1[2], r2[2]),
        max(r1[3], r2[3]),
    ]
    if r1 != r2:
        print("\t** different bounding boxes **")
        # check for any overlap at all...
        if (intersection[2] < intersection[0]) or (intersection[1] < intersection[3]):
            intersection = None
            print("\t***no overlap***")
            return
        else:
            print("\tintersection:", intersection)
            left1 = int(
                round((intersection[0] - r1[0]) / gt1[1])
            )  # difference divided by pixel dimension
            top1 = int(round((intersection[1] - r1[1]) / gt1[5]))
            col1 = (
                int(round((intersection[2] - r1[0]) / gt1[1])) - left1
            )  # difference minus offset left
            row1 = int(round((intersection[3] - r1[1]) / gt1[5])) - top1

            left2 = int(
                round((intersection[0] - r2[0]) / gt2[1])
            )  # difference divided by pixel dimension
            top2 = int(round((intersection[1] - r2[1]) / gt2[5]))
            col2 = (
                int(round((intersection[2] - r2[0]) / gt2[1])) - left2
            )  # difference minus new left offset
            row2 = int(round((intersection[3] - r2[1]) / gt2[5])) - top2

            # print '\tcol1:',col1,'row1:',row1,'col2:',col2,'row2:',row2
            if col1 != col2 or row1 != row2:
                print("ERROR: Columns and rows still do not match! ***")
            # these arrays should now have the same spatial geometry though NaNs may differ
            array1 = band1.ReadAsArray(left1, top1, col1, row1)
            array2 = band2.ReadAsArray(left2, top2, col2, row2)

    else:  # same dimensions from the get go
        col1 = raster1.RasterXSize  # = col2
        row1 = raster1.RasterYSize  # = row2
        array1 = band1.ReadAsArray()
        array2 = band2.ReadAsArray()

    return array1, array2, intersection

def _match_raster_to_vector(hazard, landuse, lu_crs, haz_crs, resolution, hazard_col):
    """Matches the resolution and extent of a raster to a vector file.

    Arguments:
        *hazard* : netCDF4 with hazard intensity per grid cell.
        *landuse* : netCDF4 with land-use information per grid cell.
        *lu_crs* : EPSG code of the land-use file.
        *haz_crs* : EPSG code of the hazard file.
        *resolution* : Desired resolution of the raster file.
        *hazard_col* : Name of the column in the hazard file that contains the hazard intensity.

    Returns:
        *hazard* : DataSet with hazard intensity per grid cell.

        *landuse* : DataSet with land-use information per grid cell.

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

        landuse, hazard, intersection = _match_rasters(exposure_file, hazard_file)

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
