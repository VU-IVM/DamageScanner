# Get all the needed modules
import xarray as xr
import numpy as np
import shapely
import pandas as pd
import geopandas as gpd
import pyproj
from tqdm import tqdm
from pathlib import PurePath
import rasterio
from exactextract import exact_extract
from pyproj import Geod
from shapely.geometry import Point, LineString

import traceback

from damagescanner.osm import read_osm_data
from damagescanner.base_values import DICT_CIS_VULNERABILITY_FLOOD

def _convert_to_meters(feature):
    """Convert coverage to meters for each individual row in the dataframe.

    Args:
        feature (gpd.GeoDataFrame): GeoDataFrame with the hazard and exposure information.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with the hazard and exposure information
        converted to meters.
    """

    line_string = feature.geometry

    # only continue if geometry is a line
    if line_string.geom_type not in ["LineString", "MultiLineString"]:
        return feature.coverage

    if line_string.geom_type == "MultiLineString":
        line_string = line_string[0]

    geod = Geod(ellps="WGS84")
    coverage_meters = []
    for cover in feature.coverage:
        new_for_length = LineString(
            [Point(line_string.coords[0]), line_string.interpolate(cover)]
        )
        coverage_meters.append(geod.geometry_length(new_for_length))

    return coverage_meters

def _get_cell_area_m2(features,hazard_resolution):
    
    geod = Geod(ellps="WGS84")

    asset_point = features.geometry.iloc[0]
    new_geom = asset_point.centroid.buffer(hazard_resolution,cap_style='square')

    line_string = shapely.shortest_line(Point(shapely.get_coordinates(new_geom)[0]),Point(shapely.get_coordinates(new_geom)[1]))


    new_for_length = LineString(
                [Point(line_string.coords[0]), line_string.interpolate(hazard_resolution)]
            )

    resolution = geod.geometry_length(new_for_length)
    
    return resolution*resolution


def _create_grid(bbox, height):
    """Create a vector-based grid

    Args:
        bbox ([type]): [description]
        height ([type]): [description]

    Returns:
        [type]: [description]
    """

    # set xmin,ymin,xmax,and ymax of the grid
    xmin, ymin = shapely.total_bounds(bbox)[0], shapely.total_bounds(bbox)[1]
    xmax, ymax = shapely.total_bounds(bbox)[2], shapely.total_bounds(bbox)[3]

    # estimate total rows and columns
    rows = int(np.ceil((ymax - ymin) / height))
    cols = int(np.ceil((xmax - xmin) / height))

    # set corner points
    x_left_origin = xmin
    x_right_origin = xmin + height
    y_top_origin = ymax
    y_bottom_origin = ymax - height

    # create actual grid
    res_geoms = []
    for countcols in range(cols):
        y_top = y_top_origin
        y_bottom = y_bottom_origin
        for countrows in range(rows):
            res_geoms.append(
                (
                    (
                        (x_left_origin, y_top),
                        (x_right_origin, y_top),
                        (x_right_origin, y_bottom),
                        (x_left_origin, y_bottom),
                    )
                )
            )
            y_top = y_top - height
            y_bottom = y_bottom - height
        x_left_origin = x_left_origin + height
        x_right_origin = x_right_origin + height

    # return grid as shapely polygons
    return shapely.polygons(res_geoms)


def _remove_duplicates(df):
    """Remove duplicates from a dataframe.

    Args:
        df (pd.DataFrame): DataFrame with duplicates.

    Returns:
        pd.DataFrame: DataFrame without duplicates.
    """

    no_duplicates = []
    for row in df.groupby(level=0):
        if len(row[1]) == 1:
            no_duplicates.append(
                [
                    row[0],
                    row[1]["coverage"].values[0],
                    row[1]["values"].values[0],
                ]
            )
        elif len(row[1]) > 1:
            # concatenate coverage column lists
            import itertools

            cov_all = list(itertools.chain.from_iterable(row[1]["coverage"].values))
            # concatenate values column lists
            val_all = list(itertools.chain.from_iterable(row[1]["values"].values))

            # append to no_duplicates with new cov_sum and val_sum
            no_duplicates.append([row[0], cov_all, val_all])

    df_no_duplicates = pd.DataFrame(
        no_duplicates, columns=["index", "coverage", "values"]
    ).set_index("index")

    return df_no_duplicates


def _reproject(hazard, features, hazard_crs):
    """
    Function to reproject to the same CRS in meters.

    Args:


    Returns:
        _type_: _description_
    """

    bounds = features.total_bounds

    bbox = shapely.box(bounds[0], bounds[1], bounds[2], bounds[3])

    centre_point = bbox.centroid
    lat = shapely.get_y(centre_point)
    lon = shapely.get_x(centre_point)

    # formula below based on :https://gis.stackexchange.com/a/190209/80697

    approximate_crs = "EPSG:" + str(
        int(32700 - np.round((45 + lat) / 90, 0) * 100 + np.round((183 + lon) / 6, 0))
    )

    # reproject if needed

    if hazard_crs.to_epsg() == features.crs.to_epsg():
        return hazard, features, approximate_crs

    else:
        features = features.to_crs(hazard_crs)
        return hazard, features, approximate_crs

    # if not pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name == "metre":
    #     hazard = hazard.rio.reproject(approximate_crs)

    # if (
    #     not pyproj.CRS.from_epsg(features.crs.to_epsg()).axis_info[0].unit_name
    #     == "metre"
    # ):
    #     features = features.to_crs(approximate_crs)

    # return hazard, features, approximate_crs


def _overlay_raster_vector(
    hazard, features, hazard_crs, hazard_col="band_data", nodata=-9999, gridded=True, disable_progress=False
):
    """
    Function to overlay a raster with a vector file.

    Args:
        hazard (xarray.DataArray): Raster file with hazard information.
        exposure (gpd.GeoDataFrame): Vector file with exposure information.
        hazard_col (str): Name of the column in the hazard file that contains the hazard intensity.
        hazard_crs (str): CRS of the hazard file.
        nodata (int, optional): Value that indicates no data in the raster file. Defaults to -9999.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with the hazard and exposure information.
    """
    # make sure the hazard data has a crs
    if hazard_crs.to_epsg() == None:
        RuntimeWarning(
            "Hazard crs is not correctly defined. We will now assume it is EPSG:4326"
        )
        hazard = hazard.rio.set_crs("EPSG:4326")
        hazard_crs = pyproj.CRS.from_epsg(4326)

    hazard[hazard_col].rio.write_nodata(nodata, inplace=True)

    area_and_line_objects = features.geom_type.isin(
        ["Polygon", "MultiPolygon", "LineString", "MultiLineString"]
    )

    point_objects = features.geom_type == "Point"

    # make sure they are both in the same coordinate system
    hazard, features, approximate_crs = _reproject(hazard, features, hazard_crs)

    # Check if the exposure data contains any area or line objects
    assert area_and_line_objects.sum() + point_objects.sum() == len(features)

    if point_objects.sum() > 0:
        # if hazard data is a rasterio object:
        if isinstance(hazard, rasterio.io.DatasetReader):
            values = np.array(
                [
                    value[0]
                    for value in rasterio.sample.sample_gen(
                        hazard,
                        [
                            (point.x, point.y)
                            for point in features[point_objects].geometry
                        ],
                    )
                ]
            )

        # if hazard data is a xarray object:
        else:
            hazard_for_points = hazard[hazard_col]
            values = hazard_for_points.sel(
                {
                    hazard_for_points.rio.x_dim: xr.DataArray(
                        features[point_objects].geometry.x
                    ),
                    hazard_for_points.rio.y_dim: xr.DataArray(
                        features[point_objects].geometry.y
                    ),
                },
                method="nearest",
            ).values[0]

        # add values to the features
        features.loc[point_objects, "values"] = values
        features.loc[point_objects, "coverage"] = [1] * len(values)

        # turn values into lists
        features.loc[point_objects, "values"] = features.loc[
            point_objects, "values"
        ].apply(lambda x: [x] if x > 0 else [])

    if not gridded:
        if area_and_line_objects.sum() > 0:
            values_and_coverage_per_area_and_line_object = exact_extract(
                hazard,
                features[area_and_line_objects],
                ["coverage", "values"],
                output="pandas",
            )

            features.loc[area_and_line_objects, "coverage"] = (
                values_and_coverage_per_area_and_line_object["coverage"].values
            )
            features.loc[area_and_line_objects, "values"] = (
                values_and_coverage_per_area_and_line_object["values"].values
            )

            # convert coverage to meters, only do this if the crs is not in meters
            if (
                not pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name
                == "metre"
            ):
                tqdm.pandas(desc="convert coverage to meters",disable=disable_progress)

                features.loc[:, "coverage"] = features.progress_apply(
                    lambda feature: _convert_to_meters(feature), axis=1
                )


    elif gridded:
        if area_and_line_objects.sum() > 0:

            if pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name == "metre":
                grid_cell_size = 50000  # in meters
            else:
                grid_cell_size = 0.5  # in degrees

            # create grid
            bbox = shapely.box(
                hazard.rio.bounds()[0],
                hazard.rio.bounds()[1],
                hazard.rio.bounds()[2],
                hazard.rio.bounds()[3],
            )

            gridded = _create_grid(bbox, grid_cell_size)

            # get all bounds
            all_bounds = gpd.GeoDataFrame(gridded, columns=["geometry"]).bounds

            collect_overlay = []

            # loop over all grids
            for bounds in tqdm(
                all_bounds.itertuples(),
                total=len(all_bounds),
                desc="Overlay raster with vector",disable=disable_progress
            ):
                try:
                    # subset hazard
                    subset_hazard = hazard.rio.clip_box(
                        minx=bounds.minx,
                        miny=bounds.miny,
                        maxx=bounds.maxx,
                        maxy=bounds.maxy,
                    )

                    sub_bbox = shapely.box(
                        bounds.minx, bounds.miny, bounds.maxx, bounds.maxy
                    )

                    # subset features
                    subset_features = gpd.clip(features, sub_bbox)  # list(bounds)[1:])
                    subset_area_and_line_objects = subset_features.geom_type.isin(
                        ["Polygon", "MultiPolygon", "LineString", "MultiLineString"]
                    ).values

                    if len(subset_features) == 0:
                        continue

                    values_and_coverage_per_area_and_line_object = exact_extract(
                        subset_hazard,
                        subset_features[subset_area_and_line_objects],
                        ["coverage", "values"],
                        output="pandas",
                    )

                    # make sure we can connect the results with the features
                    values_and_coverage_per_area_and_line_object.index = (
                        subset_features[subset_area_and_line_objects].index
                    )
                    collect_overlay.append(values_and_coverage_per_area_and_line_object)

                except:
                    get_error = traceback.format_exc()
                    error_to_ignore = "At least one of the clipped raster x,y coordinates has only one point."
                    if error_to_ignore not in get_error:
                        traceback.print_exc()

            df = pd.concat(collect_overlay).sort_index()

            # remove duplicates
            df_no_duplicates = _remove_duplicates(df)

            ## add features to the original features
            features.loc[df_no_duplicates.index, "coverage"] = df_no_duplicates[
                "coverage"
            ]

            features.loc[df_no_duplicates.index, "values"] = df_no_duplicates["values"]

            # Sometimes, with large datasets, a feature may have been excluded from the bbox
            # this has resulted in a null value for the coverage and values. We remove these features.	
            features = features[~features['values'].isnull()]

            # only keep features with values
            features = features[features["values"].apply(lambda x: len(x) > 0)]

            # convert coverage to meters, only do this if the crs is not in meters
            if (
                not pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name
                == "metre"
            ):
                tqdm.pandas(desc="convert coverage to meters",disable=disable_progress)

                features.loc[:, "coverage"] = features.progress_apply(
                    lambda feature: _convert_to_meters(feature), axis=1
                )

    return features


def _overlay_vector_vector(
    hazard, features, hazard_col="band_data", nodata=-9999, gridded=False
):
    """
    Function to overlay a vector with another vector file.

    Args:
        hazard (gpd.GeoDataFrame): Vector file with hazard information.
        exposure (gpd.GeoDataFrame): Vector file with exposure information.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with the hazard and exposure information.

    WIP: This function is not yet finished. It is currently only able to overlay point objects.
    """
    # make sure the hazard data has a crs
    if hazard_crs.to_epsg() == None:
        RuntimeWarning(
            "Hazard crs is not correctly defined. We will now assume it is EPSG:4326"
        )
        hazard = hazard.rio.set_crs("EPSG:4326")
        hazard_crs = pyproj.CRS.from_epsg(4326)

    hazard[hazard_col].rio.write_nodata(nodata, inplace=True)

    area_and_line_objects = features.geom_type.isin(
        ["Polygon", "MultiPolygon", "LineString", "MultiLineString"]
    )
    point_objects = features.geom_type == "Point"

    # Check if the exposure data contains any area or line objects
    assert area_and_line_objects.sum() + point_objects.sum() == len(features)

    if not gridded:
        # reproject if needed
        hazard, features, approximate_crs = _reproject(hazard, features, hazard_crs)

    return features


def _estimate_damage(features, curves, cell_area_m2):
    """
    Function to estimate the damage based on the hazard and exposure information.

    Args:
        exposure (gpd.GeoDataFrame): GeoDataFrame with the hazard and exposure information.
        curves (pd.DataFrame): DataFrame with the stage-damage curves of the different land-use classes.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with the hazard, exposure and damage information.
    """
    
    features["damage"] = features.progress_apply(
        lambda _object: _get_damage_per_object(_object, curves, cell_area_m2),
        axis=1,
    )


    return features


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
        coverage = np.array(asset["coverage"]) * cell_area_m2
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


def VectorExposure(hazard_file, feature_file, asset_type="roads",object_col="object_type", disable_progress=False):
    """
    Function to assess the exposure of objects.

    Args:
        hazard_file (str): Path to the hazard file.
        exposure_file (str): Path to the exposure file.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with the hazard and exposure information.
    """
    # load exposure data
    if isinstance(feature_file, PurePath):
        # if exposure_file is a shapefile, geopackage or parquet file
        if feature_file.suffix in [".shp", ".gpkg", "parquet"]:
            features = gpd.read_file(feature_file)

        # if exposure_file is an osm.pbf file
        elif feature_file.suffix == ".pbf":
            features = read_osm_data(feature_file, asset_type)
            object_col = "object_type"
        else:
            raise ValueError(
                "exposure data should either be a shapefile, geopackage, parquet or osm.pbf file"
            )

    elif isinstance(feature_file, gpd.GeoDataFrame) | isinstance(
        feature_file, pd.DataFrame
    ):
        features = gpd.GeoDataFrame(feature_file.copy())

    else:
        raise ValueError(
            "exposure data should either be a shapefile, geopackage, parquet or osm.pbf file"
        )

    if len(features) == 0:
        hazard_crs = None
        cell_area_m2 = None
        return features, object_col, hazard_crs , cell_area_m2

    # load hazard data
    if isinstance(hazard_file, PurePath):
        if hazard_file.suffix in [".tif", ".tiff", ".nc"]:
            hazard = xr.open_dataset(hazard_file, engine="rasterio")
            hazard_crs = hazard.rio.crs

            # check if crs is already in meters
            if pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name == "metre":
                cell_area_m2 = (hazard.x[1].values - hazard.x[0].values) * (hazard.y[0].values - hazard.y[1].values)
            else:
                cell_area_m2 = _get_cell_area_m2(features,abs(hazard.rio.resolution()[0]))

        elif hazard_file.suffix in [".shp", ".gpkg", ".pbf"]:
            hazard = gpd.read_file(hazard_file)
            hazard_crs = hazard.crs
            cell_area_m2 = 1
        else:
            raise ValueError(
                "hazard data should either be a geotiff, netcdf, shapefile, geopackage or parquet file"
            )
    elif isinstance(hazard_file, rasterio.io.DatasetReader):
        hazard = hazard_file.copy()
        hazard_crs = hazard.crs
        cell_area_m2 = _get_cell_area_m2(features,abs(hazard.res[0]))
    elif isinstance(hazard_file, (xr.Dataset, xr.DataArray)):   
        hazard = hazard_file.copy()
        hazard_crs = hazard.rio.crs

        # check if crs is already in meters
        if pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name == "metre":
            cell_area_m2 = (hazard.x[1].values - hazard.x[0].values) * (hazard.y[0].values - hazard.y[1].values)

        # if not, extract it more cumbersome
        else:
            cell_area_m2 = _get_cell_area_m2(features,abs(hazard.rio.resolution()[0]))

    elif isinstance(hazard, gpd.GeoDataFrame):
        hazard = hazard_file.copy()
        hazard_crs = hazard.crs
        cell_area_m2 = 1
    else:
        raise ValueError(
            f"Hazard should be a raster or GeoDataFrame object, {type(hazard)} given"
        )

    # Run exposure overlay
    if isinstance(hazard, (rasterio.io.DatasetReader, xr.Dataset, xr.DataArray)):
        features = _overlay_raster_vector(hazard, features, hazard_crs, disable_progress=disable_progress)
    elif isinstance(hazard, (gpd.GeoDataFrame, pd.DataFrame)): 
        features = _overlay_vector_vector(hazard, features) ## NOT WORKING YET

    return features, object_col, hazard_crs, cell_area_m2

def VectorScanner(
    hazard_file,
    feature_file,
    curve_path,
    maxdam_path,
    asset_type="roads",
    multi_curves=dict(),
    object_col="object_type",
    disable_progress=False,
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

        *object_col* : Specify the column name of the unique object id's.
        Default is set to **landuse**.

        *sub_types* : List of subtypes that should be included in the analysis.

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

    TO DO: add the option to use vector data for the hazard as well.

    """

    # Load hazard and exposure data, and perform the overlay
    features, object_col, hazard_crs, cell_area_m2 = VectorExposure(
        hazard_file, feature_file, asset_type, object_col, disable_progress
    )

    if len(features) == 0:
        return features

    # Load curves
    if isinstance(curve_path, pd.DataFrame):
        curves = curve_path.copy()
    elif isinstance(curve_path, np.ndarray):
        raise ValueError(
            "For the vector-based approach we use a pandas DataFrame, not a Numpy Array"
        )
    elif curve_path.parts[-1].endswith(".csv"):
        curves = pd.read_csv(curve_path, index_col=[0])

    # Load maximum damages
    if isinstance(maxdam_path, PurePath) and maxdam_path.parts[-1].endswith(".csv"):
        maxdam = pd.read_csv(maxdam_path)
        maxdam = dict(zip(maxdam["object_type"], maxdam["damage"]))
    elif isinstance(maxdam_path, pd.DataFrame):
        maxdam = dict(zip(maxdam_path["object_type"], maxdam_path["damage"]))
    elif isinstance(maxdam_path, np.ndarray):
        maxdam = dict(zip(maxdam_path[:, 0], maxdam_path[:, 1]))
    elif isinstance(maxdam_path, dict):
        maxdam = maxdam_path

    # remove features that are not part of this object type
    if asset_type in DICT_CIS_VULNERABILITY_FLOOD.keys():
        unique_objects_in_asset_type = list(DICT_CIS_VULNERABILITY_FLOOD[asset_type].keys())
        features = features[features['object_type'].isin(unique_objects_in_asset_type)]

    # connect maxdam to exposure
    try:
        features["maximum_damage"] = features.apply(
            lambda x: maxdam[x["object_type"]], axis=1
        )
    except KeyError:
        missing_object_types = [i for i in features.object_type.unique() if i not in maxdam.keys()]
        raise KeyError(
           f"Not all object types in the exposure are included in the maximum damage file: {missing_object_types}"
        )
        
    tqdm.pandas(desc="Calculating damage",disable=disable_progress)

    # Calculate damage
    if not multi_curves:
        features = _estimate_damage(features, curves, cell_area_m2)
    else:
        collect_sub_outcomes  = []
        for curve_id in multi_curves:
            curves = multi_curves[curve_id]
            collect_sub_outcomes.append(_estimate_damage(features, curves, cell_area_m2)['damage'])
        
        all_curve_damages = pd.concat(collect_sub_outcomes,axis=1)
        all_curve_damages.columns = multi_curves.keys()

        # add all curve damages to the features dataframe
        features.loc[:,all_curve_damages.columns] = all_curve_damages

        if 'damage' in features.columns:
            features = features.drop(columns='damage')

    # # Save output
    # if save == True:
    #     # requires adding output_path and scenario_name to function call
    #     # If output path is not defined, will place file in current directory
    #     output_path = _check_output_path(kwargs)
    #     scenario_name = _check_scenario_name(kwargs)
    #     path_prefix = PurePath(output_path, scenario_name)

    #     damage_fn = f'{path_prefix}_damages.csv'
    #     damaged_objects.to_csv(damage_fn)
    #     return damaged_objects

    # else:
    #     return damaged_objects

    return features
