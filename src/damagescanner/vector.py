# Get all the needed modules
import xarray as xr
import numpy as np
import shapely
import pandas as pd
import geopandas as gpd
import pyproj
from pathlib import Path
from tqdm import tqdm
from pathlib import PurePath
import rasterio
from exactextract import exact_extract
from pyproj import Geod
from shapely.geometry import Point, LineString

import traceback

from damagescanner.osm import read_osm_data
from damagescanner.config import DICT_CIS_VULNERABILITY_FLOOD


def _convert_to_meters(feature):
    """
    Convert coverage length to meters for each row with LineString geometries.

    Args:
        feature (gpd.GeoSeries): A GeoSeries row containing geometry and coverage fields.

    Returns:
        list: Coverage values in meters.
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


def _get_cell_area_m2(features, hazard_resolution):
    """
    Estimate the area (m²) of a raster grid cell using the feature centroid and resolution.

    Args:
        features (gpd.GeoDataFrame): Feature set used to center the buffer.
        hazard_resolution (float): Hazard raster resolution (in degrees).

    Returns:
        float: Grid cell area in square meters.
    """
    geod = Geod(ellps="WGS84")

    asset_point = features.geometry.iloc[0]
    new_geom = asset_point.centroid.buffer(hazard_resolution, cap_style="square")

    line_string = shapely.shortest_line(
        Point(shapely.get_coordinates(new_geom)[0]),
        Point(shapely.get_coordinates(new_geom)[1]),
    )

    new_for_length = LineString(
        [Point(line_string.coords[0]), line_string.interpolate(hazard_resolution)]
    )

    resolution = geod.geometry_length(new_for_length)

    return resolution * resolution


def _create_grid(bbox, height):
    """
    Create a regular vector grid over a bounding box.

    Args:
        bbox (shapely.Geometry): The bounding box to cover.
        height (float): Cell height (assumed square cells).

    Returns:
        list: List of shapely Polygon grid cells.
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
    """
    Merge duplicate features by concatenating coverage and values.

    Args:
        df (pd.DataFrame): DataFrame with possible duplicate index entries.

    Returns:
        pd.DataFrame: De-duplicated DataFrame.
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
    Reproject features to match hazard CRS. Also infers a suitable UTM CRS.

    Args:
        hazard (xr.Dataset): Hazard dataset.
        features (gpd.GeoDataFrame): Exposure features.
        hazard_crs (pyproj.CRS): Coordinate system of hazard layer.

    Returns:
        tuple: (hazard, reprojected_features, approximate_crs)
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
    hazard,
    features,
    hazard_crs,
    nodata=-9999,
    gridded=True,
    disable_progress=False,
):
    """
    Overlay raster hazard data on vector exposure features.

    Args:
        hazard (xr.Dataset | rasterio.io.DatasetReader): Raster hazard layer.
        features (gpd.GeoDataFrame): Vector exposure features.
        hazard_crs (pyproj.CRS): CRS of hazard data.
        nodata (int): No-data value in raster.
        gridded (bool): Whether to process in spatial chunks.
        disable_progress (bool): Disable tqdm progress bar.

    Returns:
        gpd.GeoDataFrame: Exposure features with added `coverage` and `values`.
    """

    # make sure the hazard data has a crs
    if hazard_crs.to_epsg() is None:
        RuntimeWarning(
            "Hazard crs is not correctly defined. We will now assume it is EPSG:4326"
        )
        hazard = hazard.rio.set_crs("EPSG:4326")
        hazard_crs = pyproj.CRS.from_epsg(4326)

    hazard.rio.write_nodata(nodata, inplace=True)

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
            values = hazard.sel(
                {
                    hazard.rio.x_dim: xr.DataArray(features[point_objects].geometry.x),
                    hazard.rio.y_dim: xr.DataArray(features[point_objects].geometry.y),
                },
                method="nearest",
            ).values

        # add values to the features
        features.loc[point_objects, "values"] = values
        features.loc[point_objects, "coverage"] = [1] * len(values)

        # turn values into lists
        features.loc[point_objects, "values"] = features.loc[
            point_objects, "values"
        ].apply(lambda x: [x] if x > 0 else [])

    exact_extract_kwargs = {
        "output": "pandas",
        "include_geom": False,
        "strategy": "raster-sequential",
    }

    if not gridded:
        if area_and_line_objects.sum() > 0:
            values_and_coverage_per_area_and_line_object = exact_extract(
                hazard,
                features[area_and_line_objects][["geometry"]],  # only pass the geometry
                ["coverage", "values"],
                **exact_extract_kwargs,
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
                tqdm.pandas(desc="convert coverage to meters", disable=disable_progress)

                features.loc[:, "coverage"] = features.progress_apply(
                    lambda feature: _convert_to_meters(feature), axis=1
                )

    elif gridded:
        if area_and_line_objects.sum() > 0:
            if (
                pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name
                == "metre"
            ):
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
                desc="Overlay raster with vector",
                disable=disable_progress,
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
                        subset_features[subset_area_and_line_objects][
                            ["geometry"]
                        ],  # only pass the geometry
                        ["coverage", "values"],
                        **exact_extract_kwargs,
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
            features = features[~features["values"].isnull()]

            # only keep features with values
            features = features[features["values"].apply(lambda x: len(x) > 0)]

            # convert coverage to meters, only do this if the crs is not in meters
            if (
                not pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name
                == "metre"
            ):
                tqdm.pandas(desc="convert coverage to meters", disable=disable_progress)

                features.loc[:, "coverage"] = features.progress_apply(
                    lambda feature: _convert_to_meters(feature), axis=1
                )

    return features


def _overlay_vector_vector(hazard, features, nodata=-9999, gridded=False):
    """
    Overlay a vector hazard layer onto vector exposure features.

    Args:
        hazard (gpd.GeoDataFrame): Hazard vector features.
        features (gpd.GeoDataFrame): Exposure vector features.
        nodata (int): No-data placeholder.
        gridded (bool): Chunk processing toggle (not implemented).

    Returns:
        gpd.GeoDataFrame: Original feature set (currently not modified).

    WIP: This function is not yet finished. It is currently only able to overlay point objects.

    """
    # make sure the hazard data has a crs
    if hazard.crs.to_epsg() is None:
        RuntimeWarning(
            "Hazard crs is not correctly defined. We will now assume it is EPSG:4326"
        )
        hazard = hazard.rio.set_crs("EPSG:4326")
        hazard_crs = pyproj.CRS.from_epsg(4326)

    hazard.rio.write_nodata(nodata, inplace=True)

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
    Estimate total damage per asset using vulnerability curves.

    Args:
        features (gpd.GeoDataFrame): Exposure with hazard info.
        curves (pd.DataFrame): Vulnerability curve per asset type.
        cell_area_m2 (float): Area per grid cell in m² (for polygons).

    Returns:
        gpd.GeoDataFrame: Exposure data with added `damage` column.
    """
    features["damage"] = features.progress_apply(
        lambda _object: _get_damage_per_object(_object, curves, cell_area_m2),
        axis=1,
    )

    return features


def _get_damage_per_object(asset, curves, cell_area_m2):
    """
    Compute damage for a single asset using interpolated vulnerability.

    Args:
        asset (pd.Series): Single feature with geometry, values, object_type.
        curves (pd.DataFrame): Vulnerability curves.
        cell_area_m2 (float | int): Cell area in square meters.

    Returns:
        float: Estimated damage value.
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


def VectorExposure(
    hazard_file,
    feature_file,
    asset_type="roads",
    object_col="object_type",
    disable_progress=False,
    gridded: bool = True,
):
    """
    Load and overlay vector or raster hazard with vector exposure data.

    Args:
        hazard_file (Path | xr.Dataset | rasterio.DatasetReader): Hazard input.
        feature_file (Path | GeoDataFrame | pd.DataFrame): Exposure input.
        asset_type (str): Infrastructure category (only for OSM).
        object_col (str): Name of the object type column.
        disable_progress (bool): Whether to suppress progress bars.
        gridded (bool): Whether to process in spatial chunks.

    Returns:
        tuple: (features, object_col, hazard_crs, cell_area_m2)
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
        return features, object_col, hazard_crs, cell_area_m2

    # load hazard data
    if isinstance(hazard_file, PurePath):
        if hazard_file.suffix in [".tif", ".tiff", ".nc"]:
            hazard = xr.open_dataset(hazard_file, engine="rasterio")
            hazard_crs = hazard.rio.crs

            # check if crs is already in meters
            if (
                pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name
                == "metre"
            ):
                cell_area_m2 = abs(
                    (hazard.x[1].values - hazard.x[0].values)
                    * (hazard.y[0].values - hazard.y[1].values)
                )
            else:
                cell_area_m2 = _get_cell_area_m2(
                    features, abs(hazard.rio.resolution()[0])
                )

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
        cell_area_m2 = _get_cell_area_m2(features, abs(hazard.res[0]))
    elif isinstance(hazard_file, (xr.Dataset, xr.DataArray)):
        hazard = hazard_file.copy()
        hazard_crs = hazard.rio.crs

        # check if crs is already in meters
        if pyproj.CRS.from_epsg(hazard_crs.to_epsg()).axis_info[0].unit_name == "metre":
            cell_area_m2 = abs(
                (hazard.x[1].values - hazard.x[0].values)
                * (hazard.y[0].values - hazard.y[1].values)
            )

        # if not, extract it more cumbersome
        else:
            cell_area_m2 = _get_cell_area_m2(features, abs(hazard.rio.resolution()[0]))

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
        features = _overlay_raster_vector(
            hazard,
            features,
            hazard_crs,
            disable_progress=disable_progress,
            gridded=gridded,
        )
    elif isinstance(hazard, (gpd.GeoDataFrame, pd.DataFrame)):
        features = _overlay_vector_vector(
            hazard, features, gridded=gridded
        )  ## NOT WORKING YET

    return features, object_col, hazard_crs, cell_area_m2


def VectorScanner(
    hazard_file,
    feature_file,
    curve_path,
    maxdam_path: Path | pd.DataFrame | dict | None = None,
    asset_type: str | None = None,
    multi_curves: dict = dict(),
    object_col="object_type",
    disable_progress=False,
    gridded: bool = True,
    **kwargs,
):
    """
    Perform vector-based direct damage assessment using hazard and exposure layers.

    Args:
        hazard_file (Path | xr.Dataset | gpd.GeoDataFrame): Hazard input.
        feature_file (Path | gpd.GeoDataFrame): Exposure input.
        curve_path (Path | pd.DataFrame): Vulnerability curve(s).
        maxdam_path (Path | pd.DataFrame | dict): Maximum damage values.
        asset_type (str): Infrastructure class (only for OSM).
        multi_curves (dict, optional): Multiple curve sets.
        object_col (str): Column name with object type.
        disable_progress (bool): Whether to suppress progress bars.
        gridded (bool): Whether to process in spatial chunks.

    Returns:
        gpd.GeoDataFrame: Exposure data with calculated damages.
    """
    # Load hazard and exposure data, and perform the overlay
    features, object_col, hazard_crs, cell_area_m2 = VectorExposure(
        hazard_file,
        feature_file,
        asset_type,
        object_col,
        disable_progress,
        gridded=gridded,
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
    if asset_type is not None and asset_type in DICT_CIS_VULNERABILITY_FLOOD.keys():
        unique_objects_in_asset_type = list(
            DICT_CIS_VULNERABILITY_FLOOD[asset_type].keys()
        )
        features = features[features["object_type"].isin(unique_objects_in_asset_type)]

    # connect maxdam to exposure
    if maxdam_path is None:
        assert "maximum_damage" in features.columns, (
            "If maximum_damage is not provided as argument, maximum damage must be provided in the exposure data."
        )
    else:
        try:
            features["maximum_damage"] = features.apply(
                lambda x: maxdam[x["object_type"]], axis=1
            )
        except KeyError:
            missing_object_types = [
                i for i in features.object_type.unique() if i not in maxdam.keys()
            ]
            raise KeyError(
                f"Not all object types in the exposure are included in the maximum damage file: {missing_object_types}"
            )

    tqdm.pandas(desc="Calculating damage", disable=disable_progress)

    # Calculate damage
    if not multi_curves:
        features = _estimate_damage(features, curves, cell_area_m2)
    else:
        collect_sub_outcomes = []
        for curve_id in multi_curves:
            curves = multi_curves[curve_id]
            collect_sub_outcomes.append(
                _estimate_damage(features, curves, cell_area_m2)["damage"]
            )

        all_curve_damages = pd.concat(collect_sub_outcomes, axis=1)
        all_curve_damages.columns = multi_curves.keys()

        # add all curve damages to the features dataframe
        features.loc[:, all_curve_damages.columns] = all_curve_damages

        if "damage" in features.columns:
            features = features.drop(columns="damage")

    return features
