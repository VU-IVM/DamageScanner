"""Vector-based damage estimation functions for DamageScanner."""

import itertools
import warnings
from collections.abc import Iterator
from pathlib import Path, PurePath

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import rasterio
import rasterio.sample
import shapely
import xarray as xr
from exactextract import exact_extract
from pyproj import Geod
from shapely.geometry import LineString, Point
from tqdm import tqdm

from damagescanner.osm import read_osm_data


def _crs_is_meters(crs: pyproj.CRS) -> bool:
    """Check if a CRS uses meters as its unit.

    Returns:
        True if the CRS uses meters, False otherwise.
    """
    try:
        epsg_code = crs.to_epsg()
        if epsg_code is not None:
            return pyproj.CRS.from_epsg(epsg_code).axis_info[0].unit_name == "metre"
        else:
            # Fallback: check directly from the CRS axis info
            if crs.axis_info:
                return crs.axis_info[0].unit_name == "metre"
            return False
    except Exception:
        return False


def _ensure_list(x: None | float | np.ndarray | list) -> list:
    """
    Ensure a value is a Python list.

    Converts numpy arrays to lists, handles None and NaN values.

    Args:
        x: Input value (can be list, np.ndarray, None, float, or NaN).

    Returns:
        A Python list (empty if input was None/NaN).
    """
    if x is None:
        return []
    if isinstance(x, float) and np.isnan(x):
        return []
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, list):
        return x
    return []


def _convert_to_meters(feature: pd.Series) -> list[float]:
    """
    Convert coverage length to meters for each row with LineString geometries.

    Args:
        feature: A GeoSeries row containing geometry and coverage fields.

    Returns:
        Coverage values in meters.
    """
    line_string = feature.geometry

    # only continue if geometry is a line
    if line_string.geom_type not in ["LineString", "MultiLineString"]:
        return feature.coverage

    if line_string.geom_type == "MultiLineString":
        line_string = line_string.geoms[0]

    geod = Geod(ellps="WGS84")
    coverage_meters = []
    for cover in feature.coverage:
        new_for_length = LineString(
            [Point(line_string.coords[0]), line_string.interpolate(cover)]
        )
        coverage_meters.append(geod.geometry_length(new_for_length))

    return coverage_meters


def _get_cell_area_m2(
    features: gpd.GeoDataFrame, hazard_crs: pyproj.CRS, hazard_resolution: float
) -> float:
    """
    Estimate the area (m²) of a raster grid cell using the feature centroid and resolution.

    Args:
        features: Feature set used to center the buffer.
        hazard_crs: Hazard raster CRS.
        hazard_resolution: Hazard raster resolution (in degrees).

    Returns:
        Grid cell area in square meters.
    """
    geod = Geod(ellps="WGS84")

    asset_point = features.geometry.iloc[0]
    if features.crs is not None and features.crs != hazard_crs:
        asset_point = (
            gpd.GeoSeries([asset_point], crs=features.crs).to_crs(hazard_crs).iloc[0]
        )

    # geodesic area of one cell centred on the asset; east-west extent shrinks
    # with cos(latitude), so squaring a single edge length overstates it (#82)
    cell = asset_point.centroid.buffer(hazard_resolution / 2, cap_style="square")
    return abs(geod.geometry_area_perimeter(cell)[0])


def _tile_slices(n: int, k: int) -> list[slice]:
    """Split ``n`` cells into consecutive index slices of about ``k`` cells.

    A trailing sliver narrower than two cells is merged into the previous
    slice, because rioxarray cannot derive a resolution from one coordinate.

    Args:
        n: Number of cells along the axis.
        k: Target number of cells per slice.

    Returns:
        Consecutive slices that together cover range(n) exactly once.
    """
    k = max(int(k), 2)
    starts = list(range(0, n, k))
    if len(starts) > 1 and n - starts[-1] < 2:
        starts.pop()
    ends = starts[1:] + [n]
    return [slice(a, b) for a, b in zip(starts, ends)]


def _iter_raster_tiles(
    hazard: xr.Dataset | xr.DataArray, tile_size: float
) -> Iterator[xr.Dataset | xr.DataArray]:
    """Yield tiles that partition the raster cells exactly.

    Tiling by array index (rather than by coordinates with ``rio.clip_box``)
    guarantees every cell belongs to exactly one tile, so results from
    neighbouring tiles can be concatenated without double counting.

    Args:
        hazard: Hazard raster.
        tile_size: Tile edge length in the units of the hazard CRS.

    Yields:
        Subsets of the hazard raster, one per tile.
    """
    res_x, res_y = (abs(r) for r in hazard.rio.resolution())
    x_dim, y_dim = hazard.rio.x_dim, hazard.rio.y_dim
    for xs in _tile_slices(hazard.sizes[x_dim], round(tile_size / res_x)):
        for ys in _tile_slices(hazard.sizes[y_dim], round(tile_size / res_y)):
            yield hazard.isel({x_dim: xs, y_dim: ys})


def _geographic_cell_area_m2(
    center_y: np.ndarray, res_x: float, res_y: float
) -> np.ndarray:
    """Geodesic area (m²) of geographic raster cells, given their centre latitudes.

    Args:
        center_y: Latitude of each cell centre, in degrees.
        res_x: Cell width in degrees.
        res_y: Cell height in degrees.

    Returns:
        Area of each cell in m², in the order of ``center_y``.
    """
    geod = Geod(ellps="WGS84")
    center_y = np.asarray(center_y, dtype=float)
    unique_lat, inverse = np.unique(np.round(center_y, 9), return_inverse=True)
    areas = np.empty(len(unique_lat))
    for i, lat in enumerate(unique_lat):
        lats = [lat - res_y / 2, lat - res_y / 2, lat + res_y / 2, lat + res_y / 2]
        lons = [-res_x / 2, res_x / 2, res_x / 2, -res_x / 2]
        areas[i] = abs(geod.polygon_area_perimeter(lons, lats)[0])
    return areas[inverse]


def _create_grid(bbox: shapely.Geometry, height: float) -> np.ndarray:
    """
    Create a regular vector grid over a bounding box.

    Args:
        bbox: The bounding box to cover.
        height: Cell height (assumed square cells).

    Returns:
        List of shapely Polygon grid cells.
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


def _remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge the per tile results of features that span several tiles.

    Tiles partition the raster cells exactly (see ``_iter_raster_tiles``), so
    the cell lists of different tiles never overlap and can be concatenated.

    Args:
        df: DataFrame of list columns, possibly with repeated index entries.

    Returns:
        DataFrame with one row per feature.
    """
    if not df.index.has_duplicates:
        return df
    # only features that span several tiles need merging
    dup = df.index.duplicated(keep=False)
    merged = (
        df[dup]
        .groupby(level=0)
        .agg(
            lambda col: list(
                itertools.chain.from_iterable(_ensure_list(v) for v in col)
            )
        )
    )
    return pd.concat([df[~dup], merged]).sort_index()


def _reproject(
    hazard: xr.Dataset | xr.DataArray,
    features: gpd.GeoDataFrame,
    hazard_crs: pyproj.CRS,
) -> tuple[xr.Dataset | xr.DataArray, gpd.GeoDataFrame, str]:
    """
    Reproject features to match hazard CRS. Also infers a suitable UTM CRS.

    Args:
        hazard: Hazard dataset.
        features: Exposure features.
        hazard_crs: Coordinate system of hazard layer.

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


def _overlay_raster_vector(
    hazard: xr.Dataset | xr.DataArray | rasterio.DatasetReader,
    features: gpd.GeoDataFrame,
    hazard_crs: pyproj.CRS,
    gridded: bool = True,
    disable_progress: bool = False,
    extract_strategy: str = "raster-sequential",
    return_full: bool = True,
) -> gpd.GeoDataFrame:
    """
    Overlay raster hazard data on vector exposure features.

    Args:
        hazard: Raster hazard layer.
        features: Vector exposure features.
        hazard_crs: CRS of hazard data.
        gridded: Whether to process in spatial chunks.
        disable_progress: Disable tqdm progress bar.
        extract_strategy: exactextract strategy - "feature-sequential" or "raster-sequential".
        return_full: Whether to return all features, even those with no hazard intersection.

    Returns:
        Exposure features with added `coverage` and `values`.
    """
    # make sure the hazard data has a crs
    if hazard_crs.to_epsg() is None:
        RuntimeWarning(
            "Hazard crs is not correctly defined. We will now assume it is EPSG:4326"
        )
        hazard = hazard.rio.write_crs("EPSG:4326")
        hazard_crs = pyproj.CRS.from_epsg(4326)

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
            result = hazard.sel(
                {
                    hazard.rio.x_dim: xr.DataArray(features[point_objects].geometry.x),
                    hazard.rio.y_dim: xr.DataArray(features[point_objects].geometry.y),
                },
                method="nearest",
            )

            # Convert to DataArray only if it's a Dataset
            if isinstance(result, xr.Dataset):
                result = result.to_dataarray()

            values = result.values.flatten()

        # Convert values to lists directly (to match polygon/line format)
        values_as_lists = [[v] if v > 0 else [] for v in values]
        coverage_as_lists = [[1] if v > 0 else [] for v in values]

        # Initialize columns as object dtype if they don't exist
        if "values" not in features.columns:
            features["values"] = pd.Series(
                [None] * len(features), index=features.index, dtype=object
            )
        if "coverage" not in features.columns:
            features["coverage"] = pd.Series(
                [None] * len(features), index=features.index, dtype=object
            )

        # Convert to object dtype to hold lists
        features["values"] = features["values"].astype(object)
        features["coverage"] = features["coverage"].astype(object)

        # Assign using a temporary Series with matching index
        point_idx = features.index[point_objects]
        features.loc[point_idx, "values"] = pd.Series(
            values_as_lists, index=point_idx
        ).values
        features.loc[point_idx, "coverage"] = pd.Series(
            coverage_as_lists, index=point_idx
        ).values

    output_type: str = "pandas"
    include_geom: bool = False

    # ask exactextract for each cell's latitude on geographic rasters
    geographic = not _crs_is_meters(hazard_crs)
    extract_ops = ["coverage", "values"] + (["center_y"] if geographic else [])

    if not gridded:
        if area_and_line_objects.sum() > 0:
            values_and_coverage_per_area_and_line_object = exact_extract(
                hazard,
                features[area_and_line_objects][["geometry"]],  # only pass the geometry
                extract_ops,
                include_geom=include_geom,
                output=output_type,
                strategy=extract_strategy,
            )

            # Set the index to match the features (like gridded case does)
            values_and_coverage_per_area_and_line_object.index = features[
                area_and_line_objects
            ].index

            # Initialize columns as object dtype if they don't exist
            if "values" not in features.columns:
                features["values"] = pd.Series(
                    [None] * len(features), index=features.index, dtype=object
                )
            if "coverage" not in features.columns:
                features["coverage"] = pd.Series(
                    [None] * len(features), index=features.index, dtype=object
                )

            features.loc[area_and_line_objects, "coverage"] = (
                values_and_coverage_per_area_and_line_object["coverage"]
            )
            features.loc[area_and_line_objects, "values"] = (
                values_and_coverage_per_area_and_line_object["values"]
            )
            # CHANGED (#82): carry cell latitudes onto the features
            if geographic:
                features["center_y"] = pd.Series(
                    [None] * len(features), index=features.index, dtype=object
                )
                features.loc[area_and_line_objects, "center_y"] = (
                    values_and_coverage_per_area_and_line_object["center_y"]
                )

    elif gridded:
        if area_and_line_objects.sum() > 0:
            if _crs_is_meters(hazard_crs):
                grid_cell_size = 50000  # in meters
            else:
                grid_cell_size = 0.5  # in degrees

            # CHANGED (#83): tiles by array index replace _create_grid + rio.clip_box
            collect_overlay = []
            tiles = list(_iter_raster_tiles(hazard, grid_cell_size))

            for subset_hazard in tqdm(
                tiles,
                total=len(tiles),
                desc="Overlay raster with vector",
                disable=disable_progress,
            ):
                sub_bbox = shapely.box(*subset_hazard.rio.bounds())

                candidate_idx = features.sindex.query(sub_bbox, predicate="intersects")
                subset_features = features.iloc[candidate_idx]
                subset_features = subset_features[
                    subset_features.geom_type.isin(
                        ["Polygon", "MultiPolygon", "LineString", "MultiLineString"]
                    )
                ]
                if len(subset_features) == 0:
                    continue

                values_and_coverage_per_area_and_line_object = exact_extract(
                    subset_hazard,
                    subset_features[["geometry"]],
                    extract_ops,
                    include_geom=include_geom,
                    output=output_type,
                    strategy=extract_strategy,
                )
                values_and_coverage_per_area_and_line_object.index = (
                    subset_features.index
                )
                collect_overlay.append(values_and_coverage_per_area_and_line_object)

            if len(collect_overlay) == 0:
                df: pd.DataFrame = pd.DataFrame(
                    columns=["coverage", "values"], dtype=object
                )
            else:
                df: pd.DataFrame = pd.concat(collect_overlay).sort_index()

            # remove duplicates
            df_no_duplicates = _remove_duplicates(df)

            ## add features to the original features
            features.loc[df_no_duplicates.index, "coverage"] = df_no_duplicates[
                "coverage"
            ]

            features.loc[df_no_duplicates.index, "values"] = df_no_duplicates["values"]
            # CHANGED (#82): carry cell latitudes onto the features
            if geographic:
                features["center_y"] = pd.Series(
                    [None] * len(features), index=features.index, dtype=object
                )
                if "center_y" in df_no_duplicates:
                    features.loc[df_no_duplicates.index, "center_y"] = df_no_duplicates[
                        "center_y"
                    ]

    # Sometimes, with large datasets, a feature may have been excluded from the bbox
    # this has resulted in a null value for the coverage and values. We remove these features.
    if not return_full:
        features = features[~features["values"].isnull()]

        # only keep features with values
        features = features[features["values"].apply(lambda x: len(x) > 0)]

    else:
        # If we return all features, ensure values and coverage are at least empty lists
        features["values"] = features["values"].apply(_ensure_list)
        features["coverage"] = features["coverage"].apply(_ensure_list)

    # convert polygon coverage (fraction of a cell) to m² using the geodesic
    # area of each individual cell, which shrinks with cos(latitude)
    if geographic and "center_y" in features:
        res_x, res_y = (abs(r) for r in hazard.rio.resolution())
        is_polygon = features.geom_type.isin(["Polygon", "MultiPolygon"]).to_numpy()
        poly_idx = features.index[is_polygon]
        cov = [_ensure_list(c) for c in features.loc[poly_idx, "coverage"]]
        lat = [_ensure_list(y) for y in features.loc[poly_idx, "center_y"]]
        lengths = np.fromiter((len(c) for c in cov), dtype=int, count=len(cov))
        if lengths.sum() > 0:
            # flatten all cells, compute each distinct row area once, split back
            flat_cov = np.fromiter(itertools.chain.from_iterable(cov), dtype=float)
            flat_lat = np.fromiter(itertools.chain.from_iterable(lat), dtype=float)
            flat_m2 = flat_cov * _geographic_cell_area_m2(flat_lat, res_x, res_y)
            parts = np.split(flat_m2, np.cumsum(lengths)[:-1])
            features.loc[poly_idx, "coverage"] = pd.Series(
                [p.tolist() for p in parts], index=poly_idx, dtype=object
            )
        features = features.drop(columns="center_y")

    # convert coverage to meters, only do this if the crs is not in meters
    if not _crs_is_meters(hazard_crs):
        tqdm.pandas(desc="convert coverage to meters", disable=disable_progress)

        features.loc[:, "coverage"] = features.progress_apply(
            lambda feature: _convert_to_meters(feature), axis=1
        )

    return features


def _overlay_vector_vector(
    hazard: gpd.GeoDataFrame,
    features: gpd.GeoDataFrame,
    hazard_value_col: str = "band_data",
    gridded: bool = False,
    disable_progress: bool = False,
    return_full: bool = True,
) -> gpd.GeoDataFrame:
    """
    Overlay a vector hazard layer onto vector exposure features.

    Args:
        hazard: Hazard vector features with geometry and value column.
        features: Exposure vector features.
        hazard_value_col: Column name in hazard containing intensity values.
        gridded: Chunk processing toggle (not yet implemented).
        disable_progress: Disable tqdm progress bar.
        return_full: Whether to return all features, even those with no hazard intersection.

    Returns:
        Features with added `coverage` and `values` columns.
    """
    # Store original CRS for later conversion
    original_crs = features.crs

    # Make sure the hazard data has a crs
    if hazard.crs is None or hazard.crs.to_epsg() is None:
        warnings.warn(
            "Hazard CRS is not correctly defined. We will now assume it is EPSG:4326"
        )
        hazard = hazard.set_crs("EPSG:4326")  # ty:ignore[invalid-assignment]

    hazard_crs = hazard.crs

    # Reproject features to match hazard CRS if needed
    if features.crs != hazard_crs:
        features = features.to_crs(hazard_crs)

    # Reproject both to a metric CRS for accurate area/length calculations
    if not _crs_is_meters(hazard_crs):
        # Use a common metric CRS (EPSG:3035 for Europe, or estimate from centroid)
        centroid = features.geometry.union_all().centroid
        utm_zone = int((centroid.x + 180) / 6) + 1
        hemisphere = "north" if centroid.y >= 0 else "south"
        metric_crs = (
            f"EPSG:{32600 + utm_zone if hemisphere == 'north' else 32700 + utm_zone}"
        )

        hazard_metric = hazard.to_crs(metric_crs)
        features_metric = features.to_crs(metric_crs)
    else:
        hazard_metric = hazard
        features_metric = features
        metric_crs = hazard_crs

    # Identify geometry types
    area_objects = features_metric.geom_type.isin(["Polygon", "MultiPolygon"])
    line_objects = features_metric.geom_type.isin(["LineString", "MultiLineString"])
    point_objects = features_metric.geom_type == "Point"

    area_and_line_objects = area_objects | line_objects

    # Initialize columns as object dtype to hold lists
    features["values"] = pd.Series(
        [None] * len(features), index=features.index, dtype=object
    )
    features["coverage"] = pd.Series(
        [None] * len(features), index=features.index, dtype=object
    )

    # Build spatial index for hazard
    hazard_tree = shapely.STRtree(hazard_metric.geometry.values)

    # Handle point objects
    if point_objects.sum() > 0:
        point_indices = features_metric[point_objects].index

        for idx in tqdm(
            point_indices,
            desc="Overlay points with vector hazard",
            disable=disable_progress,
        ):
            point_geom = features_metric.loc[idx, "geometry"]

            # Find hazard features containing this point
            candidate_idx = hazard_tree.query(point_geom, predicate="intersects")

            if len(candidate_idx) > 0:
                # Get the hazard values at this point
                hazard_values = hazard_metric.iloc[candidate_idx][
                    hazard_value_col
                ].values
                # Filter out NaN and zero values
                valid_mask = ~np.isnan(hazard_values) & (hazard_values > 0)
                hazard_values = hazard_values[valid_mask].tolist()

                if len(hazard_values) > 0:
                    features.at[idx, "values"] = hazard_values
                    features.at[idx, "coverage"] = [1.0] * len(hazard_values)
                else:
                    features.at[idx, "values"] = []
                    features.at[idx, "coverage"] = []
            else:
                features.at[idx, "values"] = []
                features.at[idx, "coverage"] = []

    # Handle area and line objects
    if area_and_line_objects.sum() > 0:
        area_line_indices = features_metric[area_and_line_objects].index

        for idx in tqdm(
            area_line_indices,
            desc="Overlay areas/lines with vector hazard",
            disable=disable_progress,
        ):
            feature_geom = features_metric.loc[idx, "geometry"]
            is_line = features_metric.loc[idx, "geometry"].geom_type in [
                "LineString",
                "MultiLineString",
            ]

            # Find candidate hazard features
            candidate_idx = hazard_tree.query(feature_geom, predicate="intersects")

            if len(candidate_idx) == 0:
                features.at[idx, "values"] = []
                features.at[idx, "coverage"] = []
                continue

            # Get intersecting hazard features
            candidate_hazards = hazard_metric.iloc[candidate_idx]

            values_list = []
            coverage_list = []

            for h_idx in candidate_hazards.index:
                h_geom = candidate_hazards.loc[h_idx, "geometry"]
                h_value = candidate_hazards.loc[h_idx, hazard_value_col]

                # Skip NaN or zero hazard values
                if pd.isna(h_value) or h_value <= 0:
                    continue

                # Calculate actual intersection
                try:
                    intersection = shapely.intersection(feature_geom, h_geom)
                except Exception:
                    continue

                if intersection.is_empty:
                    continue

                # Calculate coverage (length for lines, area for polygons) in meters
                if is_line:
                    coverage = shapely.length(intersection)
                else:
                    coverage = shapely.area(intersection)

                if coverage > 0:
                    values_list.append(float(h_value))
                    coverage_list.append(float(coverage))

            features.at[idx, "values"] = values_list
            features.at[idx, "coverage"] = coverage_list

    # Remove features with no values
    if not return_full:
        features = features[~features["values"].isnull()]
        features = features[features["values"].apply(lambda x: len(x) > 0)]
    else:
        # If we return all features, ensure values and coverage are at least empty lists
        features["values"] = features["values"].apply(
            lambda x: x if x is not None else []
        )
        features["coverage"] = features["coverage"].apply(
            lambda x: x if x is not None else []
        )

    # Restore original CRS
    if features.crs != original_crs:
        features = features.to_crs(original_crs)

    return features


def _estimate_damage(
    features: gpd.GeoDataFrame,
    curves: pd.DataFrame,
    object_col: str,
    cell_area_m2: float,
) -> gpd.GeoDataFrame:
    """
    Estimate total damage per asset using vulnerability curves.

    Args:
        features: Exposure with hazard info.
        curves: Vulnerability curve per asset type.
        object_col: Name of the object type column.
        cell_area_m2: Area per grid cell in m² (for polygons).

    Features whose object type has no column in ``curves`` are not assessed
    with this curve set and get NaN damage, which is distinct from a damage
    of 0 for an asset that is covered but undamaged (#58).

    Returns:
        Exposure data with added `damage` column.
    """
    applicable = features[object_col].isin(curves.columns)
    damage = pd.Series(np.nan, index=features.index, dtype=float)
    if applicable.any():
        damage[applicable] = features[applicable].progress_apply(
            lambda _object: _get_damage_per_object(
                _object, curves, object_col, cell_area_m2
            ),
            axis=1,
        )
    features["damage"] = damage

    return features


def _get_damage_per_object(
    asset: pd.Series, curves: pd.DataFrame, object_col: str, cell_area_m2: float
) -> float:
    """
    Compute damage for a single asset using interpolated vulnerability.

    Args:
        asset: Single feature with geometry, values, object_type.
        curves: Vulnerability curves.
        object_col: Name of the object type column in asset and curves.
        cell_area_m2: Cell area in square meters.

    Returns:
        Estimated damage value.

    Raises:
        ValueError: If geometry type is unsupported or if object type is not in curves.
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
            np.interp(asset["values"], curves.index, curves[asset[object_col]].values)
            * coverage
        )
        * asset["maximum_damage"]
    )


def VectorExposure(
    hazard_file: Path
    | xr.Dataset
    | xr.DataArray
    | rasterio.DatasetReader
    | gpd.GeoDataFrame,
    feature_file: Path | gpd.GeoDataFrame | pd.DataFrame | str,
    asset_type: str | None = "roads",
    object_col: str = "object_type",
    hazard_value_col: str = "band_data",
    disable_progress: bool = False,
    gridded: bool = True,
    extract_strategy: str = "raster-sequential",
    return_full: bool = True,
) -> tuple[gpd.GeoDataFrame, str, pyproj.CRS | None, float | None]:
    """
    Load and overlay vector or raster hazard with vector exposure data.

    Args:
        hazard_file: Hazard input.
        feature_file: Exposure input.
        asset_type: Infrastructure category (only for OSM).
        object_col: Name of the object type column.
        hazard_value_col: Column name in vector hazard containing intensity values.
        disable_progress: Whether to suppress progress bars.
        gridded: Whether to process in spatial chunks.
        extract_strategy: exactextract strategy - "feature-sequential" or "raster-sequential".
        return_full: Whether to return all features, even those with no hazard intersection.

    Returns:
        tuple: (features, object_col, hazard_crs, cell_area_m2)

    Raises:
        ValueError: If input files are not in expected formats or if geometry types are unsupported.
    """
    # load exposure data
    if isinstance(feature_file, Path):
        # if exposure_file is a shapefile, geopackage or parquet file
        if feature_file.suffix in [".shp", ".gpkg"]:
            features = gpd.read_file(feature_file)
        elif feature_file.suffix == ".parquet":
            features = gpd.read_parquet(feature_file)
        # if exposure_file is an osm.pbf file
        elif feature_file.suffix == ".pbf":
            if asset_type is None:
                raise ValueError(
                    "When using an OSM file as exposure, you need to specify the asset type (e.g. 'roads' or 'buildings')"
                )
            features = read_osm_data(feature_file, asset_type)
        else:
            raise ValueError(
                "exposure data should either be a shapefile, geopackage, parquet or osm.pbf file"
            )

    elif isinstance(feature_file, (gpd.GeoDataFrame, pd.DataFrame)):
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
            if hazard_file.suffix in (".tif", ".tiff"):
                assert hazard.band.size == 1, (
                    "Hazard data should only contain one band. If you have multiple bands, please select one band using the `band` argument."
                )
                hazard = hazard["band_data"].sel(band=1)
            hazard_crs = hazard.rio.crs

            # check if crs is already in meters
            if _crs_is_meters(hazard_crs):
                cell_area_m2: int | float = abs(
                    (hazard.x[1].values - hazard.x[0].values)
                    * (hazard.y[0].values - hazard.y[1].values)
                )
            else:
                # CHANGED (#82): polygon coverage is converted to m² per cell during the overlay
                cell_area_m2 = 1

        elif hazard_file.suffix in [".shp", ".gpkg"]:
            hazard = gpd.read_file(hazard_file)
            hazard_crs = hazard.crs
            cell_area_m2 = 1  # For vector hazards, coverage is already in meters
        elif hazard_file.suffix == ".parquet":
            hazard = gpd.read_parquet(hazard_file)
            hazard_crs = hazard.crs
            cell_area_m2 = 1  # For vector hazards, coverage is already in meters
        else:
            raise ValueError(
                "hazard data should either be a geotiff, netcdf, shapefile, geopackage or parquet file"
            )
    elif isinstance(hazard_file, rasterio.io.DatasetReader):
        hazard = hazard_file.copy()
        hazard_crs = hazard.crs
        cell_area_m2 = _get_cell_area_m2(
            features, hazard_crs, abs(hazard.rio.resolution()[0])
        )
    elif isinstance(hazard_file, (xr.Dataset, xr.DataArray)):
        hazard = hazard_file.copy()
        hazard_crs = hazard.rio.crs

        # check if crs is already in meters
        if _crs_is_meters(hazard_crs):
            cell_area_m2: float | int = abs(
                (hazard.x[1].values - hazard.x[0].values)
                * (hazard.y[0].values - hazard.y[1].values)
            )
        # if not, extract it more cumbersome
        else:
            # CHANGED (#82): polygon coverage is converted to m² per cell during the overlay
            cell_area_m2 = 1

    elif isinstance(hazard_file, gpd.GeoDataFrame):
        hazard = hazard_file.copy()
        hazard_crs = hazard.crs
        cell_area_m2 = 1  # For vector hazards, coverage is already in meters
    else:
        raise ValueError(
            f"Hazard should be a raster or GeoDataFrame object, {type(hazard_file)} given"
        )

    # Run exposure overlay
    if isinstance(hazard, (rasterio.io.DatasetReader, xr.Dataset, xr.DataArray)):
        features = _overlay_raster_vector(
            hazard,
            features,
            hazard_crs,
            disable_progress=disable_progress,
            gridded=gridded,
            extract_strategy=extract_strategy,
            return_full=return_full,
        )

    elif isinstance(hazard, gpd.GeoDataFrame):
        features = _overlay_vector_vector(
            hazard,
            features,
            hazard_value_col=hazard_value_col,
            gridded=gridded,
            disable_progress=disable_progress,
            return_full=return_full,
        )

    return features, object_col, hazard_crs, cell_area_m2


def VectorScanner(
    hazard_file: Path | xr.Dataset | xr.DataArray | gpd.GeoDataFrame | str,
    feature_file: Path | gpd.GeoDataFrame | str,
    curve_path: Path | pd.DataFrame | str,
    maxdam_path: Path | pd.DataFrame | dict | None = None,
    asset_type: str | None = None,
    multi_curves: dict = dict(),
    object_col: str = "object_type",
    hazard_value_col: str = "band_data",
    disable_progress: bool = False,
    gridded: bool = True,
    extract_strategy: str = "raster-sequential",
    return_full: bool = True,
) -> gpd.GeoDataFrame:
    """
    Perform vector-based direct damage assessment using hazard and exposure layers.

    Args:
        hazard_file: Hazard input.
        feature_file: Exposure input.
        curve_path: Vulnerability curve(s).
        maxdam_path: Maximum damage values.
        asset_type: Infrastructure class (only for OSM).
        multi_curves: Multiple curve sets, as {curve_id: DataFrame}. Each DataFrame
            only needs columns for the object types that curve set applies to;
            other object types get NaN in that curve set's damage column.
        object_col: Column name with object type.
        hazard_value_col: Column name in vector hazard containing intensity values.
        disable_progress: Whether to suppress progress bars.
        gridded: Whether to process in spatial chunks.
        extract_strategy: exactextract strategy - "feature-sequential" or "raster-sequential".
        return_full: Whether to return all features, even those with no hazard intersection.

    Returns:
        Exposure data with calculated damages.

    Raises:
        ValueError: If input files are of unsupported formats or if geometry types are not supported.
        KeyError: If object types in exposure are not covered by maximum damage data.
    """
    if isinstance(feature_file, str):
        feature_file = Path(feature_file)
    if isinstance(hazard_file, str):
        hazard_file = Path(hazard_file)
    if isinstance(curve_path, str):
        curve_path = Path(curve_path)
    if isinstance(maxdam_path, str):
        maxdam_path = Path(maxdam_path)

    # Load hazard and exposure data, and perform the overlay
    features: gpd.GeoDataFrame
    features, object_col, hazard_crs, cell_area_m2 = VectorExposure(
        hazard_file=hazard_file,
        feature_file=feature_file,
        asset_type=asset_type,
        object_col=object_col,
        hazard_value_col=hazard_value_col,
        disable_progress=disable_progress,
        gridded=gridded,
        extract_strategy=extract_strategy,
        return_full=return_full,
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
        maxdam = dict(zip(maxdam[object_col], maxdam["damage"]))
    elif isinstance(maxdam_path, pd.DataFrame):
        maxdam = dict(zip(maxdam_path[object_col], maxdam_path["damage"]))
    elif isinstance(maxdam_path, np.ndarray):
        if not isinstance(maxdam_path, np.ndarray) or maxdam_path.shape[1] != 2:
            raise ValueError(
                "If maximum damage is provided as a numpy array, it should be a 2D array with two columns: object type and damage value."
            )
        maxdam = dict(zip(maxdam_path[:, 0], maxdam_path[:, 1]))  # ty:ignore[invalid-argument-type]
    elif isinstance(maxdam_path, dict):
        maxdam = maxdam_path

    # remove features that are not part of this object type
    # Only filter by asset_type when using OSM data
    # is_osm_data = isinstance(feature_file, PurePath) and feature_file.suffix == ".pbf"

    # if is_osm_data and asset_type is not None and asset_type in DICT_CIS_VULNERABILITY_FLOOD.keys():
    #     unique_objects_in_asset_type = list(
    #         DICT_CIS_VULNERABILITY_FLOOD[asset_type].keys()
    #     )
    #     features = features[features[object_col].isin(unique_objects_in_asset_type)]

    # connect maxdam to exposure
    if maxdam_path is None:
        assert "maximum_damage" in features.columns, (
            "If maximum_damage is not provided as argument, maximum damage must be provided in the exposure data."
        )
    else:
        try:
            features["maximum_damage"] = features.apply(
                lambda x: maxdam[x[object_col]], axis=1
            )
        except KeyError:
            missing_object_types = [
                i for i in features[object_col].unique() if i not in maxdam.keys()
            ]
            raise KeyError(
                f"Not all object types in the exposure are included in the maximum damage file: {missing_object_types}"
            )

    tqdm.pandas(desc="Calculating damage", disable=disable_progress)

    # Check curve coverage of the object types in the exposure
    object_types = set(features[object_col].dropna().unique())
    if not multi_curves:
        missing = sorted(object_types - set(curves.columns), key=str)
        if missing:
            raise KeyError(
                f"Not all object types in the exposure have a vulnerability curve: {missing}"
            )
    else:
        # a curve set only needs columns for the object types it applies to (#58)
        covered = set().union(*(set(c.columns) for c in multi_curves.values()))
        uncovered = sorted(object_types - covered, key=str)
        if uncovered:
            warnings.warn(
                f"These object types have no curve in any curve set and get NaN damage: {uncovered}"
            )

    # Calculate damage
    if not multi_curves:
        features: gpd.GeoDataFrame = _estimate_damage(
            features, curves, object_col, cell_area_m2
        )
    else:
        collect_sub_outcomes = []
        for curve_id in multi_curves:
            curves = multi_curves[curve_id]
            collect_sub_outcomes.append(
                _estimate_damage(features, curves, object_col, cell_area_m2)["damage"]
            )

        all_curve_damages = pd.concat(collect_sub_outcomes, axis=1)
        all_curve_damages.columns = multi_curves.keys()

        # add all curve damages to the features dataframe
        features.loc[:, all_curve_damages.columns] = all_curve_damages

        if "damage" in features.columns:
            features: gpd.GeoDataFrame = features.drop(columns="damage")  # ty: ignore[invalid-assignment]

    return features
