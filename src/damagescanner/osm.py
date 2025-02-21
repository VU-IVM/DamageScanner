import os
import re
import functools
import operator
import numpy as np
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
import pandas as pd
import geopandas as gpd

from damagescanner.base_values import DICT_CIS_VULNERABILITY_FLOOD

DICT_CIS_OSM = {
    "roads": {
        "osm_keys": ["highway", "name", "maxspeed", "lanes", "surface"],
        "osm_query": {
            "highway": [
                "motorway",
                "motorway_link",
                "trunk",
                "trunk_link",
                "primary",
                "primary_link",
                "secondary",
                "secondary_link",
                "tertiary",
                "tertiary_link",
                "residential",
                "road",
                "unclassified",
                "track",
            ]
        },
    },
    "main_roads": {
        "osm_keys": ["highway", "name", "maxspeed", "lanes", "surface"],
        "osm_query": {
            "highway": [
                "primary",
                "primary_link",
                "secondary",
                "secondary_link",
                "tertiary",
                "tertiary_link",
                "trunk",
                "trunk_link",
                "motorway",
                "motorway_link",
            ]
        },
    },
    "rail": {
        "osm_keys": ["railway", "name", "gauge", "electrified", "voltage"],
        "osm_query": {"railway": ["rail", "narrow_gauge"]},
    },
    "air": {
        "osm_keys": ["aeroway", "name",""],
        "osm_query": {"aeroway": ["aerodrome","apron", "terminal", "runway"]},
    },
    "telecom": {
        "osm_keys": ["man_made", "tower_type", "name"],
        "osm_query": {
            "man_made": ["mast", "communications_tower"],
            "tower_type": ["communication"],
        },
    },
    "water_supply": {
        "osm_keys": ["man_made", "name"],
        "osm_query": {
            "man_made": [
                "water_works",
                "water_well",
                "water_tower",
                "reservoir_covered",
                "storage_tank",
            ]
        },
    },
    "waste_solid": {
        "osm_keys": ["amenity", "name"],
        "osm_query": {"amenity": ["waste_transfer_station"]},
    },
    "waste_water": {
        "osm_keys": ["man_made", "name"],
        "osm_query": {"man_made": ["wastewater_plant"]},
    },
    "education": {
        "osm_keys": ["amenity", "building", "name"],
        "osm_query": {
            "building": ["school", "kindergarten", "college", "university", "library"],
            "amenity": ["school", "kindergarten", "college", "university", "library"],
        },
    },
    "healthcare": {
        "osm_keys": ["amenity", "building", "healthcare", "name"],
        "osm_query": {
            "amenity": ["hospital", "clinic", "doctors", "dentist", "pharmacy"],
            "building": ["hospital", "clinic"],
            "healthcare": [
                "pharmacy",
                "dentist",
                "physiotherapist",
                "alternative",
                "laboratory",
                "optometrist",
                "rehabilitation",
                "blood_donation",
                "birthing_center",
            ],
        },
    },
    "power": {
        "osm_keys": ["power", "voltage", "utility", "name"],
        "osm_query": {
            "power": [
                "line",
                "cable",
                "minor_line",
                "plant",
                "generator",
                "substation",
                "transformer",
                "pole",
                "portal",
                "tower",
                "terminal",
                "switch",
                "catenary_mast",
            ]
        },
    },
    "gas": {
        "osm_keys": ["man_made", "pipeline", "utility", "name", "substance", "content"],
        "osm_query": {
            "man_made": ["pipeline", "storage_tank"],
            "pipeline": ["substation"],
            "utility": ["gas"],
            "substance": ["gas"],
            "content": ["gas"],
        },
    },
    "food": {
        "osm_keys": ["amenity", "building", "name"],
        "osm_query": {
            "amenity": ["restaurant", "fast_food", "cafe", "pub", "bar"],
            "building": ["restaurant", "fast_food", "cafe", "pub", "bar"],
        },
    },
    "oil": {
        "osm_keys": ["pipeline", "man_made", "amenity", "name", "substance"],
        "osm_query": {
            "pipeline": ["substation"],
            "man_made": ["pipeline", "petroleum_well", "oil_refinery"],
            "amenity": ["fuel"],
            "substance": ["oil"],
        },
    },
    "buildings": {
        "osm_keys": ["building", "amenity", "name"],
        "osm_query": {
            "building": [
                "yes",
                "house",
                "residential",
                "detached",
                "hut",
                "industrial",
                "shed",
                "apartments",
            ]
        },
    },
}


def _combine_columns(a, b):
    """
    Combine values from two input arguments 'a' and 'b' into a single string.
    Arguments:
    - a (str or None): Value from column 'A'.
    - b (str or None): Value from column 'B'.
    Returns:
    - str or None: A string of 'a', 'b' or combination. If both 'a' and 'b' are None, return None.
    """

    if pd.notna(a) and pd.notna(b) == False:  # if only a contains a string
        return f"{a}"
    elif pd.notna(b) and pd.notna(a) == False:  # if only b contains a string
        return f"{b}"
    elif pd.notna(a) and pd.notna(b):  # if both values contain a string
        if a == b:
            return f"{a}"
        elif a == "yes" or b == "yes":
            if a == "yes":
                return f"{b}"
            elif b == "yes":
                return f"{a}"
        else:
            return f"{a}"  # f"{a}_{b}" # assuming that value from column A contains the more detailed information
    else:
        None  # Decision point: If nones are existent, decide on what to do with Nones. Are we sure that these are education facilities? Delete them? Provide another tag to them?


def _filter_dataframe(features, column_names_lst):
    """
    Filter a GeoDataFrame by combining information from two specified columns and removing selected columns.
    Args:
        assets (geopandas.GeoDataFrame): The input GeoDataFrame containing spatial geometries and columns to filter.
        column_names_lst (list): A list of two column names whose information needs to be combined to create a new 'asset' column.
    Returns:
        geopandas.GeoDataFrame: A filtered GeoDataFrame with a new 'asset' column and selected columns dropped, and points converted to polygons.
    """
    if len(column_names_lst) == 2:
        features["object_type"] = features.apply(
            lambda row: _combine_columns(
                row[column_names_lst[0]], row[column_names_lst[1]]
            ),
            axis=1,
        )  # create new column based on tag information provided in two columns
    elif len(column_names_lst) == 3:
        features["object_type_temp"] = features.apply(
            lambda row: _combine_columns(
                row[column_names_lst[0]], row[column_names_lst[1]]
            ),
            axis=1,
        )  # create temp column based on tag information provided in two columns
        features["object_type"] = features.apply(
            lambda row: _combine_columns(
                row["object_type_temp"], row[column_names_lst[2]]
            ),
            axis=1,
        )  # create new column based on tag information provided in two columns
        column_names_lst.append("object_type_temp")
    else:
        print("Warning: column_names_lst should contain 2 or 3 items")
    features = features.drop(columns=column_names_lst, axis=1)  # drop columns

    return features


def _remove_contained_assets(features):
    """
    Process the geometry of assets, removing contained points and polygons, and converting points to polygons.
    Args:
        assets (geopandas.GeoDataFrame): Input GeoDataFrame containing asset geometries.
    Returns:
        geopandas.GeoDataFrame: Processed GeoDataFrame with updated asset geometries.
    """
    features = _remove_contained_polys(
        _remove_contained_points(features)
    )  # remove points and polygons within a (larger) polygon

    return features


def _remove_contained_points(gdf_p_mp):
    """
    from a GeoDataFrame containing points and (multi-)polygons, remove those
    points that are contained in a multipolygons entry.
    Resets the index of the dataframe.

    Parameters
    ----------
    gdf_p_mp : gpd.GeoDataFrame
        GeoDataFrame containing entries with point and (multi-)polygon geometry
    """

    gdf_p_mp = gdf_p_mp.reset_index(drop=True)

    ind_dupl = np.unique(
        gpd.sjoin(
            gdf_p_mp[gdf_p_mp.geometry.type == "Point"],
            gdf_p_mp[
                (gdf_p_mp.geometry.type == "MultiPolygon")
                | (gdf_p_mp.geometry.type == "Polygon")
            ],
            predicate="within",
        ).index
    )

    return gdf_p_mp.drop(index=ind_dupl).reset_index(drop=True)


def _remove_contained_polys(gdf):
    """
    from a GeoDataFrame containing (multi-)polygons (and potentially other
    geometries), remove those polygon entries that are already fully
    contained in another polygon entries. Removes smaller polygons within
    polygons and full duplicates, but leaves contained points untouched
    (see remove_contained_points() for this).

    Resets the index of the dataframe.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame containing entries with (multi-)polygon geometry
    """

    gdf = gdf.reset_index(drop=True)

    contained = gpd.sjoin(
        gdf[(gdf.geometry.type == "MultiPolygon") | (gdf.geometry.type == "Polygon")],
        gdf[(gdf.geometry.type == "MultiPolygon") | (gdf.geometry.type == "Polygon")],
        predicate="contains",
    )

    subset = contained[contained.index != contained.index_right]
    to_drop = set(subset.index_right) - set(subset.index)

    return gdf.drop(index=to_drop).reset_index(drop=True)


def create_point_from_polygon(gdf):
    """
    Transforms polygons into points
    Arguments:
        gdf: A geodataframe containing a column geometry
    Returns:
    - geopandas.GeoDataFrame: The updated GeoDataFrame without polygons but with only point geometries
    """
    gdf["geometry"] = gdf["geometry"].apply(
        lambda geom: MultiPolygon([geom]) if geom.geom_type == "Polygon" else geom
    )  # convert to multipolygons in case polygons are in the df
    # gdf.loc[gdf.geom_type == 'MultiPolygon','geometry'] = gdf.loc[assets.geom_type == 'MultiPolygon'].centroid #convert polygon to point
    gdf.loc[gdf.geom_type == "MultiPolygon", "geometry"] = gdf.loc[
        gdf.geom_type == "MultiPolygon"
    ].centroid  # convert polygon to point
    return gdf


def _extract_value(text, key):
    pattern = rf'"{key}"=>"([^"]+)"'
    try:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
        return None
    except:
        return None


def extract(osm_path, geom_type, osm_keys, osm_query):
    """
    Extracts data from the given osm.pbf file based on the given osm_keys and
    osm_query. The function is used by extract_osm_data() to extract critical
    infrastructure data from the osm.pbf file.
    Parameters
    ----------
    osm_path : str or Path
        location of osm.pbf file from which to parse
    geom_type : str
        one of 'points', 'lines', 'multipolygons'
    osm_keys : list
        list of osm keys to query for
    osm_query : str
        query to be executed on the osm.pbf file
    Returns
    -------
    gdf : geopandas.GeoDataFrame
        geodataframe containing the extracted data
    """
    features = gpd.read_file(osm_path, layer=geom_type, engine="pyogrio")

    if 'osm_way_id' in features.columns:
        features['osm_id'] = features['osm_id'].fillna(features['osm_way_id'])

    for key in osm_keys:
        if key not in features.columns:
            features[key] = features["other_tags"].apply(
                lambda x: _extract_value(x, key)
            )

    # build query
    collect_indices = []
    for query_key in osm_query.keys():
        collect_indices.append(
            features[features[query_key].isin(osm_query[query_key])].index.values
        )

    # get complete list
    collect_indices = functools.reduce(operator.iconcat, collect_indices, [])

    # remove duplicates from list
    collect_indices = list(set(collect_indices))
    features = features.iloc[collect_indices]

    features = features[["osm_id", "geometry"] + osm_keys]

    # make sure we keep collect the necessary feature charactericts in same column
    # and cover the points to polygons for better damage assessment

    if len(osm_keys) < 3:
        features = _filter_dataframe(features, osm_keys)
    else:
        features = _filter_dataframe(features, osm_keys[:3])

    features.rename(columns={osm_keys[0]: "object_type"}, inplace=True)

    return features


def read_osm_data(osm_path, asset_type):
    """
    A wrapper around extract() to conveniently extract map info for a
    selection of  critical infrastructure types from the given osm.pbf file.
    No need to search for osm key/value tags and relevant geometry types.
    Parameters
    ----------
    osm_path : str or Path
        location of osm.pbf file from which to parse
    asset_type : str
        one of DICT_CIS_OSM.keys(), i.e. 'education', 'healthcare',
        'water', 'telecom', 'road', 'rail', 'air', 'gas', 'oil', 'power',
        'wastewater', 'food'
    See also
    -------
    DICT_CIS_OSM for the keys and key/value tags queried for the respective
    CIs. Modify if desired.
    """

    # features consisting in points and multipolygon results:
    if asset_type in ["healthcare", "education", "food", "buildings"]:
        gdf = pd.concat(
            [
                extract(
                    osm_path,
                    "points",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                extract(
                    osm_path,
                    "multipolygons",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
            ]
        )

    # features consisting in points, multipolygons and lines:
    elif asset_type in ["gas", "oil", "water", "power"]:
        gdf = pd.concat(
            [
                extract(
                    osm_path,
                    "points",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                extract(
                    osm_path,
                    "multipolygons",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                extract(
                    osm_path,
                    "lines",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
            ]
        )

    # features consisting in multipolygons and lines:
    elif asset_type in ["air"]:
        gdf = pd.concat(
            [
                extract(
                    osm_path,
                    "multipolygons",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                extract(
                    osm_path,
                    "lines",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
            ]
        )

    # features consisting in multiple datattypes, but only lines needed:
    elif asset_type in ["rail", "roads", "main_roads"]:
        gdf = pd.concat(
            [
                extract(
                    osm_path,
                    "lines",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                )
            ]
        )

    # features consisting in all data types, but only points and multipolygon needed:
    elif asset_type in [
        "telecom",
        "wastewater",
        "waste_solid",
        "waste_water",
        "water_supply",
    ]:
        gdf = pd.concat(
            [
                extract(
                    osm_path,
                    "points",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                extract(
                    osm_path,
                    "multipolygons",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
            ]
        )

    else:
        return ImportWarning("feature not in DICT_CIS_OSM. Returning empty gdf")

    features = _remove_contained_assets(gdf)

    # remove features that are not in the asset_type list
    unique_objects_in_asset_type = list(DICT_CIS_VULNERABILITY_FLOOD[asset_type].keys())

    return features[features["object_type"].isin(unique_objects_in_asset_type)]
