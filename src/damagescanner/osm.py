import re
import functools
import operator
import numpy as np
import shapely
from shapely.geometry import (
    MultiPolygon,
    GeometryCollection,
)
import pandas as pd
import geopandas as gpd

from damagescanner.config import DICT_CIS_VULNERABILITY_FLOOD

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
        "osm_keys": ["aeroway", "name", ""],
        "osm_query": {"aeroway": ["aerodrome", "apron", "terminal", "runway"]},
    },
    "telecom": {
        "osm_keys": ["man_made", "tower:type", "name"],
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
        "osm_keys": ["power", "voltage", "utility", "name", "source"],
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
    Combine two values into a single object type string.

    Args:
        a (str or None): First attribute value.
        b (str or None): Second attribute value.

    Returns:
        str or None: Combined string or fallback based on logic.
    """
    if pd.notna(a) and not pd.notna(b):  # if only a contains a string
        return f"{a}"
    elif pd.notna(b) and not pd.notna(a):  # if only b contains a string
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
    Combine values from specified columns to create the `object_type` field.

    Args:
        features (gpd.GeoDataFrame): Input GeoDataFrame with OSM attributes.
        column_names_lst (list of str): Columns to merge into `object_type` (2 or 3 max).

    Returns:
        gpd.GeoDataFrame: DataFrame with filtered and renamed columns.
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
    Remove assets whose geometries are fully contained within others.

    Args:
        features (gpd.GeoDataFrame): GeoDataFrame with point and polygon features.

    Returns:
        gpd.GeoDataFrame: Cleaned GeoDataFrame with unique geometries.
    """
    features = _remove_contained_polys(
        _remove_contained_points(features)
    )  # remove points and polygons within a (larger) polygon

    return features


def extract_first_geom(geom):
    """
    Extract the first geometry from a GeometryCollection.

    Args:
        geom (shapely.Geometry): Shapely geometry object.

    Returns:
        shapely.Geometry: First geometry or unchanged object.
    """
    if isinstance(geom, GeometryCollection) and len(geom.geoms) > 0:
        return geom.geoms[0]

    return geom


def _remove_contained_points(gdf_p_mp):
    """
    Remove point features contained within any polygon in the dataset.

    Args:
        gdf_p_mp (gpd.GeoDataFrame): GeoDataFrame with point and polygon geometries.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame without contained points.
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
    From a GeoDataFrame containing (multi-)polygons (and potentially other
    geometries), remove those polygon entries that are already fully
    contained in another polygon entries. Removes smaller polygons within
    polygons and full duplicates, but leaves contained points untouched
    (see remove_contained_points() for this).

    Resets the index of the dataframe.

    Args:
        gdf (gpd.GeoDataFrame): GeoDataFrame with polygon geometries.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with outermost geometries.
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
    Convert multipolygon geometries into their centroid points.

    Args:
        gdf (gpd.GeoDataFrame): GeoDataFrame with polygon geometries.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with point geometries instead.
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
    """
    Parse the value of a specific key from a semi-structured OSM tag string.

    Args:
        text (str): Raw OSM `other_tags` string.
        key (str): Key to extract value for.

    Returns:
        str or None: Extracted value or None.
    """
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
    Extract specific infrastructure features from a .pbf file using OSM keys/values.

    Args:
        osm_path (str or Path): Path to .osm.pbf file.
        geom_type (str): One of 'points', 'lines', 'multipolygons'.
        osm_keys (list): Keys to extract from OSM file.
        osm_query (dict): Key-value mapping used to filter.

    Returns:
        gpd.GeoDataFrame: Extracted GeoDataFrame with `object_type` field.
    """
    features = gpd.read_file(osm_path, layer=geom_type, engine="pyogrio")

    if "osm_way_id" in features.columns:
        features["osm_id"] = features["osm_id"].fillna(features["osm_way_id"])

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
    Load and extract OSM features for a given critical infrastructure type.

    Args:
        osm_path (str or Path): Path to .osm.pbf file.
        asset_type (str): One of the keys in DICT_CIS_OSM.

    Returns:
        gpd.GeoDataFrame: Cleaned and validated exposure GeoDataFrame.

    Raises:
        ImportWarning: If asset_type is not supported.
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
        raise ImportWarning("feature not in DICT_CIS_OSM. Returning empty gdf")

    # make all geometries valid
    gdf["geometry"] = shapely.make_valid(gdf["geometry"])
    gdf = gdf[gdf.geometry.is_valid]

    # only keep assets with unique geometries
    features = _remove_contained_assets(gdf)

    # remove potential geometrycollections to avoid errors later on
    features["geometry"] = features["geometry"].apply(extract_first_geom)

    # remove features that are not in the asset_type list
    unique_objects_in_asset_type = list(DICT_CIS_VULNERABILITY_FLOOD[asset_type].keys())

    return features[features["object_type"].isin(unique_objects_in_asset_type)]
