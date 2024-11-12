import os
import re
import functools
import operator
import pandas as pd
import geopandas as gpd

DICT_CIS_OSM = {
    "roads": {
        "osm_keys": ["highway", "name", "maxspeed", "lanes", "surface"],
        "osm_query": {'highway': [
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
        ] }
    },
    "main_roads": {
        "osm_keys": ["highway", "name", "maxspeed", "lanes", "surface"],
        "osm_query": {'highway': [
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
        ] }
    },
    "rail": {
        "osm_keys": ["railway", "name", "gauge", "electrified", "voltage"],
        "osm_query": {"railway": ["rail", "narrow_gauge"]}
    },
    "air": {
        "osm_keys": ["aeroway", "name"],
        "osm_query": {"aeroway":["aerodrome", "terminal", "runway"]}
    },
    "telecom": {
        "osm_keys": ["man_made", "tower_type", "name"],
        "osm_query": {"man_made": ["mast", "communications_tower"],
                      "tower_type": ["communication"]}
    },
    "water_supply": {
        "osm_keys": ["man_made", "name"],
        "osm_query": {"man_made":[
            "water_works",
            "water_well",
            "water_tower",
            "reservoir_covered",
            "storage_tank",
        ]}
    },
    "waste_solid": {
        "osm_keys": ["amenity", "name"],
        "osm_query": {"amenity": ["waste_transfer_station"]}
    },
    "waste_water": {
        "osm_keys": ["man_made", "name"],
        "osm_query": {"man_made":["wastewater_plant"]}
    },
    "education": {
        "osm_keys": ["amenity", "building", "name"],
        "osm_query" : {'building' : ['school','kindergarten',
                                     'college','university',
                                     'library'],
                        'amenity' : ['school','kindergarten',
                                     'college','university',
                                     'library']}
        },
    "healthcare": {
        "osm_keys": ["amenity", "building", "healthcare", "name"],
        "osm_query": {"amenity": ["hospital", "clinic", "doctors", "dentist", "pharmacy"],
                        "building": ["hospital", "clinic"],
                        "healthcare": ["pharmacy", "dentist", "physiotherapist", 
                        "alternative", "laboratory", "optometrist", 
                        "rehabilitation", "blood_donation", "birthing_center"]}
},
    "power": {
        "osm_keys": ["power", "voltage", "utility", "name"],
        "osm_query": {"power" :[
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
        ]}
    },
    "gas": {
        "osm_keys": ["man_made", "pipeline", "utility", "name","substance","content"],
        "osm_query" : {"man_made": ["pipeline", "storage_tank"],
                        "pipeline": ["substation"],
                        "utility": ["gas"],
                        "substance": ["gas"],
                        "content": ["gas"]}
    },
    "food": {
        "osm_keys": ["amenity", "building", "name"],
        "osm_query": {"amenity": ["restaurant", "fast_food", "cafe", "pub", "bar"],
                        "building": ["restaurant", "fast_food", "cafe", "pub", "bar"]}
    },
    "oil": {
        "osm_keys": ["pipeline", "man_made", "amenity", "name","substance"],
        "osm_query": {"pipeline": ["substation"],
                      "man_made": ["pipeline", "petroleum_well", "oil_refinery"],
                        "amenity": ["fuel"],
                        "substance": ["oil"]}},
    "wastewater": {
        "osm_keys": ["man_made", "amenity", "name"],
        "osm_query": {"man_made": ["wastewater_plant"],
                        "amenity": ["waste_transfer_station"]}
    },
    "buildings": {
        "osm_keys": ["building", "amenity", "name"],
        "osm_query": {"building": [
            "yes",
            "house",
            "residential",
            "detached",
            "hut",
            "industrial",
            "shed",
            "apartments",
        ]}
    },
}

def extract_value(text, key):
    pattern = rf'"{key}"=>"([^"]+)"'
    try:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
        return None
    except:
        return None


def _extract(osm_path, geom_type, osm_keys, osm_query):
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
    
    for key in osm_keys:
        if key not in features.columns:
            features[key] = features['other_tags'].apply(lambda x: extract_value(x, key))
    
    # build query
    collect_indices = []
    for query_key in osm_query.keys():
        collect_indices.append(
            features[features[query_key].isin(
                osm_query[query_key])].index.values)
    
    # get complete list
    collect_indices = functools.reduce(
        operator.iconcat, collect_indices, [])

    # remove duplicates from list
    collect_indices = list(set(collect_indices))

    features = features.iloc[collect_indices]

    features = features[["osm_id", "geometry"] + osm_keys]

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
                _extract(
                    osm_path,
                    "points",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                _extract(
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
                _extract(
                    osm_path,
                    "points",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                _extract(
                    osm_path,
                    "multipolygons",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                _extract(
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
                _extract(
                    osm_path,
                    "multipolygons",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                _extract(
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
                _extract(
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
                _extract(
                    osm_path,
                    "points",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
                _extract(
                    osm_path,
                    "multipolygons",
                    DICT_CIS_OSM[asset_type]["osm_keys"],
                    DICT_CIS_OSM[asset_type]["osm_query"],
                ),
            ]
        )

    else:
        return ImportWarning("feature not in DICT_CIS_OSM. Returning empty gdf")
        gdf = gpd.GeoDataFrame()
    return gdf
