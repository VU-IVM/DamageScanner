import os
import pandas as pd
import geopandas as gpd

DICT_CIS_OSM = {
    "roads": {
        "osm_keys": ["highway", "name", "maxspeed", "lanes", "surface"],
        "osm_query": [
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
        ],
    },
    "main_roads": {
        "osm_keys": ["highway", "name", "maxspeed", "lanes", "surface"],
        "osm_query": [
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
        ],
    },
    "rail": {
        "osm_keys": ["railway", "name", "gauge", "electrified", "voltage"],
        "osm_query": ["rail", "narrow_gauge"],
    },
    "air": {
        "osm_keys": ["aeroway", "name"],
        "osm_query": ["aerodrome", "terminal", "runway"],
    },
    "telecom": {
        "osm_keys": ["man_made", "tower_type", "name"],
        "osm_query": """tower_type='communication' or man_made='mast' or man_made='communications_tower'""",
    },
    "water_supply": {
        "osm_keys": ["man_made", "name"],
        "osm_query": [
            "water_works",
            "water_well",
            "water_tower",
            "reservoir_covered",
            "storage_tank",
        ],
    },
    "waste_solid": {
        "osm_keys": ["amenity", "name"],
        "osm_query": ["waste_transfer_station"],
    },
    "waste_water": {
        "osm_keys": ["man_made", "name"],
        "osm_query": ["wastewater_plant"],
    },
    "education": {
        "osm_keys": ["amenity", "building", "name"],
        "osm_query": """building='school' or amenity='school' or
                             building='kindergarten' or 
                             amenity='kindergarten' or
                             building='college' or amenity='college' or
                             building='university' or amenity='university' or
                             building='library' or amenity='library'""",
    },
    "healthcare": {
        "osm_keys": ["amenity", "building", "healthcare", "name"],
        "osm_query": ["hospital", "clinic", "doctors", "dentist", "pharmacy"],
    },
    #  or healthcare='hospital' or
    #                  building='hospital' or building='clinic' or
    #                  amenity='clinic' or healthcare='clinic' or
    #                  amenity='doctors' or healthcare='doctors' or
    #                  amenity='dentist' or amenity='pharmacy' or
    #                  healthcare='pharmacy' or healthcare='dentist' or
    #                  healthcare='physiotherapist' or healthcare='alternative' or
    #                  healthcare='laboratory' or healthcare='optometrist' or
    #                  healthcare='rehabilitation' or healthcare='blood_donation' or
    #                  healthcare='birthing_center'
    #                  """},
    "power": {
        "osm_keys": ["power", "voltage", "utility", "name"],
        "osm_query": [
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
        ],
    },
    "gas": {
        "osm_keys": ["man_made", "pipeline", "utility", "name"],
        "osm_query": """(man_made='pipeline' and substance='gas') or
                              (pipeline='substation' and substance='gas') or
                              (man_made='storage_tank' and content='gas') or
                              utility='gas'""",
    },
    "oil": {
        "osm_keys": ["pipeline", "man_made", "amenity", "name"],
        "osm_query": """(pipeline='substation' and substance='oil') or
                              (man_made='pipeline' and substance='oil') or
                              man_made='petroleum_well' or 
                              man_made='oil_refinery' or
                              amenity='fuel'""",
    },
    "wastewater": {
        "osm_keys": ["man_made", "amenity", "name"],
        "osm_query": """amenity='waste_transfer_station' or man_made='wastewater_plant'""",
    },
    "buildings": {
        "osm_keys": ["building", "amenity", "name"],
        "osm_query": [
            "yes",
            "house",
            "residential",
            "detached",
            "hut",
            "industrial",
            "shed",
            "apartments",
        ],
    },
    # 'osm_query' : """building='yes' or building='house' or
    #                 building='residential' or building='detached' or
    #                 building='hut' or building='industrial' or
    #                 building='shed' or building='apartments'"""}
}


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
    exposure = gpd.read_file(osm_path, layer=geom_type, engine="pyogrio")

    exposure = exposure[exposure[osm_keys[0]].isin(osm_query)]

    exposure = exposure[["osm_id", "geometry"] + osm_keys]

    exposure.rename(columns={osm_keys[0]: "object_type"}, inplace=True)

    return exposure


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
