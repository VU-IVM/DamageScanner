import geopandas as gpd
import pandas as pd
import xarray as xr
import numpy as np
import shapely
import json

import osm_flex.download as dl
import osm_flex.extract as ex
from osm_flex.config import OSM_DATA_DIR, DICT_GEOFABRIK
from osm_flex.simplify import remove_contained_points

from tqdm import tqdm

from pathlib import Path

from damagescanner.vector import (
    _get_damage_per_element,
)
from .utils import fetch_and_save


def download_osm(geometry, overwrite=False):
    """Download OSM data from geofabrik
    using the geofabrik index to find the regions that intersect with the geometry

    Parameters
    ----------
    geometry : shapely.geometry

    """
    index_file = OSM_DATA_DIR / "geofabrik_region_index.geojson"
    fetch_and_save(
        "https://download.geofabrik.de/index-v1.json", index_file, overwrite=overwrite
    )

    index = gpd.read_file(index_file)
    # remove Dach region as all individual regions within dach countries are also in the index
    index = index[index["id"] != "dach"]

    # find all regions that intersect with the bbox
    intersecting_regions = index[index.intersects(geometry)]
    intersecting_regions = intersecting_regions[
        intersecting_regions["id"].apply(
            lambda x: x not in intersecting_regions["parent"].tolist()
        )
    ]

    # download all regions
    file_paths = []
    for _, row in tqdm(intersecting_regions.iterrows()):
        url = row["urls"]["pbf"]
        filepath = OSM_DATA_DIR / url.split("/")[-1]
        file_paths.append(filepath)
        fetch_and_save(url, filepath, overwrite=overwrite)

    return file_paths


def country_download(iso3):
    """Download country data from geofabrik

    Parameters
    ----------
    iso3 : str
        ISO3 country code

    Returns
    -------
    data_loc : pathlib.Path
        Path to downloaded data
    """

    dl.get_country_geofabrik(iso3)
    data_loc = OSM_DATA_DIR.joinpath(f"{DICT_GEOFABRIK[iso3][1]}-latest.osm.pbf")
    return data_loc


def read_hazard_data(data_path, hazard_type):
    """Read hazard data from data_path

    Parameters
    ----------
    data_path : pathlib.Path
        Path to data folder
    hazard_type : str
        Hazard type to read

    Returns
    -------
    list
        List of paths to hazard data
    """

    if hazard_type == "fluvial":
        hazard_data = (
            data_path / "Floods" / "Jamaica" / "fluvial_undefended"
        )  # need to make country an input
        return list(hazard_data.iterdir())

    elif hazard_type == "pluvial":
        hazard_data = (
            data_path / "Floods" / "Jamaica" / "pluvial"
        )  # need to make country an input
        return list(hazard_data.iterdir())

    elif hazard_type == "windstorm":
        hazard_data = data_path / "Windstorms"
        return list(hazard_data.iterdir())

    elif hazard_type == "earthquake":
        hazard_data = data_path / "Earthquakes"
        return list(hazard_data.iterdir())

    elif hazard_type == "landslides":
        hazard_data = data_path / "Landslides"
        return list(hazard_data.iterdir())


def read_vul_maxdam(data_path, hazard_type, infra_type):
    """Read vulnerability and maximum damage data

    Parameters
    ----------
    data_path : pathlib.Path
        Path to data folder
    hazard_type : str
        Hazard type to read
    infra_type : str
        Infrastructure type to read

    Returns
    -------
    infra_curves : pandas.DataFrame
        Vulnerability curves
    infra_maxdam : pandas.DataFrame
        Maximum damage data
    """

    vul_data = data_path / "Vulnerability"

    if hazard_type in ["pluvial", "fluvial"]:
        curves = pd.read_excel(
            vul_data
            / "Table_D2_Multi-Hazard_Fragility_and_Vulnerability_Curves_V1.0.0.xlsx",
            sheet_name="F_Vuln_Depth",
            index_col=[0],
            header=[0, 1, 2, 3, 4],
        )
    elif hazard_type == "windstorm":
        curves = pd.read_excel(
            vul_data
            / "Table_D2_Multi-Hazard_Fragility_and_Vulnerability_Curves_V1.0.0.xlsx",
            sheet_name="W_Vuln_V10m",
            index_col=[0],
            header=[0, 1, 2, 3, 4],
        )

    infra_curves = curves.loc[
        :,
        curves.columns.get_level_values("Infrastructure description")
        .str.lower()
        .str.contains(infra_type),
    ]

    maxdam = pd.read_excel(
        vul_data / "Table_D3_Costs_V1.0.0.xlsx",
        sheet_name="Cost_Database",
        index_col=[0],
    )
    infra_maxdam = maxdam.loc[
        maxdam.index.get_level_values("Infrastructure description")
        .str.lower()
        .str.contains(infra_type),
        "Amount",
    ].dropna()
    infra_maxdam = infra_maxdam[pd.to_numeric(infra_maxdam, errors="coerce").notnull()]

    return infra_curves, infra_maxdam


def read_flood_map(flood_map_path):
    """Read flood map data

    Parameters
    ----------
    flood_map_path : pathlib.Path
        Path to flood map data

    Returns
    -------
    flood_map_vector : geopandas.GeoDataFrame
        Flood map data
    """

    flood_map = xr.open_dataset(flood_map_path, engine="rasterio")

    flood_map_vector = (
        flood_map["band_data"].to_dataframe().reset_index()
    )  # transform to dataframe

    # remove data that will not be used
    flood_map_vector = flood_map_vector.loc[
        (flood_map_vector.band_data > 0) & (flood_map_vector.band_data < 100)
    ]

    # create geometry values and drop lat lon columns
    flood_map_vector["geometry"] = [
        shapely.points(x)
        for x in list(zip(flood_map_vector["x"], flood_map_vector["y"]))
    ]
    flood_map_vector = flood_map_vector.drop(["x", "y", "band", "spatial_ref"], axis=1)

    # drop all non values to reduce size
    flood_map_vector = flood_map_vector.loc[
        ~flood_map_vector["band_data"].isna()
    ].reset_index(drop=True)

    # and turn them into squares again:
    flood_map_vector.geometry = shapely.buffer(
        flood_map_vector.geometry, distance=0.00083 / 2, cap_style="square"
    ).values  # distance should be made an input still!

    return flood_map_vector


def read_windstorm_map(windstorm_map_path, bbox):
    # load data from NetCDF file
    with xr.open_dataset(windstorm_map_path) as ds:
        # convert data to WGS84 CRS
        ds.rio.write_crs(4326, inplace=True)
        ds = ds.rio.clip_box(minx=bbox[0], miny=bbox[1], maxx=bbox[2], maxy=bbox[3])
        # ds['band_data'] = ds['band_data']/0.88*1.11 #convert 10-min sustained wind speed to 3-s gust wind speed

        ds_vector = (
            ds["band_data"].to_dataframe().reset_index()
        )  # transform to dataframe

        # remove data that will not be used
        ds_vector = ds_vector.loc[
            (ds_vector.band_data > 0) & (ds_vector.band_data < 100)
        ]

        # create geometry values and drop lat lon columns
        ds_vector["geometry"] = [
            shapely.points(x) for x in list(zip(ds_vector["x"], ds_vector["y"]))
        ]
        ds_vector = ds_vector.drop(["x", "y", "band", "spatial_ref"], axis=1)
        ds_vector["geometry"] = shapely.buffer(
            ds_vector.geometry, distance=0.1 / 2, cap_style="square"
        ).values

        return ds_vector


def country_infrastructure_hazard(data_path, country_code, infra_type, hazard_type):
    # get country osm data
    data_loc = country_download(country_code)

    # get infrastructure data:
    assets = ex.extract_cis(data_loc, infra_type)

    # convert assets to epsg3857
    assets = gpd.GeoDataFrame(assets).set_crs(4326).to_crs(3857)

    if infra_type == "road":
        assets = assets.loc[assets.geometry.geom_type == "LineString"]
        assets = assets.rename(columns={"highway": "asset"})
    elif infra_type == "rail":
        assets = assets.loc[assets.geometry.geom_type == "LineString"]
        assets = assets.rename(columns={"railway": "asset"})
    elif infra_type == "education":
        assets = assets.rename(columns={"building": "asset"})
        assets = assets.reset_index(drop=True)
        assets = remove_contained_points(assets)

        # convert points to polygons
        assets.loc[assets.geom_type == "Point", "geometry"] = assets.loc[
            assets.geom_type == "Point"
        ].buffer(
            distance=np.sqrt(
                assets.loc[assets.geom_type == "MultiPolygon"].area.median()
            )
            / 2,
            cap_style="square",
        )

    elif infra_type == "air":
        assets = assets.rename(columns={"aeroway": "asset"})

    # create dicts for quicker lookup
    geom_dict = assets["geometry"].to_dict()
    type_dict = assets["asset"].to_dict()

    # read hazard data
    hazard_data_list = read_hazard_data(data_path, hazard_type)

    # read vulnerability and maxdam data:
    infra_curves, maxdams = read_vul_maxdam(data_path, hazard_type, infra_type)

    # start analysis
    print(
        f"{country_code} runs for {infra_type} for {hazard_type} for {len(hazard_data_list)} maps for {len(infra_curves.T)*len(maxdams)} combinations"
    )

    if hazard_type in ["windstorm", "earthquake", "landslide"]:
        # load country geometry file and create geometry to clip
        ne_countries = gpd.read_file(
            data_path / "natural_earth" / "ne_10m_admin_0_countries.shp"
        )
        bbox = (
            ne_countries.loc[ne_countries["ISO_A3"] == country_code]
            .geometry.envelope.values[0]
            .bounds
        )

    collect_output = {}
    for (
        single_footprint
    ) in hazard_data_list:  # tqdm(hazard_data_list,total=len(hazard_data_list)):
        hazard_name = single_footprint.parts[-1].split(".")[0]

        # load hazard map
        if hazard_type in ["pluvial", "fluvial"]:
            hazard_map = read_flood_map(single_footprint)
        elif hazard_type in ["windstorm"]:
            hazard_map = read_windstorm_map(single_footprint, bbox)
        elif hazard_type in ["earthquake"]:
            hazard_map = read_earthquake_map(single_footprint)
        elif hazard_type in ["landslide"]:
            hazard_map = read_landslide_map(single_footprint)

        # convert hazard data to epsg 3857
        hazard_map = gpd.GeoDataFrame(hazard_map).set_crs(4326).to_crs(3857)

        # overlay assets:
        overlay_assets = pd.DataFrame(
            overlay_hazard_assets(hazard_map, buffer_assets(assets)).T,
            columns=["asset", "hazard_point"],
        )

        # convert dataframe to numpy array
        hazard_numpified = hazard_map.to_numpy()

        for infra_curve in infra_curves:
            # get curves
            curve = infra_curves[infra_curve[0]]
            hazard_intensity = curve.index.values
            fragility_values = (
                np.nan_to_num(curve.values, nan=(np.nanmax(curve.values)))
            ).flatten()

            for maxdam in maxdams:
                collect_inb = {}
                for asset in tqdm(
                    overlay_assets.groupby("asset"),
                    total=len(overlay_assets.asset.unique()),
                ):  # group asset items for different hazard points per asset and get total number of unique assets
                    asset_type = type_dict[asset[0]]

                    if np.max(fragility_values) == 0:
                        collect_inb[asset_type] = 0
                    else:
                        asset_geom = geom_dict[asset[0]]

                    collect_inb[asset[0]] = asset_type, _get_damage_per_element(
                        asset,
                        hazard_numpified,
                        asset_geom,
                        hazard_intensity,
                        fragility_values,
                        maxdam,
                    )
                collect_output[hazard_name, infra_curve[0], maxdam] = (
                    pd.DataFrame(
                        pd.DataFrame(collect_inb).T.values,
                        columns=["asset_type", "damage"],
                    )
                    .groupby("asset_type")
                    .sum()
                )

    return collect_output


if __name__ == "__main__":
    p = Path("..")
    data_path = Path(Path.home().parts[0]) / "Data"
    country_code = "JAM"
    infra_type = "road"
    hazard_type = "fluvial"

    country_infrastructure_hazard(data_path, country_code, infra_type, hazard_type)
