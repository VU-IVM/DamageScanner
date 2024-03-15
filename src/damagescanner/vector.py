# Get all the needed modules
import xarray as xr
import numpy as np
import shapely
import pandas as pd
import geopandas as gpd
import pyproj
from tqdm import tqdm
from pathlib import PurePath
import osm_flex.extract as ex

from .utils import _check_output_path, _check_scenario_name


def landuse(self):
    """ """
    extracted_exposure = ex.extract(self.exposure_path, "multipolygons", ["landuse"])

    extracted_exposure.geometry = shapely.make_valid(extracted_exposure.geometry)

    return extracted_exposure


def buildings(self):
    """ """
    extracted_exposure = ex.extract(self.exposure_path, "multipolygons", ["building"])

    extracted_exposure.geometry = shapely.make_valid(extracted_exposure.geometry)

    return extracted_exposure


def cis(self, infra_type):
    """ """
    extracted_exposure = ex.extract_cis(self.exposure_path, infra_type)

    if infra_type == "road":
        extracted_exposure = extracted_exposure.loc[
            extracted_exposure.geometry.geom_type == "LineString"
        ]

    extracted_exposure.geometry = shapely.make_valid(extracted_exposure.geometry)

    return extracted_exposure


def cis_all(self, to_exclude=[]):
    """ """
    cis = [
        "healthcare",
        "education",
        "gas",
        "oil",
        "telecom",
        "water",
        "wastewater",
        "power",
        "rail",
        "road",
        "air",
    ]
    pass


def download(self, country_code="JAM"):
    """ """
    pass


def _overlay_hazard_assets(df_ds, assets):
    """
    Overlay hazard assets on a dataframe of spatial geometries.
    Arguments:
        *df_ds*: GeoDataFrame containing the spatial geometries of the hazard data.
        *assets*: GeoDataFrame containing the infrastructure assets.
    Returns:
        *geopandas.GeoSeries*: A GeoSeries containing the spatial geometries of df_ds that intersect with the infrastructure assets.
    """
    # overlay
    hazard_tree = shapely.STRtree(df_ds.geometry.values)
    if (shapely.get_type_id(assets.iloc[0].geometry) == 3) | (
        shapely.get_type_id(assets.iloc[0].geometry) == 6
    ):  # id types 3 and 6 stand for polygon and multipolygon
        return hazard_tree.query(assets.geometry, predicate="intersects")
    else:
        return hazard_tree.query(assets.buffered, predicate="intersects")


def _buffer_assets(assets, buffer_size=0.00083):
    """
    Buffer spatial assets in a GeoDataFrame.
    Arguments:
        *assets*: GeoDataFrame containing spatial geometries to be buffered.
        *buffer_size* (float, optional): The distance by which to buffer the geometries. Default is 0.00083.
    Returns:
        *GeoDataFrame*: A new GeoDataFrame with an additional 'buffered' column containing the buffered geometries.
    """
    assets["buffered"] = shapely.buffer(assets.geometry.values, distance=buffer_size)
    return assets


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


def _extract_value_other_gdf(x, gdf, col_name):
    """
    Function to extract value from column from other GeoDataFrame

    Arguments:
        *x* : row of main GeoDataFrame.

        *gdf* : geopandas GeoDataFrame from which we want to extract values.

        *col_name* : the column name from which we want to get the value.


    """
    try:
        return gdf.loc[list(gdf.sindex.intersection(x.geometry.bounds))][
            col_name
        ].values[0]
    except:
        return None


def _reproject(df_ds, current_crs="epsg:4326", approximate_crs="epsg:3035"):
    """
    Function to reproject a GeoDataFrame to a different CRS.

    Args:
        df_ds (pandas DataFrame): _description_
        current_crs (str, optional): _description_. Defaults to "epsg:4326".
        approximate_crs (str, optional): _description_. Defaults to "epsg:3035".

    Returns:
        _type_: _description_
    """
    geometries = df_ds["geometry"]
    coords = shapely.get_coordinates(geometries)
    transformer = pyproj.Transformer.from_crs(
        current_crs, approximate_crs, always_xy=True
    )
    new_coords = transformer.transform(coords[:, 0], coords[:, 1])

    return shapely.set_coordinates(geometries.copy(), np.array(new_coords).T)


def _get_damage_per_element(
    asset,
    hazard_numpified,
    asset_geom,
    hazard_intensity,
    fragility_values,
    maxdam_asset,
):
    """
    Calculate damage for a given asset based on hazard information.
    Arguments:
        *asset*: Tuple containing information about the asset. It includes:
            - Index or identifier of the asset (asset[0]).
            - Asset-specific information, including hazard points (asset[1]['hazard_point']).
        *flood_numpified*: NumPy array representing flood hazard information.
        *asset_geom*: Shapely geometry representing the spatial coordinates of the asset.
        *curve*: Pandas DataFrame representing the curve for the asset type.
        *maxdam*: Maximum damage value.
    Returns:
        *tuple*: A tuple containing the asset index or identifier and the calculated damage.
    """

    # find the exact hazard overlays:
    get_hazard_points = hazard_numpified[asset[1]["hazard_point"].values]
    get_hazard_points[shapely.intersects(get_hazard_points[:, 1], asset_geom)]

    # estimate damage
    if len(get_hazard_points) == 0:  # no overlay of asset with hazard
        return 0

    else:
        if asset_geom.geom_type == "LineString":
            overlay_meters = shapely.length(
                shapely.intersection(get_hazard_points[:, 1], asset_geom)
            )  # get the length of exposed meters per hazard cell
            return np.sum(
                (
                    np.interp(
                        np.float16(get_hazard_points[:, 0]),
                        hazard_intensity,
                        fragility_values,
                    )
                )
                * overlay_meters
                * maxdam_asset
            )  # return asset number, total damage for asset number (damage factor * meters * max. damage)
        elif asset_geom.geom_type in ["MultiPolygon", "Polygon"]:
            overlay_m2 = shapely.area(
                shapely.intersection(get_hazard_points[:, 1], asset_geom)
            )
            return np.sum(
                (
                    np.interp(
                        np.float16(get_hazard_points[:, 0]),
                        hazard_intensity,
                        fragility_values,
                    )
                )
                * overlay_m2
                * maxdam_asset
            )
        elif asset_geom.geom_type == "Point":
            return np.sum(
                (
                    np.interp(
                        np.float16(get_hazard_points[:, 0]),
                        hazard_intensity,
                        fragility_values,
                    )
                )
                * maxdam_asset
            )


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
    lat_col="y",
    lon_col="x",
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
            hazard = hazard = hazard.loc[
                (hazard.band_data > 0) & (hazard.band_data < 100)
            ]

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
            with xr.open_dataset(hazard_file) as ds:

                # get bbox of the exposure data
                bbox = exposure.to_crs(haz_crs).total_bounds

                # convert data to WGS84 CRS
                ds.rio.write_crs(haz_crs, inplace=True)
                ds.rio.set_spatial_dims(x_dim=lon_col, y_dim=lat_col, inplace=True)

                hazard_map = ds.rio.clip_box(
                    minx=bbox[0], miny=bbox[1], maxx=bbox[2], maxy=bbox[3]
                )

            # convert to dataframe
            hazard = hazard_map[hazard_col].to_dataframe().reset_index()

            # drop all non values and below zeros to reduce size. Might cause issues for cold waves
            # hazard = hazard.loc[~(hazard[hazard_col].isna() | hazard[hazard_col]<=0)].reset_index(drop=True)
            hazard = hazard.loc[~hazard[hazard_col].isna()].reset_index(drop=True)

            # create geometry values and drop lat lon columns
            hazard["geometry"] = shapely.points(
                np.array(list(zip(hazard[lon_col], hazard[lat_col])))
            )
            hazard = hazard[[hazard_col, "geometry"]]

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
        exposure.geometry = _reproject(
            exposure, current_crs=f"epsg:{exp_crs}", approximate_crs=f"epsg:{haz_crs}"
        )
    elif exp_crs_unit_name == "degree" and haz_crs_unit_name == "degree":
        exposure.geometry = _reproject(
            exposure, current_crs=f"epsg:{exp_crs}", approximate_crs="epsg:3857"
        )

    if haz_crs_unit_name == "degree" and exp_crs_unit_name == "metre":
        hazard.geometry = _reproject(
            hazard, current_crs=haz_crs, approximate_crs=f"epsg:{exp_crs}"
        )
    elif haz_crs_unit_name == "degree" and exp_crs_unit_name == "degree":
        hazard.geometry = _reproject(
            hazard, current_crs=haz_crs, approximate_crs="epsg:3857"
        )

    # rename inundation colum, set values to centimeters if required and to integers:
    hazard = hazard.rename(columns={hazard_col: "haz_val"})
    if not centimeters:
        hazard["haz_val"] = hazard.haz_val * 100
        hazard["haz_val"] = hazard.haz_val.astype(int)

    # set element column to element_type:
    if "landuse" in exposure.columns:
        exposure = exposure.rename({"landuse": "element_type"}, axis=1)

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

    # create dicts for quicker lookup
    geom_dict = exposure["geometry"].to_dict()
    type_dict = exposure["element_type"].to_dict()

    # overlay hazard and exposure data
    if (shapely.get_type_id(exposure.iloc[0].geometry) == 3) | (
        shapely.get_type_id(exposure.iloc[0].geometry) == 6
    ):
        overlay_assets = pd.DataFrame(
            _overlay_hazard_assets(hazard, exposure).T,
            columns=["asset", "hazard_point"],
        )
    else:
        overlay_assets = pd.DataFrame(
            _overlay_hazard_assets(hazard, _buffer_assets(exposure)).T,
            columns=["asset", "hazard_point"],
        )

    # convert dataframe to numpy array
    hazard_numpified = hazard.to_numpy()

    # prepare calculations
    hazard_intensity = curves.index.values

    # Perform calculation
    collect_damage = {}
    for asset in tqdm(
        overlay_assets.groupby("asset"), total=len(overlay_assets.asset.unique())
    ):  # group asset items for different hazard points per asset and get total number of unique assets
        element_type = type_dict[asset[0]]
        element_geom = geom_dict[asset[0]]

        curve = curves[element_type].values
        fragility_values = (np.nan_to_num(curve, nan=(np.nanmax(curve)))).flatten()
        maxdam_element = maxdam[element_type]

        collect_damage[asset[0]] = element_type, _get_damage_per_element(
            asset,
            hazard_numpified,
            element_geom,
            hazard_intensity,
            fragility_values,
            maxdam_element,
        )

    # Merge results
    damaged_objects = (
        pd.DataFrame(
            pd.DataFrame(collect_damage).T.values, columns=["element_type", "damage"]
        )
        .groupby("element_type")
        .sum()
    )

    # Save output
    if save == True:
        # requires adding output_path and scenario_name to function call
        # If output path is not defined, will place file in current directory
        output_path = _check_output_path(kwargs)
        scenario_name = _check_scenario_name(kwargs)
        path_prefix = PurePath(output_path, scenario_name)

        damage_fn = f"{path_prefix}_damages.csv"
        damaged_objects.to_csv(damage_fn)
        return damaged_objects

    else:
        return damaged_objects
