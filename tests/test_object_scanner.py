from damagescanner.core import object_scanner
from .setup import data_path, output_folder
import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
import rioxarray as rxr


def test_object_scanner_polygon():
    hazard = rasterio.open(data_path / "hazard" / "inundation_map.tif")
    objects = gpd.read_file(data_path / "landuse" / "landuse.shp")
    objects = objects[objects["landuse"] == "residential"]
    objects["object_type"] = "residential"
    objects["maximum_damage"] = 1600
    curves = pd.read_csv(data_path / "curves" / "curves.csv")[["111"]].rename(
        columns={"111": "residential"}
    )

    objects["damage_rasterio"] = object_scanner(
        objects=objects,
        hazard=hazard,
        curves=curves,
    )
    objects.to_file(output_folder / "damaged_objects_rasterio.gpkg", driver="GPKG")

    hazard = rxr.open_rasterio(data_path / "hazard" / "inundation_map.tif")
    objects["damage_xarray"] = object_scanner(
        objects=objects,
        hazard=hazard,
        curves=curves,
    )
    objects.to_file(output_folder / "damaged_objects_xarray.gpkg", driver="GPKG")

    assert objects["damage_rasterio"].equals(objects["damage_xarray"])

    objects = objects.clip(hazard.rio.bounds())
    objects["maximum_damage"] = 1
    hazard.values[:] = (
        curves.index.max()
    )  # settings hazard to maximum value of curve in the entire area
    damage = object_scanner(
        objects=objects,
        hazard=hazard,
        curves=curves,
    )
    # thus the damage should be equal to the area of the objects times the maximum value of the curve
    np.testing.assert_allclose(objects.area * curves["residential"].values[-1], damage)

    objects["maximum_damage"] = 10

    # maximum damage of 10, should also result in 10 times the damage
    np.testing.assert_allclose(
        damage * 10,
        object_scanner(
            objects=objects,
            hazard=hazard,
            curves=curves,
        ),
    )


def test_object_scanner_line():
    hazard = rxr.open_rasterio(data_path / "hazard" / "inundation_map.tif")
    hazard.values[:] = 5
    objects = gpd.read_file(data_path / "landuse" / "kampen.osm.pbf", layer="lines")
    objects = objects[
        objects["highway"].isin(
            [
                "residential",
                "secondary",
                "tertiary",
                "unclassified",
                "track",
            ]
        )
    ]
    objects["object_type"] = "road"
    objects["maximum_damage"] = 1
    objects = objects.to_crs(hazard.rio.crs)
    curves = pd.DataFrame(
        index=np.linspace(0, 5, 10),
        data=np.linspace(0, 0.6, 10),
        columns=["road"],
    )
    damage = object_scanner(
        objects=objects,
        hazard=hazard,
        curves=curves,
    )
    np.testing.assert_allclose(objects.length * curves["road"].values[-1], damage)

    # maximum damage of 10, should also result in 10 times the damage
    objects["maximum_damage"] = 10
    np.testing.assert_allclose(
        damage * 10,
        object_scanner(
            objects=objects,
            hazard=hazard,
            curves=curves,
        ),
    )


def test_object_scanner_mix():
    hazard = rasterio.open(data_path / "hazard" / "inundation_map.tif")
    hazard_xarray = rxr.open_rasterio(data_path / "hazard" / "inundation_map.tif")

    objects = gpd.read_file(data_path / "landuse" / "landuse.shp")
    objects = objects[objects["landuse"] == "residential"]
    objects = objects.clip(hazard_xarray.rio.bounds())

    objects["object_type"] = "residential"
    objects["maximum_damage"] = 1600

    # make point data from polygons
    objects.loc[objects.index[:10], "geometry"] = objects.iloc[:10].centroid

    curves = pd.read_csv(data_path / "curves" / "curves.csv")[["111"]].rename(
        columns={"111": "residential"}
    )

    objects["damage_rasterio"] = object_scanner(
        objects=objects,
        hazard=hazard,
        curves=curves,
    )
    objects.to_file(output_folder / "damaged_objects_rasterio.gpkg", driver="GPKG")

    objects["damage_xarray"] = object_scanner(
        objects=objects,
        hazard=hazard_xarray,
        curves=curves,
    )
    objects.to_file(output_folder / "damaged_objects_xarray.gpkg", driver="GPKG")

    assert objects["damage_rasterio"].equals(objects["damage_xarray"])

    objects["maximum_damage"] = 1
    hazard_xarray.values[:] = (
        curves.index.max()
    )  # settings hazard to maximum value of curve in the entire area
    damage = object_scanner(
        objects=objects,
        hazard=hazard_xarray,
        curves=curves,
    )
    # thus the damage should be equal to the maximum damage of the objects
    np.testing.assert_allclose(curves["residential"].values[-1], damage[:10])
    np.testing.assert_allclose(
        curves["residential"].values[-1] * objects.area[10:], damage[10:]
    )

    objects["maximum_damage"] = 10

    # maximum damage of 10, should also result in 10 times the damage
    np.testing.assert_allclose(
        damage * 10,
        object_scanner(
            objects=objects,
            hazard=hazard_xarray,
            curves=curves,
        ),
    )
