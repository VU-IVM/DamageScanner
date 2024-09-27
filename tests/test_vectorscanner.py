from damagescanner.core import VectorScanner
from .setup import data_path, output_folder
import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
import rioxarray as rxr


def test_VectorScanner_polygon():
    hazard = rasterio.open(data_path / "hazard" / "inundation_map.tif")
    objects = gpd.read_file(data_path / "landuse" / "landuse.shp")
    objects = objects[objects["landuse"] == "residential"]
    curve = pd.read_csv(data_path / "curves" / "curves.csv")["111"]
    maximum_damage = 1600

    objects["damage_rasterio"] = VectorScanner(
        objects=objects,
        hazard=hazard,
        curve=curve,
        maximum_damage=maximum_damage,
    )
    objects.to_file(output_folder / "damaged_objects_rasterio.gpkg", driver="GPKG")

    hazard = rxr.open_rasterio(data_path / "hazard" / "inundation_map.tif")
    objects["damage_xarray"] = VectorScanner(
        objects=objects,
        hazard=hazard,
        curve=curve,
        maximum_damage=maximum_damage,
    )
    objects.to_file(output_folder / "damaged_objects_xarray.gpkg", driver="GPKG")

    assert objects["damage_rasterio"].equals(objects["damage_xarray"])

    objects = objects.clip(hazard.rio.bounds())
    hazard.values[:] = (
        curve.index.max()
    )  # settings hazard to maximum value of curve in the entire area
    damage = VectorScanner(
        objects=objects,
        hazard=hazard,
        curve=curve,
        maximum_damage=1,  # maximum damage of 1
    )
    objects.area * curve.values[-1], damage
    # thus the damage should be equal to the area of the objects times the maximum value of the curve
    np.testing.assert_allclose(objects.area * curve.values[-1], damage)

    # maximum damage of 10, should also result in 10 times the damage
    np.testing.assert_allclose(
        damage * 10,
        VectorScanner(
            objects=objects,
            hazard=hazard,
            curve=curve,
            maximum_damage=10,
        ),
    )


def test_VectorScanner_line():
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
    objects = objects.to_crs(hazard.rio.crs)
    curve = pd.Series(index=np.linspace(0, 5, 10), data=np.linspace(0, 0.6, 10))
    damage = VectorScanner(
        objects=objects,
        hazard=hazard,
        curve=curve,
        maximum_damage=1,
    )
    np.testing.assert_allclose(objects.length * curve.values[-1], damage)

    # maximum damage of 10, should also result in 10 times the damage
    np.testing.assert_allclose(
        damage * 10,
        VectorScanner(
            objects=objects,
            hazard=hazard,
            curve=curve,
            maximum_damage=10,
        ),
    )
