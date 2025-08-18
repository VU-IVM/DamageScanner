import math

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray as rxr  # noqa: F401
import xarray as xr
from damagescanner.vector import VectorScanner
from pytest import fixture
from shapely.geometry import Polygon


@fixture
def flood_raster():
    x = np.zeros((10, 10), dtype=np.float32)
    x[:5, :5] = 0.5
    x[:5, 5:] = 3
    x[5:, :5] = 0
    x[5:, 5:] = np.nan
    flood_raster = xr.DataArray(
        x,
        coords={
            "x": np.arange(10, dtype=np.int32) + 0.5,  # 0.5 for center of the pixel
            "y": np.arange(10, dtype=np.int32) + 0.5,  # 0.5 for center of the pixel
        },
        dims=["y", "x"],
    )
    flood_raster.attrs["crs"] = "EPSG:32631"
    flood_raster.attrs["_FillValue"] = np.nan
    return flood_raster


@fixture
def buildings():
    data = {
        "object_type": [
            "residential",
            "residential",
            "commercial",
            "residential",
            "residential",
            "residential",
            "residential",
            "residential",
        ],
        "maximum_damage": [100, 200, 300, 100, 100, 100, 100, 100],
    }
    polygons = [
        Polygon([(0, 0), (0, 1), (1, 1), (1, 0)]),
        Polygon([(0, 0), (0, 2), (2, 2), (2, 0)]),
        Polygon([(0, 0), (0, 2), (2, 2), (2, 0)]),
        Polygon([(9, 0), (9, 1), (10, 1), (10, 0)]),
        Polygon([(4.5, 0), (4.5, 1), (5.5, 1), (5.5, 0)]),
        Polygon([(0, 9), (0, 10), (1, 10), (1, 9)]),
        Polygon([(9, 9), (9, 10), (10, 10), (10, 9)]),
        Polygon([(4.5, 4.5), (4.5, 5.5), (5.5, 5.5), (5.5, 4.5)]),
    ]
    data["geometry"] = polygons

    gdf = gpd.GeoDataFrame(data, crs="EPSG:32631", geometry="geometry")
    return gdf


@fixture
def vulnerability_curves():
    return pd.DataFrame(
        {
            "residential": [0.0, 0.2, 0.3],
            "commercial": [0.0, 0.4, 0.6],
        },
        index=[0, 1, 2],
    )


def test_vector_scanner(flood_raster, buildings, vulnerability_curves):
    damage = VectorScanner(
        feature_file=buildings,
        hazard_file=flood_raster,
        curve_path=vulnerability_curves,
        gridded=False,
    )
    assert math.isclose(
        damage["damage"].iloc[0], 10.0
    )  # 1m2, .5 hazard severity, residential, max_damage 100 > 0.1 damage ratio > damage of 10
    assert math.isclose(
        damage["damage"].iloc[1], 80.0
    )  # 4m2, .5 hazard severity, residential, max_damage 200 > 0.1 damage ratio > damage of 80
    assert math.isclose(
        damage["damage"].iloc[2], 240.0
    )  # 4m2, .5 hazard severity, commercial, max_damage 300 > 0.2 damage ratio > damage of 240
    assert math.isclose(
        damage["damage"].iloc[3], 30.0
    )  # 1m2, 3 hazard severity, residential, max_damage 100 > 0.3 damage ratio (max of curve) > damage of 30

    # 1m2, half 0.5 hazard severity + half 3 hazard severity, residential, max_damage 100
    # > 0.1 + 0.3 damage ratio (max of curve) > damage of 30
    assert math.isclose(damage["damage"].iloc[4], 20.0)

    assert math.isclose(
        damage["damage"].iloc[5], 0.0
    )  # 1m2, no hazard severity, residential, max_damage 100
    assert math.isclose(
        damage["damage"].iloc[6], 0.0
    )  # 1m2, hazard is nan, residential, max_damage 100

    assert math.isclose(
        damage["damage"].iloc[7], 10.0
    )  # 1m2, 0.25 x 0, 0.25 x nan, 0.25 x 0.5, 0.25 x 3 hazard severity, residential, max_damage 100 > damage of 10
