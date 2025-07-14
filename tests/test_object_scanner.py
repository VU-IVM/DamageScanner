from damagescanner import DamageScanner
from .setup import data_path
import geopandas as gpd
import pandas as pd
import rasterio
import rioxarray as rxr


def test_object_scanner_polygon():
    objects = gpd.read_file(data_path / "kampen" / "exposure" / "landuse.gpkg")
    maximum_damage = pd.read_csv(data_path / "kampen" / "vulnerability" / "maxdam.csv")
    curves = pd.read_csv(data_path / "kampen" / "vulnerability" / "curves.csv")

    maximum_damage["object_type"] = maximum_damage["object_type"].astype(str)
    objects["object_type"] = objects["landuse"].map(
        {
            "grass": "230",
            "farmland": "220",
            "forest": "360",
            "industrial": "440",
            "residential": "350",
            "meadow": "220",
            "allotments": "310",
            "retail": "140",
            "cemetery": "530",
            "commercial": "110",
            "construction": "110",
            "farmyard": "120",
            "orchard": "220",
            "reservoir": "550",
            "scrub": "630",
            "static_caravan": "110",
        }
    )

    kampen = DamageScanner(
        feature_data=objects,
        hazard_data=data_path / "kampen" / "hazard" / "1in100_inundation_map.tif",
        curves=curves,
        maxdam=maximum_damage,
    )
    objects["damage_path"] = kampen.calculate()["damage"]

    kampen = DamageScanner(
        feature_data=objects,
        hazard_data=rxr.open_rasterio(
            data_path / "kampen" / "hazard" / "1in100_inundation_map.tif"
        ),
        curves=curves,
        maxdam=maximum_damage,
    )
    objects["damage_rioxarray"] = kampen.calculate()["damage"]

    kampen = DamageScanner(
        feature_data=objects,
        hazard_data=rasterio.open(
            data_path / "kampen" / "hazard" / "1in100_inundation_map.tif", "r"
        ),
        curves=curves,
        maxdam=maximum_damage,
    )
    objects["damage_rasterio"] = kampen.calculate()["damage"]

    assert objects["damage_path"].equals(objects["damage_rioxarray"])
    assert objects["damage_path"].equals(objects["damage_rasterio"])


# def test_object_scanner_line():
#     hazard = rxr.open_rasterio(data_path / "hazard" / "inundation_map.tif")
#     hazard.values[:] = 5
#     objects = gpd.read_file(data_path / "landuse" / "kampen.osm.pbf", layer="lines")
#     objects = objects[
#         objects["highway"].isin(
#             [
#                 "residential",
#                 "secondary",
#                 "tertiary",
#                 "unclassified",
#                 "track",
#             ]
#         )
#     ]
#     objects["object_type"] = "road"
#     objects["maximum_damage"] = 1
#     objects = objects.to_crs(hazard.rio.crs)
#     curves = pd.DataFrame(
#         index=np.linspace(0, 5, 10),
#         data=np.linspace(0, 0.6, 10),
#         columns=["road"],
#     )
#     damage = object_scanner(
#         objects=objects,
#         hazard=hazard,
#         curves=curves,
#     )
#     np.testing.assert_allclose(objects.length * curves["road"].values[-1], damage)

#     # maximum damage of 10, should also result in 10 times the damage
#     objects["maximum_damage"] = 10
#     np.testing.assert_allclose(
#         damage * 10,
#         object_scanner(
#             objects=objects,
#             hazard=hazard,
#             curves=curves,
#         ),
#     )


# def test_object_scanner_mix():
#     hazard = rasterio.open(data_path / "hazard" / "inundation_map.tif")
#     hazard_xarray = rxr.open_rasterio(data_path / "hazard" / "inundation_map.tif")

#     objects = gpd.read_file(data_path / "landuse" / "landuse.shp")
#     objects = objects[objects["landuse"] == "residential"]
#     objects = objects.clip(hazard_xarray.rio.bounds())

#     objects["object_type"] = "residential"
#     objects["maximum_damage"] = 1600

#     # make point data from polygons
#     objects.loc[objects.index[:10], "geometry"] = objects.iloc[:10].centroid

#     curves = pd.read_csv(data_path / "curves" / "curves.csv")[["111"]].rename(
#         columns={"111": "residential"}
#     )

#     objects["damage_rasterio"] = object_scanner(
#         objects=objects,
#         hazard=hazard,
#         curves=curves,
#     )
#     objects.to_file(output_folder / "damaged_objects_rasterio.gpkg", driver="GPKG")

#     objects["damage_xarray"] = object_scanner(
#         objects=objects,
#         hazard=hazard_xarray,
#         curves=curves,
#     )
#     objects.to_file(output_folder / "damaged_objects_xarray.gpkg", driver="GPKG")

#     assert objects["damage_rasterio"].equals(objects["damage_xarray"])

#     objects["maximum_damage"] = 1
#     hazard_xarray.values[:] = (
#         curves.index.max()
#     )  # settings hazard to maximum value of curve in the entire area
#     damage = object_scanner(
#         objects=objects,
#         hazard=hazard_xarray,
#         curves=curves,
#     )
#     # thus the damage should be equal to the maximum damage of the objects
#     np.testing.assert_allclose(curves["residential"].values[-1], damage[:10])
#     np.testing.assert_allclose(
#         curves["residential"].values[-1] * objects.area[10:], damage[10:]
#     )

#     objects["maximum_damage"] = 10

#     # maximum damage of 10, should also result in 10 times the damage
#     np.testing.assert_allclose(
#         damage * 10,
#         object_scanner(
#             objects=objects,
#             hazard=hazard_xarray,
#             curves=curves,
#         ),
#     )
