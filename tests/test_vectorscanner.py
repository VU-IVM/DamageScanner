from damagescanner.core import VectorScanner, calculate_damage_per_object
from .setup import data_path, output_folder
import geopandas as gpd
import pandas as pd
import rasterio
import rioxarray as rxr
from exactextract import exact_extract


def test_calculate_damage_per_object():
    objects = gpd.read_file(data_path / "landuse" / "landuse.shp")
    objects = objects[objects['landuse'].isin([
        'residential',
        'grass',
        'industrial',
        'farmland',
        'forest',
        'orchard',
    ])]
    objects['landuse'] = objects['landuse'].map({
        'residential': "111",
        'grass': "230",
        'industrial': "150",
        'farmland': "230",
        'forest': "360",
        'orchard': "240",
    })
    hazard = rasterio.open(data_path / "hazard" / "inundation_map.tif")
    curves = pd.read_csv(data_path / "curves" / "curves.csv", dtype={"landuse": str}).set_index('Unnamed: 0')
    maximum_damage = pd.read_csv(data_path / "curves" / "maxdam.csv", dtype={"landuse": str}).set_index("landuse")

    objects['damage_rasterio'] = calculate_damage_per_object(
        objects=objects,
        hazard=hazard,
        curves=curves,
        maximum_damage=maximum_damage,
        column='landuse',
    )
    objects.to_file(output_folder / "damaged_objects_rasterio.gpkg", driver="GPKG")

    hazard = rxr.open_rasterio(data_path / "hazard" / "inundation_map.tif")
    objects['damage_xarray'] = calculate_damage_per_object(
        objects=objects,
        hazard=hazard,
        curves=curves,
        maximum_damage=maximum_damage,
        column='landuse',
    )
    objects.to_file(output_folder / "damaged_objects_xarray.gpkg", driver="GPKG")

    assert objects['damage_rasterio'].equals(objects['damage_xarray'])


def test_exact_extract():
    exposure_file = gpd.read_file(data_path / "landuse" / "landuse.shp")
    hazard_file = data_path / "hazard" / "inundation_map.tif"

    average = exact_extract(
        hazard_file, exposure_file, ["mean"], output="pandas", include_geom=True
    )
    average.to_file(output_folder / "inundation_map.gpkg")
