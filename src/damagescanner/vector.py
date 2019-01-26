import geopandas
import numpy 
import rasterio
from rasterio.mask import mask
from shapely.geometry import mapping
from rasterio.features import shapes


def intersect(x,landuse,lu_sindex,landuse_col):
    try:
        return landuse.loc[list(lu_sindex.intersection(x.geometry.bounds))][landuse_col].values[0]
    except:
        return None

def get_losses(x,damage_curves,damage_values):
    return np.interp(x.raster_val/100,list(damage_curves.index),list(damage_curves[x.landuse]))*damage_values[x.landuse]
    
def get_exposure_losses(inun_file,landuse,lu_sindex,damage_curves,damage_values):

    geoms = [mapping(geom) for geom in landuse.geometry]

    with rasterio.open(inun_file) as src:
        out_image, out_transform = mask(src, geoms, crop=True)
        out_image[out_image == 999] = -1
        out_image[out_image <= 0] = -1

        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v)
            in enumerate(
            shapes(out_image[0,:,:], mask=None, transform=out_transform)))

        gdf = geopandas.GeoDataFrame.from_features(list(results),crs='epsg:28992')
        gdf = gdf.loc[gdf.raster_val > 0]
        gdf = gdf.loc[gdf.raster_val < 5000]

        gdf['landuse'] = gdf.apply(lambda x : intersect(x,landuse,lu_sindex),axis=1)
        gdf['area_m2'] = gdf.area
        
        gdf['damaged'] = gdf.apply(lambda x : get_losses(x,damage_curves,damage_values),axis=1)
        
        return gdf.groupby('landuse').sum()