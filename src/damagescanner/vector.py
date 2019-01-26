import geopandas
#import ogr
import numpy 

from shapely.wkt import loads

def fetch_landuse(osm_path):
#    driver=ogr.GetDriverByName('OSM')
    data = driver.Open(osm_path)
                       
    sql_lyr = data.ExecuteSQL("SELECT osm_id,landuse from multipolygons where landuse is not null")

    osm_data = []
    for feature in sql_lyr:
            if feature.GetField('landuse') is not None:
                shapely_geo = loads(feature.geometry().ExportToWkt()) 
                if shapely_geo is None:
                    continue
                data_type=feature.GetField('landuse')
                osm_id = feature.GetField('osm_id')

                osm_data.append([osm_id,data_type,shapely_geo])
                       
    return geopandas.GeoDataFrame(osm_data,columns=['osm_id','landuse','geometry'],
                            crs={'init': 'epsg:4326'})
    

def intersect(x,landuse,lu_sindex,landuse_col):
    try:
        return landuse.loc[list(lu_sindex.intersection(x.geometry.bounds))][landuse_col].values[0]
    except:
        return None

def get_losses(x,damage_curves,damage_values,**kwargs):
    return numpy.interp(x.raster_val/100,list(damage_curves.index),list(damage_curves[x.landuse]))*damage_values[x.landuse]
    
