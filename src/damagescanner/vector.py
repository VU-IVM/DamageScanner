import geopandas
import pandas
import ogr
import numpy 
from tqdm import tqdm
from shapely.wkt import loads

def fetch_landuse(osm_path):
    """
    Function to extract land-use polygons from OpenStreetMap
    
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.
        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique land-use polygons.
    
    """
    
    driver=ogr.GetDriverByName('OSM')
    data = driver.Open(osm_path)
                       
    sql_lyr = data.ExecuteSQL("SELECT osm_id,landuse from multipolygons where landuse is not null")

    osm_data = []
    if data is not None:
        for feature in sql_lyr:
            try:
                if feature.GetField('landuse') is not None:
                    shapely_geo = loads(feature.geometry().ExportToWkt()) 
                    if shapely_geo is None:
                        continue
                    data_type=feature.GetField('landuse')
                    osm_id = feature.GetField('osm_id')
    
                    osm_data.append([osm_id,data_type,shapely_geo])
            except:
                print("WARNING: skipped landuse shape")                       
    else:
        print("ERROR: Nonetype error when requesting SQL. Check required.")    

    if len(osm_data) > 0:
        return geopandas.GeoDataFrame(osm_data,columns=['osm_id','landuse','geometry'],
                                crs={'init': 'epsg:4326'})
    else:
        print("WARNING: No landuse shapes or No Memory. returning empty GeoDataFrame") 
        return geopandas.GeoDataFrame(columns=['osm_id','landuse','geometry'],crs={'init': 'epsg:4326'})

    
def fetch_buildings(osm_path):
    """
    Function to extract building polygons from OpenStreetMap
    
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.
        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique building polygons.
    
    """
    driver=ogr.GetDriverByName('OSM')
    data = driver.Open(osm_path)

    buildings=[]    
    if data is not None:
        sql_lyr = data.ExecuteSQL("SELECT osm_id,amenity,building from multipolygons where building is not null")
        for feature in sql_lyr:
            try:
                if feature.GetField('building') is not None:
                    osm_id = feature.GetField('osm_id')
                    shapely_geo = loads(feature.geometry().ExportToWkt()) 
                    if shapely_geo is None:
                        continue
                    building=feature.GetField('building')
                    amenity=feature.GetField('amenity')

                    buildings.append([osm_id,building,amenity,shapely_geo])
            except:
                    print("WARNING: skipped building")
    else:
        print("ERROR: Nonetype error when requesting SQL. Check required.")    

    if len(buildings) > 0:
        return geopandas.GeoDataFrame(buildings,columns=['osm_id','building','amenity','geometry'],crs={'init': 'epsg:4326'})
    else:
        print("WARNING: No buildings or No Memory. returning empty GeoDataFrame") 
        return geopandas.GeoDataFrame(columns=['osm_id','building','amenity','geometry'],crs={'init': 'epsg:4326'})


def fetch_roads(osm_path):
    """
    Function to extract road linestrings from OpenStreetMap
    
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.
        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique road linestrings.
    
    """
    
    driver=ogr.GetDriverByName('OSM')
    data = driver.Open(osm_path)

    roads=[]    
    if data is not None:
        sql_lyr = data.ExecuteSQL("SELECT osm_id,highway FROM lines WHERE highway IS NOT NULL")
        for feature in sql_lyr:
            try:
                if feature.GetField('highway') is not None:
                    osm_id = feature.GetField('osm_id')
                    shapely_geo = loads(feature.geometry().ExportToWkt()) 
                    if shapely_geo is None:
                        continue
                    highway=feature.GetField('highway')
                    roads.append([osm_id,highway,shapely_geo])
            except:
                    print("WARNING: skipped a road")
    else:
        print("ERROR: Nonetype error when requesting SQL. Check required.")    

    if len(roads) > 0:
        return geopandas.GeoDataFrame(roads,columns=['osm_id','highway','geometry'],crs={'init': 'epsg:4326'})
    else:
        print("WARNING: No buildings or No Memory. returning empty GeoDataFrame") 
        return geopandas.GeoDataFrame(columns=['osm_id','highway','geometry'],crs={'init': 'epsg:4326'})
    
def fetch_railway(osm_path):
    """
    Function to extract railway linestrings from OpenStreetMap
    
    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.
        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique land-use polygons.
    
    """
    
    driver=ogr.GetDriverByName('OSM')
    data = driver.Open(osm_path)

    railways=[]    
    if data is not None:
        sql_lyr = data.ExecuteSQL("SELECT osm_id,service,railway FROM lines WHERE railway IS NOT NULL")
        for feature in sql_lyr:
            try:
                if feature.GetField('railway') is not None:
                    osm_id = feature.GetField('osm_id')
                    shapely_geo = loads(feature.geometry().ExportToWkt()) 
                    if shapely_geo is None:
                        continue
                    railway=feature.GetField('railway')
                    railways.append([osm_id,railway,shapely_geo])
            except:
                    print("warning: skipped building")
    else:
        print("ERROR: Nonetype error when requesting SQL. Check required.")    

    if len(railways) > 0:
        return geopandas.GeoDataFrame(railways,columns=['osm_id','building','geometry'],crs={'init': 'epsg:4326'})
    else:
        print("WARNING: No buildings or No Memory. returning empty GeoDataFrame") 
        return geopandas.GeoDataFrame(columns=['osm_id','building','geometry'],crs={'init': 'epsg:4326'})

def remove_overlap_openstreetmap(gdf):
    """
    Function to remove overlap in polygons in from OpenStreetMap.
    
    Arguments:
        *gdf* : a geopandas GeoDataFrame with all unique railway linestrings.
        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with (almost) non-overlapping polygons.
    
    """
    
    gdf['sq_area'] = gdf.area

    new_landuse = []
    for use in tqdm(gdf.itertuples(index=False),total=len(gdf),desc='Get unique shapes'):
        use_geom = use.geometry
        matches = gdf.loc[list(gdf.sindex.intersection(use.geometry.bounds))]
        for match in matches.itertuples(index=False):
            if use.sq_area > match.sq_area:
                use_geom = use_geom.difference(match.geometry)
        new_landuse.append([use.osm_id,use.landuse,use_geom])

    new_gdf  =  geopandas.GeoDataFrame(pandas.DataFrame(new_landuse,columns=['osm_id','landuse','geometry'])) 
    new_gdf.crs = {'init' : 'epsg:4326'}
    return new_gdf


def extract_value_other_gdf(x,gdf,col_name):
    """
    Function to extract value from column from other GeoDataFrame
    
    Arguments:
        *x* : row of main GeoDataFrame.
        
        *gdf* : geopandas GeoDataFrame from which we want to extract values.
        
        *col_name* : the column name from which we want to get the value.
        
    
    """
    try:
        return gdf.loc[list(gdf.sindex.intersection(x.geometry.bounds))][col_name].values[0]
    except:
        return None

def get_losses(x,damage_curves,damage_values):
    """
    Function to estimate the damages.
    
    Arguments:
        *x* : row of main GeoDataFrame
        
        *damage_curves*: pandas DataFrame of curves. Inundation depths should be the index.
        
        *damage_values*: dictionary with maximum damage values.
        
    Returns:
        
        Total damage for the given land-use object.
    
    """
    
    return numpy.interp(x.depth,list(damage_curves.index),list(damage_curves[x.landuse]))*damage_values[x.landuse]
    
