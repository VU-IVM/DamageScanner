import shapely
import geopandas
import pandas
from osgeo import ogr,gdal
import os
import numpy as np
from tqdm import tqdm
import pyproj


def query_b(geo_type, key_col, **val_constraint):
    """
    Builds an SQL query from the values passed to the function.
    
    Arguments:
            geo_type : Type of geometry (osm layer) to search for.
            key_col : A list of keys/columns that should be selected from the layer.
            val_constraint : A dictionary of constraints for the values. e.g. WHERE 'value'>20 or 'value'='constraint'
    Returns:
        string: A SQL query string.
    """
    query = "SELECT osm_id"
    
    for column in key_col:
        query += ", " + column

    query += " FROM " + geo_type + " WHERE "

    if val_constraint:
        constraint_clauses = []
        for key, values in val_constraint.items():
            for value in values:
                constraint_clauses.append(f"{key} {value}")
        query += " AND ".join(constraint_clauses) + " AND "
        
    query += f"{key_col[0]} IS NOT NULL"
    return query


def retrieve(osm_path, geo_type, key_col, **val_constraint):
    """
    Function to extract specified geometry and keys/values from OpenStreetMap.

    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.     
        *geo_type* : Type of Geometry to retrieve. e.g. lines, multipolygons, etc.
        *key_col* : These keys will be returned as columns in the dataframe.
        ***val_constraint: A dictionary specifying the value constraints.  
        A key can have multiple values (as a list) for more than one constraint for key/value.  
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all columns, geometries, and constraints specified.    
    """
    driver = ogr.GetDriverByName('OSM')
    data = driver.Open(osm_path)
    query = query_b(geo_type, key_col, **val_constraint)
    sql_lyr = data.ExecuteSQL(query)
    features = []

    # cl = columns 
    cl = ['osm_id'] 
    for a in key_col: cl.append(a)
    if data is not None:
        for feature in sql_lyr:
            try:
                if feature.GetField(key_col[0]) is not None:
                    shapely_geo = shapely.from_wkt(feature.geometry().ExportToWkt()) 
                    if shapely_geo is None:
                        continue
                    # field will become a row in the dataframe.
                    field = []
                    for i in cl: field.append(feature.GetField(i))
                    field.append(shapely_geo)   
                    features.append(field)
            except:
                print("WARNING: skipped OSM feature")   
    else:
        print("ERROR: Nonetype error when requesting SQL. Check required.")    
    cl.append('geometry')                   
    if len(features) > 0:
        return geopandas.GeoDataFrame(features, columns=cl, crs={'init': 'epsg:4326'})
    else:
        print("WARNING: No features or No Memory. returning empty GeoDataFrame") 
        return geopandas.GeoDataFrame(columns=['osm_id','geometry'], crs={'init': 'epsg:4326'})

def landuse(osm_path):
    """
    Function to extract land-use polygons from OpenStreetMap    

    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique land-use polygons.    
    """    
    return(retrieve(osm_path,'multipolygons',['landuse']))

def buildings(osm_path):
    """
    Function to extract building polygons from OpenStreetMap    

    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique building polygons.    
    """
    return retrieve(osm_path, 'multipolygons',['building','amenity'])

def roads(osm_path):
    """
    Function to extract road linestrings from OpenStreetMap  

    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique road linestrings.
    """   
    return retrieve(osm_path,'lines',['highway']) 
 
def railway(osm_path):
    """
    Function to extract railway linestrings from OpenStreetMap

    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.       
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique land-use polygons.
    """ 
    return retrieve(osm_path,'lines',['railway','service'],**{"service":[" IS NOT NULL"]})

def ferries(osm_path):
    """
    Function to extract road linestrings from OpenStreetMap

    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique road linestrings.
    """
    return retrieve(osm_path,'lines',['route'],**{"route":["='ferry'",]})

def electricity(osm_path):
    """
    Function to extract railway linestrings from OpenStreetMap

    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique land-use polygons.   
    """    
    return retrieve(osm_path,'lines',['power','voltage'],**{'voltage':[" IS NULL"],})

def mainRoads(osm_path):
    """
    Function to extract main road linestrings from OpenStreetMap  

    Arguments:
        *osm_path* : file path to the .osm.pbf file of the region 
        for which we want to do the analysis.        
    Returns:
        *GeoDataFrame* : a geopandas GeoDataFrame with all unique main road linestrings.   
    """ 
    return retrieve(osm_path,'lines',['highway','oneway','lanes','maxspeed'],**{'highway':["='primary' or ","='trunk' or ","='motorway' or ","='trunk_link' or ",
                    "='primary_link' or ", "='secondary' or ","='tertiary' or ","='tertiary_link'"]})


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

def reproject(df_ds,current_crs="epsg:4326",approximate_crs = "epsg:3035"):
    """
    Function to reproject a GeoDataFrame to a different CRS.

    Args:
        df_ds (pandas DataFrame): _description_
        current_crs (str, optional): _description_. Defaults to "epsg:4326".
        approximate_crs (str, optional): _description_. Defaults to "epsg:3035".

    Returns:
        _type_: _description_
    """    
    geometries = df_ds['geometry']
    coords = shapely.get_coordinates(geometries)
    transformer=pyproj.Transformer.from_crs(current_crs, approximate_crs,always_xy=True)
    new_coords = transformer.transform(coords[:, 0], coords[:, 1])
    
    return shapely.set_coordinates(geometries.copy(), np.array(new_coords).T) 
    
def get_damage_per_object(obj,df_ds,objects,curves,maxdam):
    """"
    Function to calculate the damage per object.

    Args:
        obj (tuple): The object for which we want to calculate the damage.
        df_ds (pandas DataFrame): The dataframe with the hazard points.
        objects (pandas DataFrame): The dataframe with the objects.
        curves (pandas DataFrame): The dataframe with the vulnerability curves.
        maxdam (dictionary): The dictionary with the maximum damage per object type.

    Returns:
        tuple: The damage per object.
    """  

    # find the exact hazard overlays:
    get_hazard_points = df_ds.iloc[obj[1]['hazard_point'].values].reset_index()
    get_hazard_points = get_hazard_points.loc[shapely.intersects(get_hazard_points.geometry.values,objects.iloc[obj[0]].geometry)]

    # get the object type and the object geometry
    object_type = objects.iloc[obj[0]].obj_type
    object_geom = objects.iloc[obj[0]].geometry
    
    # get the maximum damage for the object type
    maxdam_object = maxdam[object_type]

    # get the vulnerability curves for the object type
    hazard_intensity = curves[object_type].index.values
    fragility_values = curves[object_type].values
                
    # if there are no hazard points, return 0 damage
    if len(get_hazard_points) == 0:
        return obj[0],0
    else:
        
        # run the analysis for lines (object_type = 1, e.g. railway lines)
        if shapely.get_type_id(object_geom) == 1:
            get_hazard_points['overlay_meters'] = shapely.length(shapely.intersection(get_hazard_points.geometry.values,object_geom))
            return obj[0],np.sum((np.interp(get_hazard_points.haz_val.values,hazard_intensity,fragility_values))*get_hazard_points.overlay_meters*maxdam_object)
        
        # run the analysis for polygons (object_type = 2, e.g. buildings)
        elif (shapely.get_type_id(object_geom) == 3) | (shapely.get_type_id(object_geom) == 6):
            get_hazard_points['overlay_m2'] = shapely.area(shapely.intersection(get_hazard_points.geometry.values,object_geom))
            return obj[0],get_hazard_points.apply(lambda x: np.interp(x.haz_val, hazard_intensity, fragility_values)*maxdam_object*x.overlay_m2,axis=1).sum()     
        
        # run the analysis for points (object_type = 0, e.g. hospitals)
        else:
            print(shapely.get_type_id(object_geom))
            return obj[0],np.sum((np.interp(get_hazard_points.haz_val.values,hazard_intensity,fragility_values))*maxdam_object)