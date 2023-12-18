import shapely
import geopandas as gpd
import pandas
import numpy as np
from tqdm import tqdm
import pyproj


def match_raster_to_vector(hazard,landuse,lu_crs,haz_crs,resolution,hazard_col):
        
    """ Matches the resolution and extent of a raster to a vector file.

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
    hazard = hazard.rename({'Latitude': 'y','Longitude': 'x'})

    # Set the x and y dimensions in the hazard variable to 'x' and 'y' respectively
    hazard.rio.set_spatial_dims(x_dim="x",y_dim="y", inplace=True)

    # Set the crs of the landuse variable to lu_crs
    landuse.rio.write_crs(lu_crs,inplace=True)

    # Reproject the landuse variable from EPSG:4326 to EPSG:3857
    landuse = landuse.rio.reproject("EPSG:3857",resolution=resolution)

    # Get the minimum longitude and latitude values in the landuse variable
    min_lon = landuse.x.min().to_dict()['data']
    min_lat = landuse.y.min().to_dict()['data']

    # Get the maximum longitude and latitude values in the landuse variable
    max_lon = landuse.x.max().to_dict()['data']
    max_lat = landuse.y.max().to_dict()['data']

    # Create a bounding box using the minimum and maximum latitude and longitude values
    area = gpd.GeoDataFrame([shapely.box(min_lon,min_lat,max_lon, max_lat)],columns=['geometry'])

    # Set the crs of the bounding box to EPSG:3857
    area.crs = 'epsg:3857'

    # Convert the crs of the bounding box to EPSG:4326
    area = area.to_crs('epsg:4326')

    # Clip the hazard variable to the extent of the bounding box
    hazard = hazard.rio.clip(area.geometry.values, area.crs)

    # Reproject the hazard variable to EPSG:3857 with the desired resolution
    hazard = hazard.rio.reproject("EPSG:3857",resolution=resolution)

    # Clip the hazard variable again to the extent of the bounding box
    hazard = hazard.rio.clip(area.geometry.values, area.crs)

    # If the hazard variable has fewer columns and rows than the landuse variable, reproject 
    # the landuse variable to match the hazard variable
    if (len(hazard.x)<len(landuse.x)) & (len(hazard.y)<len(landuse.y)):
        landuse= landuse.rio.reproject_match(hazard)

    # If the hazard variable has more columns and rows than the landuse variable, 
    # reproject the hazard variable to match the landuse variable

    elif (len(hazard.x)>len(landuse.x)) & (len(hazard.y)>len(landuse.y)):
        hazard = hazard.rio.reproject_match(landuse)

    # Convert the hazard and landuse variable to a numpy array
    landuse = landuse['band_data'].to_numpy()[0,:,:]
    hazard = hazard[hazard_col].to_numpy()[0,:,:]

    return hazard,landuse

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

    new_gdf  =  gpd.GeoDataFrame(pandas.DataFrame(new_landuse,columns=['osm_id','landuse','geometry'])) 
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