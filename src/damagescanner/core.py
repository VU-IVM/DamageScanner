"""DamageScanner - a directe damage assessment toolkit
"""

# DIRECT PHYSICAL DAMAGE MODULE

# Created by: Hans de Moel & Elco Koks
# Date: Januari 2013
#
# DESCRIPTION:
# Script of the DamageScanner in Python. The damagescanner calculates
# potential flood damages based on inundation depth and land use using
# so-called stage-damage curves. The DamageScanner was originally developed
# for the 'Netherlands Later' project (Klijn et al., 2007). 
# The original land-use classes were based on the Land-Use Scanner in order 
# to evaluate the effect of future land-use change on flood damages. The 
# stage-damage relations and the maxiumun damages associated with these
# land-use classes were derived from the HIS-SSM (Kok et al., 2005).
#

# Get all the needed modules
import os
import rasterio
import numpy
import pandas 
import geopandas
from tqdm import tqdm
from rasterio.mask import mask
from shapely.geometry import mapping
from rasterio.features import shapes

from damagescanner.vector import get_losses

def RasterScanner(landuse_map,inun_map,curve_path,maxdam_path,centimeters=False,save=False,**kwargs):
    """
    Raster-based implementation of a direct damage assessment.
    
    # INPUT PARAMETERS:
     landuse_map       Land-use map. Make sure the land-use categories
                       correspond with the curves and maximum damages (see
                       below). Furthermore, the resolution and extend of the
                       land-use map has to be exactly the same as the
                       inundation map
     inun_map          Map with inundation depth per grid cell. Make sure that
                       the unit of the inundation map corresponds with the unit of the first
                       column of the curves file
     curve_path        File with the stage-damage curves of the different
                       land-use classes.% 
     maxdam_path       Vector with the maximum damages per land-use class (in
                       euro/m2)

    # OUTPUT PARAMETERS:
     damagebin         Table with the land-use class numbers (1st) and the damage
                       for that land-use class (2nd) 
     damagemap         Map displaying the damage per grid cell of the area
    
    """      
        
    # load land-use map
    if isinstance(landuse_map,str):
        with rasterio.open(landuse_map) as src:
            landuse = src.read()[0,:,:]
            transform = src.transform
    else:
        landuse = landuse_map.copy()
    
    landuse_in = landuse.copy()
    
    # Load inundation map
    if isinstance(inun_map,str):
        with rasterio.open(inun_map) as src:
            inundation = src.read()[0,:,:]
    else:
        inundation = inun_map.copy()
        
    # check if land-use and inundation map have the same shape. 
    if landuse.shape != inundation.shape:
        print('ERROR: landuse and inundation maps are not the same shape. Fix this first')
        return None

    # set cellsize:
    if isinstance(landuse_map,str) & isinstance(inun_map,str):
        cellsize = src.res[0]*src.res[1]
    else:
        cellsize = kwargs['cellsize']

    # Load curves
    if isinstance(curve_path, pandas.DataFrame):
        curves = curve_path.values   
    elif isinstance(curve_path, numpy.ndarray):
        curves = curve_path
    elif curve_path.endswith('.csv'):
        curves = pandas.read_csv(curve_path).values

    #Load maximum damages
    if isinstance(maxdam_path, pandas.DataFrame):
        maxdam = maxdam_path.values #dict(zip(maxdam['landuse'],maxdam['damage']))
    elif isinstance(maxdam_path, numpy.ndarray):
        maxdam = maxdam_path
    elif maxdam_path.endswith('.csv'):
        maxdam = pandas.read_csv(maxdam_path,skiprows=1).values#dict(zip(pd.read_csv(maxdam_path)['landuse'],pd.read_csv(maxdam_path)['damage']))
        
    # Speed up calculation by only considering feasible points
    if centimeters:
        inundation[inundation>1000] = 0
    else:
        inundation[inundation>10] = 0
    inun = inundation * (inundation>=0) + 0
    inun[inun>=curves[:,0].max()] = curves[:,0].max()
    area = inun > 0
    waterdepth = inun[inun>0]
    landuse = landuse[inun>0]

    # Calculate damage per land-use class for structures
    numberofclasses = len(maxdam)
    alldamage = numpy.zeros(landuse.shape[0])
    damagebin = numpy.empty((numberofclasses, 2,)) * numpy.nan
    for i in range(0,numberofclasses):
        n = maxdam[i,0]
        damagebin[i,0] = n
        wd = waterdepth[landuse==n]
        alpha = numpy.interp(wd,((curves[:,0])),curves[:,i+1])
        damage = alpha*(maxdam[i,1]*cellsize)
        damagebin[i,1] = sum(damage)
        alldamage[landuse==n] = damage

    # create the damagemap
    damagemap = numpy.zeros((area.shape[0],area.shape[1]), dtype='int32')
    damagemap[area] = alldamage

    # create pandas dataframe with output
    loss_df = pandas.DataFrame(damagebin.astype(int),columns=['landuse','losses']).groupby('landuse').sum()
    
    if save:
        # requires adding output_path and scenario_name to function call
        if 'output_path' in kwargs:
            output_path = kwargs['output_path']
            if not os.path.exists(output_path):
                os.mkdir(output_path)
        if 'scenario_name' in kwargs:
            scenario_name = kwargs['scenario_name']
        if 'crs' in kwargs:
            crs = kwargs['crs']
        loss_df.to_csv(os.path.join(output_path,'{}_losses.csv'.format(scenario_name)))

        with rasterio.open(os.path.join(output_path,'{}_damagemap.tif'.format(scenario_name)), 'w', 
                                        driver='GTiff', height=damagemap.shape[0],
           width=damagemap.shape[1], count=1, dtype=damagemap.dtype, crs=crs, transform=transform,compress="LZW",) as dst:
                           dst.write(damagemap, 1)
                
    if 'in_millions' in kwargs:
        loss_df = loss_df/1e6
    
    # return output
    return loss_df,damagemap,landuse_in,inundation 


def VectorScanner(landuse,inun_file,curve_path,maxdam_path,landuse_col='landuse',centimeters=False,**kwargs):
    """
    Vector based implementation of a direct damage assessment
    """
    # load land-use map
    if isinstance(landuse,str):
        landuse = geopandas.from_file(landuse)
    elif isinstance(landuse, geopandas.GeoDataFrame):
        landuse = landuse.copy()
    elif isinstance(landuse, pandas.DataFrame):
        landuse = geopandas.GeoDataFrame(landuse,geometry='geometry')
    else:
        print('ERROR: landuse should either be a shapefile, a GeoDataFrame or a pandas Dataframe with a geometry column')
    
    if isinstance(inun_file,str):
         with rasterio.open(inun_file) as src:
            if src.crs.to_dict() != landuse.crs:
                landuse = landuse.to_crs(epsg=src.crs.to_epsg())
                geoms = [mapping(geom) for geom in landuse.geometry]
    
            out_image, out_transform = mask(src, geoms, crop=True)
            
            # if inundation map is not in centimeters (and assumed to be in meters), multiple by 100
            if not centimeters:
                out_image = out_image*100
            out_image[out_image > 1000] = -1
            out_image[out_image <= 0] = -1
    else:
        print('ERROR: inundation file should be a GeoTiff, or any other georeferenced format that Rasterio can read')

    # Load curves
    if isinstance(curve_path, pandas.DataFrame):
        curves = curve_path.copy()   
    elif isinstance(curve_path, numpy.ndarray):
       print('ERROR: for the vector-based approach we use a pandas DataFrame, not a Numpy Array')
    elif curve_path.endswith('.csv'):
        curves = pandas.read_csv(curve_path,index_col=[0])

    #Load maximum damages
    if isinstance(maxdam_path, pandas.DataFrame):
        maxdam = dict(zip(maxdam_path['landuse'],maxdam_path['damage']))
    elif isinstance(maxdam_path, numpy.ndarray):
        maxdam = dict(zip(maxdam_path[:,0],maxdam_path[:,1]))
    elif isinstance(maxdam_path, dict):
        maxdam = maxdam_path
    elif maxdam_path.endswith('.csv'):
        maxdam = dict(zip(pandas.read_csv(maxdam_path)['landuse'],pandas.read_csv(maxdam_path)['damage']))

    # convert raster to polygon
    results = (
        {'properties': {'raster_val': v}, 'geometry': s}
        for i, (s, v)
        in enumerate(
        shapes(out_image[0,:,:], mask=None, transform=out_transform)))

    gdf = geopandas.GeoDataFrame.from_features(list(results),crs=src.crs)
    gdf = gdf.loc[gdf.raster_val > 0]
    if centimeters:
        gdf = gdf.loc[gdf.raster_val < 1000]
    else:
        gdf = gdf.loc[gdf.raster_val < 10]
        
    # Split GeoDataFrame to make sure we have a unique shape per land use and inundation depth
    unique_df = []
    for row in tqdm(gdf.itertuples(index=False),total=len(gdf),desc='Get unique shapes'):
        matches = landuse.loc[list(landuse.sindex.intersection(row.geometry.bounds))]
        for match in matches.itertuples(index=False):
            if match.geometry.intersects(row.geometry):
                unique_df.append([row.raster_val,match.landuse,
                                 row.geometry.buffer(0).intersection(match.geometry)])

    # Create new dataframe
    new_gdf  = geopandas.GeoDataFrame(pandas.DataFrame(unique_df,columns=['depth',
                                                                          'landuse','geometry']))
    
    # And remove empty geometries where there was no intersection in the end
    new_gdf = new_gdf.loc[new_gdf.geometry.geom_type != 'GeometryCollection']

    # Get area of shape
    new_gdf['area_m2'] = new_gdf.area
    
    # And estimate the losses
    tqdm.pandas(desc='Estimate damages')
    new_gdf['damaged'] = new_gdf.progress_apply(lambda x : get_losses(x,curves,maxdam),axis=1)

    # And return the GeoDataFrame with all unique shapes
    return gdf.groupby(landuse_col).sum()[['area_m2','damaged']]