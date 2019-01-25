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
import rasterio
import numpy
import pandas 

def RasterScanner(landuse_map,inun_map,curve_path,maxdam_path,cellsize=100,save=False,**kwargs):
    """
    Raster-based implementation of a direct damage assessment.
    
    # INPUT PARAMETERS:
    # LANDUSE           Land-use map. Make sure the land-use categories
    #                   correspond with the curves and maximum damages (see
    #                   below). Furthermore, the resolution and extend of the
    #                   land-use map has to be exactly the same as the
    #                   inundation map
    # INUNDATION        Map with inundation depth per grid cell. Make sure that
    #                   the unit of the inundation map corresponds with the unit of the first
    #                   column of the curves file
    # CURVES            File with the stage-damage curves of the different
    #                   land-use classes.% 
    # MAXDAM            Vector with the maximum damages per land-use class (in
    #                   euro/m2)

    # OUTPUT PARAMETERS:
    # DAMAGEBIN         Table with the land-use class numbers (1st), the damage
    #                   for that land-use class (2nd) and the amount of cells of that given
    #                   land-use class in the inundated area (3rd)
    # DAMAGEMAP         Map displaying the damage per grid cell of the area
    #
    """   
        
    # load land-use map
    if landuse_map.endswith('.tif'):
        with rasterio.open(landuse_map) as src:
            landuse = src.read()[0,:,:]
    else:
        landuse = landuse_map.copy()
    
    # Load inundation map
    if inun_map.endswith('.tif'):
        with rasterio.open(inun_map) as src:
            inundation = src.read()[0,:,:]
    else:
        inundation = inun_map.copy()
        
    # check if land-use and inundation map have the same shape. 
    if landuse.shape != inundation.shape:
        print('ERROR: landuse and inundation maps are not the same shape. Fix this first')
        return None

    # Load curves
    if isinstance(curve_path, pandas.DataFrame):
        curves = curve_path.values    
    elif curve_path.endswith('.csv'):
        curves = pandas.read_csv(curve_path,index_col=0).values

    #Load maximum damages
    if isinstance(maxdam_path, pandas.DataFrame):
        maxdam = maxdam_path.values #dict(zip(maxdam['landuse'],maxdam['damage']))
    elif maxdam_path.endswith('.csv'):
        maxdam = pandas.read_csv(maxdam_path,skiprows=1).values#dict(zip(pd.read_csv(maxdam_path)['landuse'],pd.read_csv(maxdam_path)['damage']))
    
    
    # Speed up calculation by only considering feasible points
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
    
    # return output
    return loss_df,damagemap

def VectorScanner(LANDUSE,INUNDATION,CURVES,MAXDAM,CELLSIZE,BASEMAP,SCEN_FLOOD,ADAPTATION,RETURN_PERIOD):
    """
    Vector based implementation of a direct damage assessment
    """


