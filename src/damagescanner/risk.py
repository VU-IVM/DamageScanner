# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 11:54:09 2019

@author: cenv0574
"""
import os
import pandas
from scipy import integrate
from tqdm import tqdm

from damagescanner.core import RasterScanner,VectorScanner

def monetary_risk(RPS,loss_list):
    """
    Calculates the monetary risk based on the return periods and the losses.

    Arguments:
        *RPS* : List of return periods (in years) for which the losses are calculated.
        *loss_list* : List of losses (in euro) for each return period.

    Returns:
        *total_risk* : Returns the total risk for the area
    """
    
    return integrate.simps(y=loss_list[::-1], x=RPS[::-1])

def RasterBased(landuse_ras,inundation_path,curve_path,maxdam_path,per_landuse=False):
    """
    Raster-based implementation of a risk assessment.
    
    Arguments:
        *landuse_map* : GeoTiff with land-use information per grid cell. Make sure 
        the land-use categories correspond with the curves and maximum damages 
        (see below). Furthermore, the resolution and extend of the land-use map 
        has to be exactly the same as the inundation map.
     
        *inun_map* : GeoTiff with inundation depth per grid cell. Make sure 
        that the unit of the inundation map corresponds with the unit of the 
        first column of the curves file.
     
        *curve_path* : File with the stage-damage curves of the different 
        land-use classes. Can also be a pandas DataFrame or numpy Array.
     
        *maxdam_path* : File with the maximum damages per land-use class 
        (in euro/m2). Can also be a pandas DataFrame or numpy Array.

    Optional Arguments:
        *per_landuse* : Set to **True** if you would like the output er land-use class.

    Returns:
        *total_risk* : Returns the total risk for the area
        
    """    
    
    # Get list of inundation maps and make sure they are in ascending order
    list_inundation_maps = os.listdir(inundation_path)
    list_inundation_maps.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    list_inundation_maps = [os.path.join(inundation_path,x) for x in list_inundation_maps]
    
    RPS = [1/int(''.join(filter(str.isdigit, x))) for x in list_inundation_maps]
    return_periods =  ['rp_'+''.join(filter(str.isdigit, x)) for x in list_inundation_maps]

    # loop through the inundation maps and calculate the losses for each map
    damage_list = {}
    for inun_map in tqdm(list_inundation_maps,total=len(list_inundation_maps),desc='Return Periods'):
        return_period = 'rp_'+''.join(filter(str.isdigit, inun_map))
        damage_list[return_period] = RasterScanner(landuse_ras,inun_map,curve_path,maxdam_path)[0]

    damage_table = pandas.concat(damage_list).unstack(level=0)
    damage_table.columns = damage_table.columns.droplevel(0)
    damage_table = damage_table[return_periods]

    # calculate the total risk and return it
    if per_landuse:
        return pandas.DataFrame(damage_table.apply(lambda x: monetary_risk(RPS,list(x)),axis=1),columns=['tot_risk'])
    else:    
        return monetary_risk(RPS,list(damage_table.sum(axis=0)))/1e6


def VectorBased(landuse_vec,inundation_path,curve_path,maxdam_path,per_landuse=False):
    """
    Vector based implementation of a direct damage assessment
    
    Arguments:
        *landuse_map* : Shapefile, Pandas DataFrame or Geopandas GeoDataFrame 
        with land-use information of the area.
     
        *inun_map* : GeoTiff with inundation depth per grid cell. Make sure 
        that the unit of the inundation map corresponds with the unit of the 
        first column of the curves file. 
     
        *curve_path* : File with the stage-damage curves of the different 
        land-use classes. Can also be a pandas DataFrame (but not a numpy Array).
     
        *maxdam_path* : File with the maximum damages per land-use class 
        (in euro/m2). Can also be a pandas DataFrame (but not a numpy Array).

    Optional Arguments:
        *per_landuse* : Set to **True** if you would like the output er land-use class.
        
        *landuse_col* : Specify the column name of the unique landuse id's. 
        Default is set to **landuse**.
        
    Returns:
        *total_risk* : Returns the total risk for the area
    """
     
    # Get list of inundation maps and make sure they are in ascending order
    list_inundation_maps = os.listdir(inundation_path)
    list_inundation_maps.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    list_inundation_maps = [os.path.join(inundation_path,x) for x in list_inundation_maps]
    
    RPS = [1/int(''.join(filter(str.isdigit, x))) for x in list_inundation_maps]
    return_periods =  ['rp_'+''.join(filter(str.isdigit, x)) for x in list_inundation_maps]

    # loop through the inundation maps and calculate the losses for each map
    damage_list = {}
    for inun_map in tqdm(list_inundation_maps,total=len(list_inundation_maps),desc='Return Periods'):
        return_period = 'rp_'+''.join(filter(str.isdigit, inun_map))
        damage_list[return_period] =  pandas.DataFrame(
            VectorScanner(landuse_vec,inun_map,curve_path,maxdam_path,print_tqdm=False).groupby('landuse')['tot_dam'].sum())      
        
    damage_table = pandas.concat(damage_list).unstack(level=0)
    damage_table.columns = damage_table.columns.droplevel(0)
    damage_table = damage_table[return_periods]

    # calculate the total risk and return it
    if per_landuse:
        return pandas.DataFrame(damage_table.apply(lambda x: monetary_risk(RPS,list(x)),axis=1),columns=['tot_risk'])
    else:    
        return monetary_risk(RPS,list(damage_table.sum(axis=0)))/1e6