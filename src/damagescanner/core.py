"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2019 Elco Koks. All versions released under the MIT license.
"""

# Get all the needed modules
import os
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
from affine import Affine
from tqdm import tqdm
from rasterio.mask import mask
from shapely.geometry import mapping
from rasterio.features import shapes
from os.path import join as p_join
import warnings

from damagescanner.vector import get_losses
from damagescanner.raster import match


def check_output_path(given_args):
    """Ensures given output path exists.

    Arguments:
        *given_args* : dict, of keyword arguments.

    Returns:
        *str* : output_path, which may be empty string ('')
    """
    output_path = given_args.get('output_path', '')
    if output_path != '' and not os.path.exists(output_path):
        os.mkdir(output_path)
    return output_path


def check_scenario_name(given_args):
    """Ensures given output path exists.

    Arguments:
        *given_args* : dict, of keyword arguments.

    Returns:
        *str* : scenario_name
    """
    scenario_name = given_args.get('scenario_name', False)
    if not scenario_name:
        raise ValueError("Required `scenario_name` not defined.")

    return scenario_name


def RasterScanner(landuse_map,
                  inun_map,
                  curve_path,
                  maxdam_path,
                  dtype = np.int32,
                  save=False,
                  **kwargs):
    """
    Raster-based implementation of a direct damage assessment.
    
    Arguments:
        *landuse_map* : GeoTiff with land-use information per grid cell. Make sure 
        the land-use categories correspond with the curves and maximum damages 
        (see below). Furthermore, the resolution and extend of the land-use map 
        has to be exactly the same as the inundation map.
     
        *inun_map* : GeoTiff with inundation depth per grid cell. Make sure 
        that the unit of the inundation map corresponds with the unit of the 
        first column of the curves file.
     
        *curve_path* : File with the stage-damage curves of the different 
        land-use classes. Values should be given as ratios, i.e. between 0 and 1.
        Can also be a pandas DataFrame or numpy Array.
     
        *maxdam_path* : File with the maximum damages per land-use class 
        (in euro/m2). Can also be a pandas DataFrame or numpy Array.
                        
        *dtype*: Set the dtype to the requires precision. This will affect the output damage raster as well 
     
    Optional Arguments:
        *save* : Set to True if you would like to save the output. Requires 
        several **kwargs**
        
    kwargs:
        *nan_value* : if nan_value is provided, will mask the inundation file. 
        This option can significantly fasten computations
        
        *cell_size* : If both the landuse and inundation map are numpy arrays, 
        manually set the cell size.
        
        *output_path* : Specify where files should be saved.
        
        *scenario_name*: Give a unique name for the files that are going to be saved.
        
        *in_millions*: Set to True if all values should be set in millions.
        
        *crs*: Specify crs if you only read in two numpy array
        
        *transform*: Specify transform if you only read in numpy arrays in order to save the result raster
        
    Raises:
        *ValueError* : on missing kwarg options

    Returns:    
     *damagebin* : Table with the land-use class numbers (1st column) and the 
     damage for that land-use class (2nd column).
     
     *damagemap* : Map displaying the damage per grid cell of the area.
     
    """
    # load land-use map
    if isinstance(landuse_map, str):
        with rasterio.open(landuse_map) as src:
            landuse = src.read()[0, :, :]
            transform = src.transform
    else:
        landuse = landuse_map.copy()

    landuse_in = landuse.copy()

    # Load inundation map
    if isinstance(inun_map, str):
        with rasterio.open(inun_map) as src:
            inundation = src.read()[0, :, :]
            transform = src.transform
    else:
        inundation = inun_map.copy()

    # check if land-use and inundation map have the same shape.
    if landuse.shape != inundation.shape:
        warnings.warn(
            "WARNING: landuse and inundation maps are not the same shape. Let's fix this first!"
        )
        landuse, inundation, intersection = match(landuse_map, inun_map)

        # create the right affine for saving the output
        transform = Affine(transform[0], transform[1], intersection[0],
                           transform[3], transform[4], intersection[1])

    # set cellsize:
    if isinstance(landuse_map, str) | isinstance(inun_map, str):
        cellsize = src.res[0] * src.res[1]
    else:
        try:
            cellsize = kwargs['cellsize']
        except KeyError:
            raise ValueError("Required `cellsize` not given.")

    # Load curves
    if isinstance(curve_path, pd.DataFrame):
        curves = curve_path.values
    elif isinstance(curve_path, np.ndarray):
        curves = curve_path
    elif curve_path.endswith('.csv'):
        curves = pd.read_csv(curve_path).values
    
    if ((curves>1).all()) or ((curves<0).all()):
        raise ValueError("Stage-damage curve values must be between 0 and 1")

    # Load maximum damages
    if isinstance(maxdam_path, pd.DataFrame):
        maxdam = maxdam_path.values
    elif isinstance(maxdam_path, np.ndarray):
        maxdam = maxdam_path
    elif maxdam_path.endswith('.csv'):
        maxdam = pd.read_csv(maxdam_path).values
    
    if maxdam.shape[0] != (curves.shape[1]-1):
        raise ValueError("Dimensions between maximum damages and the number of depth-damage curve do not agree")
  
    # Speed up calculation by only considering feasible points
    if kwargs.get('nan_value'):
        nan_value = kwargs.get('nan_value')
        inundation[inundation == nan_value] = 0

    inun = inundation * (inundation >= 0) + 0   
    inun[inun >= curves[:, 0].max()] = curves[:, 0].max()   
    area = inun > 0
    waterdepth = inun[inun > 0]
    landuse = landuse[inun > 0]

    # Calculate damage per land-use class for structures
    numberofclasses = len(maxdam)
    alldamage = np.zeros(landuse.shape[0])
    damagebin = np.zeros((
        numberofclasses,
        2,
    ))
    for i in range(0, numberofclasses):
        n = maxdam[i, 0]
        damagebin[i, 0] = n
        wd = waterdepth[landuse == n]
        alpha = np.interp(wd, ((curves[:, 0])), curves[:, i + 1])
        damage = alpha * (maxdam[i, 1] * cellsize)
        damagebin[i, 1] = sum(damage)
        alldamage[landuse == n] = damage
    
    # create the damagemap
    damagemap = np.zeros((area.shape[0], area.shape[1]), dtype=dtype)
    damagemap[area] = alldamage

    # create pandas dataframe with output
    loss_df = pd.DataFrame(damagebin.astype(dtype),
                           columns=['landuse',
                                    'losses']).groupby('landuse').sum()

    if save:
        crs = kwargs.get('crs', src.crs)
        transform = kwargs.get('transform', transform)

        # requires adding output_path and scenario_name to function call
        # If output path is not defined, will place file in current directory
        output_path = check_output_path(kwargs)
        scenario_name = check_scenario_name(kwargs)
        path_prefix = p_join(output_path, scenario_name)

        loss_fn = '{}_losses.csv'.format(path_prefix)
        loss_df.to_csv(loss_fn)

        dmap_fn = '{}_damagemap.tif'.format(path_prefix)
        rst_opts = {
            'driver': 'GTiff',
            'height': damagemap.shape[0],
            'width': damagemap.shape[1],
            'count': 1,
            'dtype': dtype,
            'crs': crs,
            'transform': transform,
            'compress': "LZW"
        }
        with rasterio.open(dmap_fn, 'w', **rst_opts) as dst:
            dst.write(damagemap, 1)

    if 'in_millions' in kwargs:
        loss_df = loss_df / 1e6

    # return output
    return loss_df, damagemap, landuse_in, inundation


def VectorScanner(landuse,
                  inun_file,
                  curve_path,
                  maxdam_path,
                  landuse_col='landuse',
                  inun_col='inun_val',
                  centimeters=False,
                  save=False,
                  **kwargs):
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
        *centimeters* : Set to True if the inundation map and curves are in 
        centimeters
        
        *landuse_col* : Specify the column name of the unique landuse id's. 
        Default is set to **landuse**.
        
        *inun_col* : Specify the column name of the inundation depth 
        Default is set to **inun_val**.       
        
        *save* : Set to True if you would like to save the output. Requires 
        several **kwargs**
        
    kwargs:
        *output_path* : Specify where files should be saved.
        
        *scenario_name*: Give a unique name for the files that are going to be saved.
        
        *print_tqdm*: Set to **False** when progress output is undesired.
    
    Raises:
        *ValueError* : on missing kwargs
    
    Returns:    
     *damagebin* : Table with the land-use class names (1st column) and the 
     damage for that land-use class (2nd column).
     
    """
    # load land-use map
    if isinstance(landuse, str):
        landuse = gpd.read_file(landuse)
    elif isinstance(landuse, gpd.GeoDataFrame):
        landuse = landuse.copy()
    elif isinstance(landuse, pd.DataFrame):
        landuse = gpd.GeoDataFrame(landuse, geometry='geometry')
    else:
        print(
            'ERROR: landuse should either be a shapefile, a GeoDataFrame or a pandas Dataframe with a geometry column'
        )

    if isinstance(inun_file, str):
        # try to open it first as raster, if that doesnt work, we assume its a shapefile or GeoPackage that can be opened by geopandas
        try:
            with rasterio.open(inun_file) as src:
                if src.crs.to_dict() != landuse.crs:
                    landuse = landuse.to_crs(src.crs.to_dict())

                geoms = [mapping(geom) for geom in landuse.geometry]

                out_image, out_transform = mask(src, geoms, crop=True)
                out_image = np.array(out_image, dtype=int)

                # if inundation map is not in centimeters (and assumed to be in meters), multiple by 100
                if not centimeters:
                    out_image = out_image * 100
                out_image[out_image > 1000] = -1
                out_image[out_image <= 0] = -1
                gdf = None
        except:
            gdf = gpd.read_file(inun_file)
            out_image = None

    elif isinstance(inun_file, gpd.GeoDataFrame):
        gdf = inun_file.copy()
    elif isinstance(inun_file, pd.DataFrame):
        gdf = gpd.GeoDataFrame(inun_file, geometry='geometry')
    else:
        raise ValueError(
            'ERROR: inundation file should be a GeoTiff,  a shapefile, a GeoDataFrame \
              or any other georeferenced format that can be read by rasterio or geopandas'
        )

    # rename inundation colum, set values to centimeters if required and to integers
    if isinstance(gdf, gpd.GeoDataFrame):
        gdf = gdf.rename(columns={inun_col:'inun_val'})
        if not centimeters:
            gdf['inun_val'] = gdf.inun_val*100
            gdf['inun_val'] = gdf.inun_val.astype(int)

    # Load curves
    if isinstance(curve_path, pd.DataFrame):
        curves = curve_path.copy()
    elif isinstance(curve_path, np.ndarray):
        raise ValueError(
            'ERROR: for the vector-based approach we use a pandas DataFrame, not a Numpy Array'
        )
    elif curve_path.endswith('.csv'):
        curves = pd.read_csv(curve_path, index_col=[0])

    # Load maximum damages
    if isinstance(maxdam_path, str) and maxdam_path.endswith('.csv'):
        maxdam_path = pd.read_csv(maxdam_path)

    if isinstance(maxdam_path, pd.DataFrame):
        maxdam = dict(zip(maxdam_path[landuse_col], maxdam_path['damage']))
    elif isinstance(maxdam_path, np.ndarray):
        maxdam = dict(zip(maxdam_path[:, 0], maxdam_path[:, 1]))
    elif isinstance(maxdam_path, dict):
        maxdam = maxdam_path

    # convert raster to polygon
    if isinstance(out_image,np.ndarray):
        results = ({
            'properties': {
                'inun_val': v
            },
            'geometry': s
        } for i, (s, v) in enumerate(
            shapes(out_image[0, :, :], mask=None, transform=out_transform)))

        gdf = gpd.GeoDataFrame.from_features(list(results), crs=src.crs)

    # cut down to feasible area
    gdf = gdf.loc[gdf.inun_val > 0]
    gdf = gdf.loc[gdf.inun_val < 1000]

    # Check if we need to turn off tqdm:
    tqdm_print = kwargs.get('print_tqdm', True)

    # Split GeoDataFrame to make sure we have a unique shape per land use and inundation depth
    unique_df = []
    for row in tqdm(gdf.itertuples(index=False),
                    total=len(gdf),
                    desc='Get unique shapes',
                    disable=not tqdm_print):
        hits = landuse.loc[list(
            landuse.sindex.intersection(row.geometry.bounds))]
        row_buff = row.geometry.buffer(0)
        for hit in hits.itertuples(index=False):
            hit_buff = hit.geometry.buffer(0)
            if hit_buff.intersects(row_buff):
                unique_df.append([
                    row.inun_val, hit.landuse,
                    row_buff.intersection(hit_buff)
                ])

    # Create new dataframe
    tmp_df = pd.DataFrame(unique_df,
                          columns=['depth', landuse_col, 'geometry'])
    new_gdf = gpd.GeoDataFrame(tmp_df)

    # And remove empty geometries where there was no intersection in the end
    new_gdf = new_gdf.loc[new_gdf.geometry.geom_type != 'GeometryCollection']

    # Get area of shape
    new_gdf['area_m2'] = new_gdf.area

    # And estimate the losses
    if tqdm_print:
        tqdm.pandas(desc='Estimate damages')
        func = new_gdf.progress_apply
    else:
        func = new_gdf.apply

    new_gdf['damaged'] = func(lambda x: get_losses(x, curves, maxdam), axis=1)

    # Write the damages back to the original land-use shapes
    d_sindex = new_gdf.sindex

    landuse.reset_index()
    landuse['area_m2'] = landuse.area

    loss_dict = {}
    for x in tqdm(landuse.itertuples(),
                  total=len(landuse),
                  desc='Damage per object',
                  disable=not tqdm_print):
        hits = new_gdf.iloc[list(d_sindex.intersection(x.geometry.bounds))]
        damage = 0
        area_flooded = 0
        inun_levels = []
        x_buff = x.geometry.buffer(0)
        for hit in hits.itertuples(index=False):
            hit_buff = hit.geometry.buffer(0)
            if (x_buff.intersection(hit_buff)).area / hit_buff.area > 0.95:
                damage += hit.damaged
                area_flooded += hit.area_m2
                inun_levels.append(hit.depth)

        if len(inun_levels) == 0:
            loss_dict[x.Index] = 0, 0, 0, 0, 0
        else:
            loss_dict[x.Index] = damage, area_flooded, min(inun_levels), max(
                inun_levels), np.mean(inun_levels)

    tgt_cols = ['tot_dam', 'area_flooded', 'min_inun', 'max_inun', 'mean_inun']
    loss_df = pd.DataFrame.from_dict(loss_dict,
                                     orient='index',
                                     columns=tgt_cols)
    loss_gdf = gpd.GeoDataFrame(
        landuse.merge(loss_df, left_index=True, right_index=True))

    # If save is set to True, save original land-use map with damage values per shape.
    if save:
        # requires adding output_path and scenario_name to function call
        output_path = check_output_path(kwargs)
        scenario_name = check_scenario_name(kwargs)
        loss_gdf.to_file(
            p_join(output_path, 'damages_{}.shp'.format(scenario_name)))

    # And return the GeoDataFrame with damage statistics per unique object or shape
    return loss_gdf
