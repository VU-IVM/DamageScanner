"""Test core objects/concepts
"""
# pylint: disable=C0103
from geopandas import GeoDataFrame
from pandas.testing import assert_frame_equal
from pytest import fixture
from shapely.geometry import Point, LineString, MultiPoint

import os

#import damagescanner
from damagescanner import RasterScanner

data_path = '..'

#def test_answer():
#    inun_map = os.path.join(data_path,'data','inundation','inundation_map.tif')
#    landuse_map = os.path.join(data_path,'data','landuse','landuse_map.tif')
#    
#    curve_path = os.path.join(data_path,'data','curves','curves.csv')
#    maxdam_path = os.path.join(data_path,'data','curves','maxdam.csv')
#    
#    assert RasterScanner(landuse_map,inun_map,curve_path,maxdam_path)[0].sum().values[0] == 3550789426