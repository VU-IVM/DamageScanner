"""Test core objects/concepts
"""
# pylint: disable=C0103
from geopandas import GeoDataFrame
from pandas.testing import assert_frame_equal
from pytest import fixture
from shapely.geometry import Point, LineString, MultiPoint

#import damagescanner
#from damagescanner import RasterScanner

def inc(x):
    return x + 1

def test_answer():
    assert inc(3) == 4