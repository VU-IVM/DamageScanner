"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2023 Elco Koks. All versions released under the MIT license.
"""

# Get all the needed modules
import os
import rasterio
import xarray as xr
import numpy as np
import shapely
import pandas as pd
import geopandas as gpd
from affine import Affine
import pyproj
from tqdm import tqdm
import warnings
from pathlib import Path,PurePath

from vector import VectorScanner
from raster import RasterScanner

class DamageScanner(object):
    """DamageScanner - a directe damage assessment toolkit	
    """
    def prepare(self, assessment_type, data_path):
        """Prepare the input for a damage assessment
        """

        self.assessment_type = assessment_type
        
        self.data_path = data_path       

        # Collect exposure data
        try:
            if self.assessment_type == 'raster':
                self.exposure_data =  list(p.resolve() for p in Path(data_path / 'exposure').glob("**/*") if p.suffix in {".tif", ".tiff", ".nc"})[0]
            elif self.assessment_type == 'vector':
                self.exposure_data =  list(p.resolve() for p in Path(data_path / 'exposure').glob("**/*") if p.suffix in {".gpkg", ".shp", ".geoparquet",".feather"})[0]
        except:
            raise ImportError(": No exposure data found. Either use the download() function to download OpenStreetMap data, or manually add data to the exposure folder")
        
        # Collect hazard data
        self.hazard_data = list(p.resolve() for p in Path(data_path / 'hazard').glob("**/*") if p.suffix in {".tif", ".tiff", ".nc"})

        if len(self.hazard_data) == 0:
            raise ImportError("NOTE: No hazard data found. Add data to the hazard folder first")

        # Collect vulnerability curves
        self.curves = list(p.resolve() for p in Path(data_path / 'vulnerability').glob("**/*") if 'curve' in str(p))[0]

        # Collect maxdam information
        self.maxdam = list(p.resolve() for p in Path(data_path / 'vulnerability').glob("**/*") if 'dam' in str(p))[0]
        
    def single_event(self,hazard_type='',file_name='',save_output=False):
        """Damage assessment for a single event. Can be a specific hazard event, or a specific single hazard footprint 
        """

        if not hasattr(self, 'assessment_type'):
            raise ImportError('ERROR: Please run .prepare() first to set up the assessment')

        if hazard_type == '':
            if len(self.hazard_data) == 1:
                self.hazard_map = self.hazard_data[0]
            else:
                raise RuntimeError('Multiple events require the .multiple_events() function!')     
                
        
        if self.assessment_type == 'raster':

            return RasterScanner(exposure_file = self.exposure_data,
                          hazard_file = self.hazard_map,
                          curve_path = self.curves,
                          maxdam_path = self.maxdam,
                          save=save_output)

    def multiple_events(self):
        """Damage assessment for multiple events. Also the input for the risk assessment
        """
        pass

    def risk(self):
        """Risk assessment
        """
        pass

    def buildings(self):
        """
        """
        pass

    def cis(self,infra_type):
        """
        """
        pass

    def cis_all(self,to_exclude=[]):
        pass