"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2023 Elco Koks. All versions released under the MIT license.
"""

# Get all the needed modules
import shapely
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from pathlib import Path
import osm_flex.extract as ex

from vector import VectorScanner
from raster import RasterScanner

class DamageScanner(object):
    """DamageScanner - a directe damage assessment toolkit	
    """
    def prepare(self, assessment_type, data_path, curves = pd.DataFrame(), maxdam = dict(),osm_pbf=False,**kwargs):
        """Prepare the input for a damage assessment
        """

        self.assessment_type = assessment_type
        
        self.data_path = data_path       

        # Collect exposure data
        try:
            if self.assessment_type == 'raster':
                self.exposure_data =  list(p.resolve() for p in Path(data_path / 'exposure').glob("**/*") if p.suffix in {".tif", ".tiff", ".nc"})[0]
            elif self.assessment_type == 'vector':
                get_data = list(p.resolve() for p in Path(data_path / 'exposure').glob("**/*") if p.suffix in {".gpkg", ".shp", ".geoparquet",".feather",".pbf"})
                if osm_pbf:
                    self.exposure_data =  list(p.resolve() for p in Path(data_path / 'exposure').glob("**/*") if p.suffix in {".pbf"})[0]
                    self.osm_pbf = osm_pbf
                elif len(get_data) > 0:
                    self.exposure_data =  [p for p in get_data if p.suffix in {".gpkg", ".shp", ".geoparquet",".feather"}][0]
                    self.osm_pbf = False
                elif len([p for p in get_data if p.suffix in {".gpkg", ".shp", ".geoparquet",".feather"}]) == 0:
                    self.exposure_data =  list(p.resolve() for p in Path(data_path / 'exposure').glob("**/*") if p.suffix in {".pbf"})[0]     
                    self.osm_pbf = osm_pbf
        except:
            raise ImportError(": No exposure data found. Either use the download() function to download OpenStreetMap data, or manually add data to the exposure folder")
        
        # Collect hazard data
        self.hazard_data = list(p.resolve() for p in Path(data_path / 'hazard').glob("**/*") if p.suffix in {".tif", ".tiff", ".nc"})

        if len(self.hazard_data) == 0:
            raise ImportError("NOTE: No hazard data found. Add data to the hazard folder first")

        # Collect vulnerability curves
        if len(curves) > 0:
            self.curves = curves
        else:
            self.curves = list(p.resolve() for p in Path(data_path / 'vulnerability').glob("**/*") if 'curve' in str(p).lower())[0]

        # Collect maxdam information
        if len(maxdam) > 0:
            self.maxdam = maxdam
        else:
            self.maxdam = list(p.resolve() for p in Path(data_path / 'vulnerability').glob("**/*") if (('dam' in str(p).lower()) |
                                                                                                       ('cost' in str(p).lower()))) [0]
       
        
    def single_event(self,hazard_type='',hazard_file='',return_period='',save_output=False,vector_landuse=True,**kwargs):
        """Damage assessment for a single event. Can be a specific hazard event, or a specific single hazard footprint 
        """

        if not hasattr(self, 'assessment_type'):
            raise ImportError('ERROR: Please run .prepare() first to set up the assessment')

        if hazard_file != '':
            self.hazard_map = [x for x in self.hazard_data if hazard_file in str(x)][0]
        else:      
            if hazard_type == '':
                if len(self.hazard_data) == 1:
                    self.hazard_map = self.hazard_data[0]
                else:
                    raise RuntimeError('Multiple events require the .multiple_events() function!')     
            elif hazard_type != '':
                if hazard_type == 'flood':
                    hazard_characteristics = ['inun','flood','fluvial','pluvial','coastal','fd','pd']
                elif hazard_type == 'windstorm':
                    hazard_characteristics = ['storm','wind','cyclone','tc']
                elif hazard_type == 'earthquake':
                    hazard_characteristics = ['eq','earthquake']
                elif hazard_type == 'landslide':
                    hazard_characteristics = ['ls','landslide']
                
                hazard_maps = [x for x in self.hazard_data if any(y in str(x).lower() for y in hazard_characteristics)]
                
                if return_period != '':
                    self.hazard_map = [x for x in hazard_maps if return_period in str(x)].lower()[0]
                
                elif return_period == '':
                    self.hazard_map = hazard_maps[0]
                
                else:
                    raise RuntimeError('Specify filename or flood type and return period to be able to run a single event analysis. Multiple events require the .multiple_events() function.')           
            
        if self.assessment_type == 'raster':

            return RasterScanner(
                            exposure_file = self.exposure_data,
                            hazard_file = self.hazard_map,
                            curve_path = self.curves,
                            maxdam_path = self.maxdam,
                            save=save_output)

        elif self.assessment_type == 'vector':

            # specificy essential data input characteristics
            if 'cell_size' in kwargs:
                self.cell_size = kwargs.get('cell_size')#0.01666,
            else:
                self.cell_size = 5
    
            if 'exp_crs' in kwargs:
                self.exp_crs = kwargs.get('exp_crs')
            else:
                self.exp_crs = 4326
    
            if 'haz_crs' in kwargs:
                self.haz_crs = kwargs.get('haz_crs')
            else:
                self.haz_crs = 4326
    
            if 'object_col' in kwargs:
                self.object_col = kwargs.get('object_col')
            else:
                self.object_col = 'landuse'
    
            if 'hazard_col' in kwargs:
                self.hazard_col = kwargs.get('hazard_col')
            else:
                self.hazard_col = 'inun_val'
    
            if 'lat_col' in kwargs:
                self.lat_col = kwargs.get('lat_col')
            else:
                self.lat_col = 'y'
    
            if 'lon_col' in kwargs:
                self.lon_col = kwargs.get('lon_col')
            else:
                self.lon_col = 'x'
    
            if 'centimers' in kwargs:
                self.centimers = kwargs.get('centimers')#0.01666,#5,
            else:
                self.centimers = False
            
            if self.osm_pbf:
                if 'buildings' in kwargs:
                    self.exposure_data = DamageScanner.buildings(self)                  
                    self.exposure_data = self.exposure_data.rename({'building':'element_type'},axis=1)       
                    
                elif 'roads' in kwargs:
                    self.exposure_data = DamageScanner.cis(self,infra_type='road')
                    self.exposure_data = self.exposure_data.rename({'highway':'element_type'},axis=1)       
                
                elif vector_landuse:
                    self.exposure_data = DamageScanner.landuse(self)
                    self.exposure_data = self.exposure_data.rename({'landuse':'element_type'},axis=1)       

            
            return VectorScanner(
                            exposure_file = self.exposure_data,
                            hazard_file = self.hazard_map,
                            curve_path = self.curves,
                            maxdam_path = self.maxdam,
                            cell_size = self.cell_size, #0.01666,#5,
                            exp_crs = self.exp_crs, #28992,
                            haz_crs = self.haz_crs, #4326,
                            object_col= self.object_col, #'landuse',
                            hazard_col= self.hazard_col,
                            lat_col = self.lat_col,
                            lon_col = self.lon_col,
                            centimeters= self.centimers, #False,
                            save=save_output)  

                
    def multiple_events(self):
        """Damage assessment for multiple events. Also the input for the risk assessment
        """
        pass

    def risk(self):
        """Risk assessment
        """
        pass

    def landuse(self):
        """
        """
        extracted_exposure = ex.extract(
                            	self.exposure_data, 'multipolygons', ['landuse'])    
        
        extracted_exposure.geometry = shapely.make_valid(extracted_exposure.geometry)        

        return extracted_exposure
    
    def buildings(self):
        """
        """
        extracted_exposure = ex.extract(
                            	self.exposure_data, 'multipolygons', ['building'])    
        
        extracted_exposure.geometry = shapely.make_valid(extracted_exposure.geometry)        

        return extracted_exposure

    def cis(self,infra_type):
        """
        """
        extracted_exposure = ex.extract_cis(self.exposure_data, infra_type)

        if infra_type == 'road':
            extracted_exposure = extracted_exposure.loc[extracted_exposure.geometry.geom_type == 'LineString']
        
        extracted_exposure.geometry = shapely.make_valid(extracted_exposure.geometry)        

        return extracted_exposure
        
    def cis_all(self,to_exclude=[]):
        """
        """
        cis = ['healthcare', 'education', 'gas', 'oil', 'telecom', 'water', 'wastewater', 'power', 'rail', 'road', 'air']
        pass

    def download(self, country_code='JAM'):
        """
        """        
        pass