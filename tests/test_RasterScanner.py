"""Test core objects/concepts
"""
import unittest
import numpy as np
import pandas as pd
from damagescanner.core import RasterScanner  
from pathlib import Path

class TestRasterScanner(unittest.TestCase):

    def test_raster_scanner(self):
        # Define paths to example files
        data_path = Path(__file__).parent.parent / 'data'
        landuse_file = str(data_path / 'landuse' / 'landuse_map.tif')
        hazard_file = str(data_path / 'hazard' / 'inundation_map.tif')
        curve_path = str(data_path / 'curves' / 'curves.csv')
        maxdam_path = str(data_path / 'curves' / 'maxdam.csv')

        # Call the RasterScanner function
        damage_df, damagemap, landuse_in, hazard = RasterScanner(
            landuse_file=landuse_file,
            hazard_file=hazard_file,
            curve_path=curve_path,
            maxdam_path=maxdam_path,
            lu_crs=28992,
            haz_crs=4326,
            dtype=np.int32,
            save=False
            # Add any other required parameters or kwargs here
        )

        # Perform assertions to check the correctness of the results
        self.assertIsInstance(damage_df, pd.DataFrame)
        self.assertIsInstance(damagemap, np.ndarray)
        self.assertIsInstance(landuse_in, np.ndarray)
        self.assertIsInstance(hazard, np.ndarray)

        # Add more specific assertions based on your requirements

if __name__ == "__main__":
    TESTS = unittest.TestLoader().loadTestsFromTestCase(TestRasterScanner)
    unittest.TextTestRunner(verbosity=2).run(TESTS)