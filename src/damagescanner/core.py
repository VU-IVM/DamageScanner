"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2023 Elco Koks. All versions released under the MIT license.
"""

# Get all the needed modules
import shapely
import xarray as xr
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from pathlib import Path
from scipy import integrate


from vector import VectorScanner, VectorExposure
from raster import RasterScanner

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


class DamageScanner(object):
    """DamageScanner - a directe damage assessment toolkit"""

    def __init__(self, hazard_data, exposure_data, curves, maxdam):
        """Prepare the input for a damage assessment"""

        # Collect the input data
        self.hazard_data = hazard_data
        self.exposure_data = exposure_data
        self.curves = curves
        self.maxdam = maxdam

        # Convert the input to a Path object if it is a string
        if isinstance(self.exposure_data, str):
            self.exposure_path = Path(exposure_data)

        if isinstance(self.hazard_data, str):
            self.hazard_data = Path(hazard_data)

        if isinstance(self.curves, str):
            self.curves = Path(curves)

        if isinstance(self.maxdam, str):
            self.maxdam = Path(maxdam)

        # Check the type of the exposure data
        if isinstance(self.exposure_data, Path):
            if self.exposure_data.suffix in [".tif", ".tiff", ".nc"]:
                self.assessment_type = "raster"

            elif self.exposure_data.suffix in [
                ".shp",
                ".gpkg",
                ".pbf",
                ".geofeather",
                ".geoparquet",
            ]:
                self.assessment_type = "vector"

                if self.exposure_data.suffix == ".pbf":
                    self.osm = True
            else:
                raise ImportError(
                    ": The exposure data should be a a shapefile, geopackage, \
                        geoparquet, osm.pbf, geotiff or netcdf file."
                )

        else:
            if isinstance(self.exposure_data, (xr.DataArray, xr.Dataset)):
                self.assessment_type = "raster"
            elif isinstance(self.exposure_data, (gpd.GeoDataFrame, pd.DataFrame)):
                self.assessment_type = "vector"

        # Collect vulnerability curves
        if isinstance(curves, (pd.DataFrame, Path)):
            self.curves = curves
        else:
            raise ImportWarning(
                "Prepare the vulnerability curves as a pandas DataFrame or a \
                     as a directory path to a csv file"
            )

        # Collect maxdam information
        if isinstance(maxdam, (pd.DataFrame, Path)):
            self.maxdam = maxdam
        else:
            raise ImportWarning(
                "Prepare the maximum damages as a pandas DataFrame or or a\
                     as a directory path to a csv file"
            )

    def exposure(self, **kwargs):
        """Exposure data"""

        if self.assessment_type == "raster":
            return xr.open_rasterio(self.exposure_data)

        elif self.assessment_type == "vector":
            # specificy essential data input characteristics
            if "asset_type" in kwargs:
                self.asset_type = kwargs.get("asset_type")
            else:
                self.asset_type = "landuse"

            exposed_assets = VectorExposure(
                hazard_file=self.hazard_data,
                exposure_file=self.exposure_data,
                asset_type=self.asset_type,
            )[0]

            return exposed_assets

    def calculate(self, save_output=False, **kwargs):
        """Damage assessment. Can be a specific hazard event, or a specific \
            single hazard footprint, or a list of events/footprints."""

        if not hasattr(self, "assessment_type"):
            raise ImportError("Please prepare the input data first")

        if self.assessment_type == "raster":
            return RasterScanner(
                exposure_file=self.exposure_data,
                hazard_file=self.hazard_data,
                curve_path=self.curves,
                maxdam_path=self.maxdam,
                save=save_output,
            )

        elif self.assessment_type == "vector":
            # specificy essential data input characteristics
            if "asset_type" in kwargs:
                self.asset_type = kwargs.get("asset_type")
            else:
                self.asset_type = "landuse"

            return VectorScanner(
                hazard_file=self.hazard_data,
                exposure_file=self.exposure_data,
                curve_path=self.curves,
                maxdam_path=self.maxdam,
                asset_type=self.asset_type,  #'landuse',
                sub_types=kwargs.get("subtypes", None),
                save=save_output,
            )

    def risk(hazard_dict, exposure, curves, maxdam):
        """
        Calculate the risk for a list of hazard events
        """

        RP_list = list(hazard_dict.keys())

        risk = {}
        for key, value in tqdm(
            hazard_dict.items(), total=len(hazard_dict), desc="Risk Calculation"
        ):
            risk[key] = DamageScanner(value, exposure, curves, maxdam).calculate()[0]

        df_risk = pd.concat(risk, axis=1)

        RPS = [1 / x for x in RP_list]

        risk = pd.DataFrame(
            df_risk.apply(
                lambda x: integrate.simpson(y=x[RP_list][::-1], x=RPS[::-1]), axis=1
            ),
            columns=["tot_risk"],
        )

        return risk


if __name__ == "__main__":
    ####################################################################################################

    # Kampen

    data_path = Path("..") / ".." / "data" / "kampen"

    # define the input data
    exposure = data_path / "exposure" / "landuse_map.tif"
    hazard = data_path / "hazard" / "1in100_inundation_map.tif"
    curves = data_path / "vulnerability" / "curves.csv"
    maxdam = data_path / "vulnerability" / "maxdam.csv"

    # initiate the damage scanner and calculate the damages
    # print(DamageScanner(hazard, exposure, curves, maxdam).calculate()[0])

    # define the input data
    exposure = data_path / "exposure" / "landuse.shp"
    hazard = data_path / "hazard" / "1in100_inundation_map.tif"
    curves = data_path / "vulnerability" / "curves_osm.csv"
    maxdam = data_path / "vulnerability" / "maxdam_osm.csv"

    # initiate the damage scanner and calculate the damages
    # print(DamageScanner(hazard, exposure, curves, maxdam).calculate())

    ####################################################################################################
    # Kampen risk assessment
    data_path = Path("..") / ".." / "data" / "kampen"

    # define the input data
    hazard_dict = {
        10: data_path / "hazard" / "1in10_inundation_map.tif",
        50: data_path / "hazard" / "1in50_inundation_map.tif",
        100: data_path / "hazard" / "1in100_inundation_map.tif",
        500: data_path / "hazard" / "1in500_inundation_map.tif",
        1000: data_path / "hazard" / "1in1000_inundation_map.tif",
    }

    exposure = data_path / "exposure" / "landuse_map.tif"
    curves = data_path / "vulnerability" / "curves.csv"
    maxdam = data_path / "vulnerability" / "maxdam.csv"

    # calculate the risk
    # print(DamageScanner.risk(hazard_dict, exposure, curves, maxdam))

    ####################################################################################################
    # Jamaica

    data_path = Path("..") / ".." / "data" / "jamaica"

    # define the input data
    exposure = data_path / "exposure" / "jamaica-latest.osm.pbf"
    hazard = data_path / "hazard" / "FD_1in1000.tif"
    curves = data_path / "vulnerability" / "curves_osm.csv"
    maxdam = data_path / "vulnerability" / "maxdam_osm.csv"

    # estimate exposure

    exposed_features = DamageScanner(hazard, exposure, curves, maxdam).exposure(
        asset_type="power"
    )

    #exposed_features.to_parquet("main_roads.parquet")
    print(exposed_features[["geometry","coverage","values"]])

    #initiate the damage scanner and calculate the damages
    # print(
    #     DamageScanner(hazard, exposure, curves, maxdam).calculate(asset_type="main_roads").damage.sum()

    # )
