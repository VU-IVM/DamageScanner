"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2023 Elco Koks. All versions released under the MIT license.
"""

# Get all the needed modules
import shapely
import pandas as pd
from pathlib import Path

from .vector import VectorScanner, buildings, landuse, cis
from .raster import RasterScanner


class DamageScanner(object):
    """DamageScanner - a directe damage assessment toolkit"""

    def __init__(self, data_path, exposure_data, hazard_data, curves=None, maxdam=None):
        """Prepare the input for a damage assessment"""

        # specify the basics
        self.data_path = data_path
        self.osm = False
        self.default_curves = False
        self.default_maxdam = False

        # Specify paths to exposure data
        self.exposure_path = Path(data_path / "exposure" / exposure_data)

        if self.exposure_path.suffix in [".tif", ".tiff", ".nc"]:
            self.assessment_type = "raster"

        elif self.exposure_path.suffix in [
            ".shp",
            ".gpkg",
            ".pbf",
            ".geofeather",
            ".geoparquet",
        ]:
            self.assessment_type = "vector"

            if self.exposure_path.suffix == ".pbf":
                self.osm = True
        else:
            raise ImportError(
                ": The exposure data should be a a shapefile, geopackage, geoparquet, osm.pbf, geotiff or netcdf file."
            )

        # Specify path to hazard data
        if not isinstance(hazard_data, list):
            self.hazard_path = Path(data_path / "hazard" / hazard_data)
        else:
            self.hazard_path = [Path(data_path / "hazard" / x) for x in hazard_data]

        # Collect vulnerability curves
        if curves is None:

            self.curves = "https://zenodo.org/records/10203846/files/Table_D2_Multi-Hazard_Fragility_and_Vulnerability_Curves_V1.0.0.xlsx?download=1"

            raise ImportWarning(
                "You have decided to choose the default set of curves. Make sure the exposure data aligns with these curves: https://zenodo.org/records/10203846"
            )
            self.default_curves = True

        elif isinstance(curves, pd.DataFrame):
            self.curves = curves
        else:
            self.curves = Path(data_path / "vulnerability" / curves)

        # Collect maxdam information
        if maxdam is None:
            self.maxdam = "https://zenodo.org/records/10203846/files/Table_D3_Costs_V1.0.0.xlsx?download=1"

            raise ImportWarning(
                "You have decided to choose the default set of maximum damages. Make sure the exposure data aligns with these maximum damages: https://zenodo.org/records/10203846"
            )
            self.default_maxdam = True

        elif isinstance(maxdam, dict):
            self.maxdam = maxdam
        else:
            self.maxdam = Path(data_path / "vulnerability" / maxdam)

    # def exposure(self):

    # def vulnerability(self):

    def calculate(self, hazard_type=None, save_output=False, **kwargs):
        """Damage assessment. Can be a specific hazard event, or a specific single hazard footprint, or a list of events/footprints."""

        if not hasattr(self, "assessment_type"):
            raise ImportError("Please run .prepare() first to set up the assessment")

        if self.assessment_type == "raster":

            return RasterScanner(
                exposure_file=self.exposure_path,
                hazard_file=self.hazard_path,
                curve_path=self.curves,
                maxdam_path=self.maxdam,
                save=save_output,
            )

        elif self.assessment_type == "vector":

            if self.default_curves:
                if hazard_type == "flood":
                    sheet_name = "F_Vuln_Depth"
                elif hazard_type == "windstorm":
                    sheet_name = "W_Vuln_V10m"

            # specificy essential data input characteristics
            if "cell_size" in kwargs:
                self.cell_size = kwargs.get("cell_size")  # 0.01666,
            else:
                self.cell_size = 5

            if "exp_crs" in kwargs:
                self.exp_crs = kwargs.get("exp_crs")
            else:
                self.exp_crs = 4326

            if "haz_crs" in kwargs:
                self.haz_crs = kwargs.get("haz_crs")
            else:
                self.haz_crs = 4326

            if "object_col" in kwargs:
                self.object_col = kwargs.get("object_col")
            else:
                self.object_col = "landuse"

            if "hazard_col" in kwargs:
                self.hazard_col = kwargs.get("hazard_col")
            else:
                self.hazard_col = "inun_val"

            if "lat_col" in kwargs:
                self.lat_col = kwargs.get("lat_col")
            else:
                self.lat_col = "y"

            if "lon_col" in kwargs:
                self.lon_col = kwargs.get("lon_col")
            else:
                self.lon_col = "x"

            if "centimers" in kwargs:
                self.centimers = kwargs.get("centimers")  # 0.01666,#5,
            else:
                self.centimers = False

            if self.osm:
                if "buildings" in kwargs:
                    self.exposure_data = buildings(self)
                    self.exposure_data = self.exposure_data.rename(
                        {"building": "element_type"}, axis=1
                    )

                elif "roads" in kwargs:
                    self.exposure_data = cis(self, infra_type="road")
                    self.exposure_data = self.exposure_data.rename(
                        {"highway": "element_type"}, axis=1
                    )

                elif "landuse" in kwargs:
                    self.exposure_data = landuse(self)
                    self.exposure_data = self.exposure_data.rename(
                        {"landuse": "element_type"}, axis=1
                    )
                else:
                    raise RuntimeError(
                        "When using OSM data, you need to specify the object type (e.g. road, buildings, landuse)"
                    )
            else:
                self.exposure_data = self.exposure_path

            return VectorScanner(
                exposure_file=self.exposure_data,
                hazard_file=self.hazard_path,
                curve_path=self.curves,
                maxdam_path=self.maxdam,
                cell_size=self.cell_size,  # 0.01666,#5,
                exp_crs=self.exp_crs,  # 28992,
                haz_crs=self.haz_crs,  # 4326,
                object_col=self.object_col,  #'landuse',
                hazard_col=self.hazard_col,
                lat_col=self.lat_col,
                lon_col=self.lon_col,
                centimeters=self.centimers,  # False,
                save=save_output,
            )

    def risk(self):
        """Risk assessment"""
        pass


if __name__ == "__main__":

    data_path = Path("..")

    kampen = DamageScanner(
        data_path=data_path / "data" / "kampen",
        exposure_data="landuse_map.tif",
        hazard_data="inundation_map.tif",
        curves="curves.csv",
        maxdam="maxdam.csv",
    )

    output = kampen.calculate()
