"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2023 Elco Koks. All versions released under the MIT license.
"""

# Get all the needed modules
import rasterio
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from pathlib import Path
from scipy import integrate

from damagescanner.vector import VectorScanner, VectorExposure
from damagescanner.raster import RasterScanner

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)


class DamageScanner(object):
    """
    DamageScanner - a direct damage assessment toolkit.

    This class provides tools to assess direct physical damage from hazard events,
    using raster or vector exposure data and vulnerability curves.
    It supports both single-hazard footprints and risk-based multi-scenario assessments.
    """

    def __init__(self, hazard_data, feature_data, curves, maxdam):
        """
        Initialize the DamageScanner class with hazard, exposure, curve, and max damage data.

        Args:
            hazard_data (str | Path | xarray.DataArray | xarray.Dataset): Path to raster hazard file or xarray object.
            feature_data (str | Path | pd.DataFrame | gpd.GeoDataFrame): Exposure data, either raster or vector.
            curves (str | Path | pd.DataFrame): Vulnerability curves as DataFrame or CSV file path.
            maxdam (str | Path | pd.DataFrame): Maximum damage values per asset type.
        """
        # Collect the input data
        self.hazard_data = hazard_data
        self.feature_data = feature_data
        self.curves = curves
        self.maxdam = maxdam

        # Convert the input to a Path object if it is a string
        if isinstance(self.feature_data, str):
            self.feature_data = Path(feature_data)

        if isinstance(self.hazard_data, str):
            self.hazard_data = Path(hazard_data)

        if isinstance(self.curves, str):
            self.curves = Path(curves)

        if isinstance(self.maxdam, str):
            self.maxdam = Path(maxdam)

        # Check the type of the exposure data
        if isinstance(self.feature_data, Path):
            if self.feature_data.suffix in [".tif", ".tiff", ".nc"]:
                self.assessment_type = "raster"

            elif self.feature_data.suffix in [
                ".shp",
                ".gpkg",
                ".pbf",
                ".geofeather",
                ".geoparquet",
            ]:
                self.assessment_type = "vector"

                if self.feature_data.suffix == ".pbf":
                    self.osm = True
            else:
                raise ImportError(
                    ": The exposure data should be a a shapefile, geopackage, \
                        geoparquet, osm.pbf, geotiff or netcdf file."
                )

        else:
            if isinstance(self.feature_data, (xr.DataArray, xr.Dataset)):
                self.assessment_type = "raster"
            elif isinstance(self.feature_data, (gpd.GeoDataFrame, pd.DataFrame)):
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

    def exposure(self, disable_progress=False, output_path=None, **kwargs):
        """
            Run the exposure analysis to identify features affected by the hazard footprint.

            This method analyzes the input data to determine which features are exposed to a hazard.
            It supports both raster and vector data types for the exposure analysis. If the data type
            is vector, additional keyword arguments can be specified to customize the analysis.

            Args:
                disable_progress (bool, optional): If True, disables progress bars during processing.
                    Defaults to False.
                output_path (str, optional): The file path to save the exposure data. The file format
                    is determined by the file extension (e.g., '.parquet', '.csv', '.gpkg', '.shp').
                    If None, the data is not saved. Defaults to None.
                **kwargs: Optional keyword arguments:
                    - asset_type (str): The type of asset to evaluate (only for vector data).

            Returns:
                geopandas.GeoDataFrame | xarray.DataArray: A GeoDataFrame containing the affected
                assets if the input data is vector, or an xarray.DataArray if the input data is raster.

            Notes:
                - If `output_path` is provided, the method saves the exposure data to the specified path.
                - The file format is inferred from the file extension of `output_path`. If the extension
                is not recognized, the data is saved as a Parquet file by default.
        """
        if self.assessment_type == "raster":
            return xr.open_rasterio(self.exposure_data)

        elif self.assessment_type == "vector":
            # specificy essential data input characteristics
            if "asset_type" in kwargs:
                self.asset_type = kwargs.get("asset_type")
            else:  ## DO WE WANT THIS?!?! or should this always be defined?!
                self.asset_type = "landuse"

            exposed_assets = VectorExposure(
                hazard_file=self.hazard_data,
                feature_file=self.feature_data,
                asset_type=self.asset_type,
                disable_progress=disable_progress,
            )[0]

            # save output when exposed assets are empty
            if output_path:
                # Determine the file format based on the file extension
                file_extension = output_path.split('.')[-1].lower()
                format_mapping = {
                    'parquet': exposed_assets.to_parquet,
                    'csv': exposed_assets.to_csv,
                    'gpkg': exposed_assets.to_file,
                    'shp': exposed_assets.to_file,
                }

                # Default to parquet if the extension is not recognized
                save_function = format_mapping.get(file_extension, exposed_assets.to_parquet)
                save_function(output_path)

                print(f'Exposure data saved to {output_path}')

            return exposed_assets

    def calculate(self, disable_progress=False, output_path=None, **kwargs):
        """
        Perform a damage calculation using the provided inputs.

        Applies vulnerability curves and maximum damage values to the exposed features
        or raster grid to calculate expected damage.

        Args:
            disable_progress (bool, optional): If True, disables progress bars. Defaults to False.
            output_path (str, optional): Path to save the calculation results. The file format
                is determined by the file extension (e.g., '.csv', '.parquet' for vector data,
                '.tif' for raster data). If None, the data is not saved. Defaults to None.
            **kwargs:
                asset_type (str, optional): Infrastructure class to evaluate.
                multi_curves (dict, optional): Mapping of asset types to curve sets.
                subtypes (list, optional): Used for subtype analysis.

        Returns:
            pd.DataFrame | xr.DataArray: Estimated damages for each asset or grid cell.
        """
        if not hasattr(self, "assessment_type"):
            raise ImportError("Please prepare the input data first")

        if self.assessment_type == "raster":
            damage_df, damagemap = RasterScanner(
                exposure_file=self.feature_data,
                hazard_file=self.hazard_data,
                curve_path=self.curves,
                maxdam_path=self.maxdam,
            )

            # Extract CRS and transform from feature_data
            if isinstance(self.feature_data, (str, Path)):
                # Assume it's a file path
                if self.feature_data.endswith('.nc'):
                    # Open with xarray if it's a NetCDF file
                    feature_data_xr = xr.open_dataset(self.feature_data)
                    crs = feature_data_xr.rio.crs  # Requires rioxarray extension
                    transform = feature_data_xr.rio.transform()
                else:
                    # Open with rasterio for other raster formats
                    with rasterio.open(self.feature_data) as src:
                        crs = src.crs
                        transform = src.transform
            elif isinstance(self.feature_data, (xr.DataArray, xr.Dataset)):
                # Directly use the xarray object
                crs = self.feature_data.rio.crs  # Requires rioxarray extension
                transform = self.feature_data.rio.transform()
            else:
                raise ValueError("Unsupported feature_data format")

        elif self.assessment_type == "vector":
            # Specify essential data input characteristics
            if "asset_type" in kwargs:
                self.asset_type = kwargs.get("asset_type")
            else:
                self.asset_type = None

            damage_df = VectorScanner(
                hazard_file=self.hazard_data,
                feature_file=self.feature_data,
                curve_path=self.curves,
                maxdam_path=self.maxdam,
                asset_type=self.asset_type,
                multi_curves=kwargs.get("multi_curves", None),
                sub_types=kwargs.get("subtypes", None),
                disable_progress=disable_progress,
            )

            # For vector data, CRS and transform are not directly applicable
            crs = None
            transform = None

        if output_path:
            file_extension = output_path.split('.')[-1].lower()
            if self.assessment_type == 'vector':
                format_mapping = {
                    'csv': damage_df.to_csv,
                    'parquet': damage_df.to_parquet,
                }
                save_function = format_mapping.get(file_extension, damage_df.to_csv)
                save_function(output_path, **kwargs)
            elif self.assessment_type == 'raster':
                # Save the damage_df as CSV
                damage_df_path = output_path.replace(f'.{file_extension}', '_damages.csv')
                damage_df.to_csv(damage_df_path)
                print(f'Damage summary saved to {damage_df_path}')

                # Save the damagemap as GeoTIFF
                dmap_fn = output_path
                rst_opts = {
                    "driver": "GTiff",
                    "height": damagemap.shape[0],
                    "width": damagemap.shape[1],
                    "count": 1,
                    "dtype": damagemap.dtype,
                    "crs": crs,
                    "transform": transform,
                    "compress": "LZW",
                }
                with rasterio.open(dmap_fn, "w", **rst_opts) as dst:
                    dst.write(damagemap, 1)
                print(f'Damage map saved to {dmap_fn}')

        return damage_df if self.assessment_type == "vector" else (damage_df, damagemap)
    

    def risk(self, hazard_dict, output_path=None, **kwargs):
        """
        Perform a risk assessment across multiple hazard return periods.

        Integrates damages from each return period and computes expected annual damages.
        Supports both single and multi-curve inputs for infrastructure types.

        Args:
            hazard_dict (dict): Dictionary mapping return periods to hazard raster paths.
            output_path (str, optional): Path to save the risk assessment results. The file format
                is determined by the file extension (e.g., '.csv', '.parquet'). If None, the data is not saved.
            **kwargs:
                asset_type (str, optional): Infrastructure class to evaluate.
                multi_curves (dict, optional): Mapping of asset types to curve sets.

        Returns:
            pd.DataFrame | None: A GeoDataFrame with risk values for each asset, or None if no results.
        """
        RP_list = list(hazard_dict.keys())

        risk = {}
        for key, hazard_map in tqdm(
            hazard_dict.items(), total=len(hazard_dict), desc="Risk Calculation"
        ):
            if self.assessment_type == "raster":
                risk[key] = DamageScanner(
                    hazard_map, self.feature_data, self.curves, self.maxdam
                ).calculate(disable_progress=True)[0]
            else:
                if kwargs.get("asset_type", None) is not None:
                    risk[key] = DamageScanner(
                        hazard_map, self.feature_data, self.curves, self.maxdam
                    ).calculate(
                        disable_progress=True, asset_type=kwargs.get("asset_type")
                    )
                elif kwargs.get("multi_curves", None) is not None:
                    risk[key] = DamageScanner(
                        hazard_map, self.feature_data, self.curves, self.maxdam
                    ).calculate(
                        disable_progress=True, multi_curves=kwargs.get("multi_curves")
                    )
                else:
                    risk[key] = DamageScanner(
                        hazard_map, self.feature_data, self.curves, self.maxdam
                    ).calculate(disable_progress=True)

        # Collect the risk for each RP
        df_risk = pd.concat(risk, axis=1)

        if (len(df_risk) == 0) or (df_risk.isnull().all().all()):
            return None

        # Get the dataframe of the largest RP
        largest_rp = df_risk.loc[:, pd.IndexSlice[RP_list[-1], :]]

        if kwargs.get("multi_curves", None) is None:
            # only keep the damage values
            df_risk = df_risk.loc[:, pd.IndexSlice[RP_list, "damage"]].fillna(0)

            RPS = [1 / x for x in RP_list]

            risk = pd.DataFrame(
                df_risk.apply(
                    lambda x: np.trapezoid(y=x[RP_list][::-1], x=RPS[::-1]), axis=1
                ),
                columns=["tot_risk"],
            )

            # save output when tot_risk returns negative values
            if risk.tot_risk.min() < 0:
                df_risk.to_csv("df_risk.csv")
                risk.to_csv("risk.csv")

            # Save the risk to the largest RP
            largest_rp.columns = largest_rp.columns.get_level_values(1)
            largest_rp = largest_rp.drop("damage", axis=1)
            largest_rp.loc[:, "risk"] = risk.values

            # Save the results if output_path is provided
            if output_path:
                file_extension = output_path.split('.')[-1].lower()
                format_mapping = {
                    'csv': largest_rp.to_csv,
                    'parquet': largest_rp.to_parquet,
                }
                save_function = format_mapping.get(file_extension, largest_rp.to_csv)
                save_function(output_path)
                print(f'Risk assessment results saved to {output_path}')

            # return the risk in a concise dataframe
            return largest_rp[["osm_id", "object_type", "geometry", "risk"]]

        else:
            multi_curves = kwargs.get("multi_curves")

            # only keep the damage values
            df_risk = df_risk.loc[
                :, pd.IndexSlice[RP_list, multi_curves.keys()]
            ].fillna(0)

            RPS = [1 / x for x in RP_list]

            # estimate risks
            collect_risks = {}

            for curve in multi_curves.keys():
                subrisk = df_risk.loc[:, pd.IndexSlice[:, curve]]
                collect_risks[curve] = subrisk.apply(
                    lambda x: np.trapz(y=x[RP_list][::-1], x=RPS[::-1]), axis=1
                ).values

                # save output when tot_risk returns negative values
                if any(subrisk.min() < 0):
                    df_risk.to_csv("df_risk.csv")
                    subrisk.to_csv("risk.csv")

            all_risks = pd.DataFrame.from_dict(collect_risks)

            largest_rp.columns = largest_rp.columns.get_level_values(1)
            largest_rp = largest_rp.drop(multi_curves.keys(), axis=1)
            largest_rp.loc[:, multi_curves.keys()] = all_risks.values

            # Save the results if output_path is provided
            if output_path:
                file_extension = output_path.split('.')[-1].lower()
                format_mapping = {
                    'csv': largest_rp.to_csv,
                    'parquet': largest_rp.to_parquet,
                }
                save_function = format_mapping.get(file_extension, largest_rp.to_csv)
                save_function(output_path)
                print(f'Risk assessment results saved to {output_path}')

            # return the risk in a concise dataframe
            return largest_rp[
                ["osm_id", "object_type", "geometry"] + list(multi_curves.keys())
            ]

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
    features = data_path / "exposure" / "jamaica-latest.osm.pbf"
    hazard = data_path / "hazard" / "FD_1in1000.tif"
    curves = data_path / "vulnerability" / "curves_osm.csv"
    maxdam = data_path / "vulnerability" / "maxdam_osm.csv"

    # estimate exposure
    asset_types = [
        "main_roads",
        "rail",
        "air",
        "telecom",
        "water_supply",
        "waste_solid",
        "waste_water",
        "education",
        "healthcare",
        "power",
        "gas",
        "oil",
        "buildings",
    ]

    for asset_type in asset_types:
        exposed_features = DamageScanner(hazard, features, curves, maxdam).exposure(
            asset_type=asset_type
        )

        # exposed_features.to_parquet("main_roads.parquet")
        print(exposed_features[["object_type", "coverage", "values"]])

    # #initiate the damage scanner and calculate the damages
    print(
        DamageScanner(hazard, features, curves, maxdam)
        .calculate(asset_type="main_roads")
        .damage.sum()
    )
