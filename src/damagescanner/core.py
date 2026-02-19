"""DamageScanner - a directe damage assessment toolkit.

Copyright (C) 2023 Elco Koks. All versions released under the MIT license.
"""

from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from tqdm import tqdm

from damagescanner.raster import RasterScanner
from damagescanner.vector import VectorExposure, VectorScanner


class DamageScanner(object):
    """
    DamageScanner - a direct damage assessment toolkit.

    This class provides tools to assess direct physical damage from hazard events,
    using raster or vector exposure data and vulnerability curves.
    It supports both single-hazard footprints and risk-based multi-scenario assessments.
    """

    def __init__(
        self,
        hazard_data: str | Path | xr.DataArray | xr.Dataset,
        feature_data: str
        | Path
        | pd.DataFrame
        | gpd.GeoDataFrame
        | xr.DataArray
        | xr.Dataset,
        curves: str | Path | pd.DataFrame,
        maxdam: str | Path | pd.DataFrame,
    ) -> None:
        """
        Initialize the DamageScanner class with hazard, exposure, curve, and max damage data.

        Args:
            hazard_data: Path to raster hazard file or xarray object.
            feature_data: Exposure data, either raster or vector.
            curves: Vulnerability curves as DataFrame or CSV file path.
            maxdam: Maximum damage values per asset type.

        Raises:
            ImportError: If the input data formats are not supported or if required data is missing.
            ImportWarning: If the curves or maxdam data is not in the expected format.
        """
        # Convert the input to a Path object if it is a string
        if isinstance(feature_data, str):
            feature_data = Path(feature_data)

        if isinstance(hazard_data, str):
            hazard_data = Path(hazard_data)

        if isinstance(curves, str):
            curves = Path(curves)

        if isinstance(maxdam, str):
            maxdam = Path(maxdam)

        # Collect the input data
        self.hazard_data = hazard_data
        self.feature_data = feature_data
        self.curves = curves
        self.maxdam = maxdam

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
                    "The exposure data should be a a shapefile, geopackage, \
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

    def exposure(
        self,
        disable_progress: bool = False,
        output_path: Path | str | None = None,
        **kwargs: Any,
    ) -> gpd.GeoDataFrame | xr.DataArray:
        """
        Run the exposure analysis to identify features affected by the hazard footprint.

        This method analyzes the input data to determine which features are exposed to a hazard.
        It supports both raster and vector data types for the exposure analysis. If the data type
        is vector, additional keyword arguments can be specified to customize the analysis.

        Args:
            disable_progress: If True, disables progress bars during processing.
            output_path: The file path to save the exposure data. The file format
                is determined by the file extension (e.g., '.parquet', '.csv', '.gpkg', '.shp').
                If None, the data is not saved.
            **kwargs: Optional keyword arguments:
                - asset_type (str): The type of asset to evaluate (only for vector data).
                - return_full (bool): Whether to return all features, even those with no hazard intersection.
                    Defaults to True.

        Returns:
            A GeoDataFrame containing the affected assets if the input data is vector,
            or an xarray.DataArray if the input data is raster.

        Notes:
            - If `output_path` is provided, the method saves the exposure data to the specified path.
            - The file format is inferred from the file extension of `output_path`. If the extension
            is not recognized, the data is saved as a Parquet file by default.

        Raises:
            ImportError: If the input data formats are not supported or if required data is missing.
        """
        if isinstance(output_path, str):
            output_path = Path(output_path)

        if self.assessment_type == "raster":
            return xr.open_rasterio(self.exposure_data)  # ty:ignore[unresolved-attribute]

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
                return_full=kwargs.get("return_full", True),
            )[0]

            # save output when exposed assets are empty
            if output_path:
                # Determine the file format based on the file extension
                file_extension = output_path.suffix[1:].lower()
                format_mapping = {
                    "parquet": exposed_assets.to_parquet,
                    "csv": exposed_assets.to_csv,
                    "gpkg": exposed_assets.to_file,
                    "shp": exposed_assets.to_file,
                }

                # Default to parquet if the extension is not recognized
                save_function = format_mapping.get(
                    file_extension, exposed_assets.to_parquet
                )
                save_function(output_path)

                print(f"Exposure data saved to {output_path}")

            return exposed_assets
        else:
            raise ImportError("Please prepare the input data first")

    def calculate(
        self,
        disable_progress: bool = False,
        output_path: str | None = None,
        **kwargs: Any,
    ) -> (
        pd.DataFrame
        | tuple[pd.DataFrame, np.ndarray, np.ndarray, xr.Dataset | np.ndarray]
    ):
        """
        Perform a damage calculation using the provided inputs.

        Applies vulnerability curves and maximum damage values to the exposed features
        or raster grid to calculate expected damage.

        Args:
            disable_progress: If True, disables progress bars.
            output_path: Path to save the calculation results. The file format
                is determined by the file extension (e.g., '.csv', '.parquet' for vector data,
                '.tif' for raster data). If None, the data is not saved.
            **kwargs:
                object_col (str): Column name containing object/landuse types.
                Defaults to "object_type". (vector only)
                asset_type (str, optional): Infrastructure class to evaluate.
                multi_curves (dict, optional): Mapping of asset types to curve sets.
                return_full (bool): Whether to return all features, even those with no hazard intersection.
                    Defaults to True. (vector only)

        Returns:
            For vector: DataFrame with damages.
            For raster: tuple of (damage_df, damagemap, landuse, hazard).

        Raises:
            ImportError: If the input data formats are not supported or if required data is missing.
            ValueError: If the geometry type is unsupported or if object types in the data are not found in the curves.
        """
        if not hasattr(self, "assessment_type"):
            raise ImportError("Please prepare the input data first")

        if self.assessment_type == "raster":
            ds_results = RasterScanner(
                exposure_file=self.feature_data,
                hazard_file=self.hazard_data,
                curve_path=self.curves,
                maxdam_path=self.maxdam,
            )

            damage_df = ds_results[0]
            damagemap = ds_results[1]

            # Extract CRS and transform from feature_data
            if isinstance(self.feature_data, (str, Path)):
                feature_path = Path(self.feature_data)  # Ensure it's a Path
                if feature_path.suffix == ".nc":
                    # Open with xarray if it's a NetCDF file
                    feature_data_xr = xr.open_dataset(feature_path)
                    crs = feature_data_xr.rio.crs
                    transform = feature_data_xr.rio.transform()
                else:
                    # Open with rasterio for other raster formats
                    with rasterio.open(feature_path) as src:
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
                object_col=kwargs.get("object_col", "object_type"),
                disable_progress=disable_progress,
                return_full=kwargs.get("return_full", True),
            )

            # For vector data, CRS and transform are not directly applicable
            crs = None
            transform = None

        if output_path:
            file_extension = output_path.split(".")[-1].lower()
            if self.assessment_type == "vector":
                format_mapping = {
                    "csv": damage_df.to_csv,
                    "parquet": damage_df.to_parquet,
                }
                save_function = format_mapping.get(file_extension, damage_df.to_csv)
                save_function(output_path, **kwargs)

            elif self.assessment_type == "raster":
                # Save the damage_df as CSV
                damage_df_path = output_path.replace(
                    f".{file_extension}", "_damages.csv"
                )
                damage_df.to_csv(damage_df_path)
                print(f"Damage summary saved to {damage_df_path}")

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
                print(f"Damage map saved to {dmap_fn}")

        return damage_df if self.assessment_type == "vector" else ds_results

    def risk(
        self,
        hazard_dict: dict[int | float, str | Path],
        output_path: str | None = None,
        **kwargs: Any,
    ) -> pd.DataFrame | None:
        """
        Perform a risk assessment across multiple hazard return periods.

        Integrates damages from each return period and computes expected annual damages.
        Supports both single and multi-curve inputs for infrastructure types.

        Args:
            hazard_dict: Dictionary mapping return periods to hazard raster paths.
            output_path: Path to save the risk assessment results. The file format
                is determined by the file extension (e.g., '.csv', '.parquet'). If None, the data is not saved.
            **kwargs:
                asset_type (str, optional): Infrastructure class to evaluate.
                multi_curves (dict, optional): Mapping of asset types to curve sets.

        Returns:
            A GeoDataFrame with risk values for each asset, or None if no results.
        """
        RP_list = list(hazard_dict.keys())

        risk = {}
        for key, hazard_map in tqdm(
            hazard_dict.items(), total=len(hazard_dict), desc="Risk Calculation"
        ):
            if self.assessment_type == "raster":
                risk[key] = DamageScanner(
                    hazard_map, self.feature_data, self.curves, self.maxdam
                ).calculate(disable_progress=True, **kwargs)[0]
            else:
                risk[key] = DamageScanner(
                    hazard_map, self.feature_data, self.curves, self.maxdam
                ).calculate(disable_progress=True, **kwargs)

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
                file_extension = output_path.split(".")[-1].lower()
                format_mapping = {
                    "csv": largest_rp.to_csv,
                    "parquet": largest_rp.to_parquet,
                }
                save_function = format_mapping.get(file_extension, largest_rp.to_csv)
                save_function(output_path)
                print(f"Risk assessment results saved to {output_path}")

            # return the risk in a concise dataframe
            if self.assessment_type == "raster":
                return largest_rp
            else:
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
                    lambda x: np.trapezoid(y=x[RP_list][::-1], x=RPS[::-1]), axis=1
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
                file_extension = output_path.split(".")[-1].lower()
                format_mapping = {
                    "csv": largest_rp.to_csv,
                    "parquet": largest_rp.to_parquet,
                }
                save_function = format_mapping.get(file_extension, largest_rp.to_csv)
                save_function(output_path)
                print(f"Risk assessment results saved to {output_path}")

            # return the risk in a concise dataframe
            return largest_rp[
                ["osm_id", "object_type", "geometry"] + list(multi_curves.keys())
            ]
