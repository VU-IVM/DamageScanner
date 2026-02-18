"""This file is part of OSM-flex.

Copyright (C) 2023 OSM-flex contributors listed in AUTHORS.
OSM-flex is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, version 3.
OSM-flex is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
-----
Downloading utilities for OSM-flex and DamageScanner.

This module provides methods to download OSM data (in `.pbf` or `.shp.zip` format)
from the Geofabrik or Planet OSM mirrors. It handles ISO3-to-region mapping,
file saving, and logging.

Functions:
    - get_country_geofabrik: Download OSM data for a specific country using ISO3.
    - get_region_geofabrik: Download OSM data for an entire region (e.g., Europe).
    - get_planet_file: Download the full global OpenStreetMap planet file (~60 GB).

See Also:
    - DICT_GEOFABRIK: ISO3-to-region mapping
    - GEOFABRIK_URL / PLANET_URL: Source locations
"""

import logging
import urllib.request
from pathlib import Path
from urllib.parse import urljoin

from damagescanner.config import DICT_GEOFABRIK, GEOFABRIK_URL, OSM_DATA_DIR, PLANET_URL

LOGGER = logging.getLogger(__name__)

# =============================================================================
#  DOWNLOAD METHODS
# =============================================================================


def _create_gf_download_url(iso3: str, file_format: str) -> str:
    """Generates a download URL for a country from Geofabrik.

    Args:
        iso3: ISO3 code of the country (e.g., 'NLD' for Netherlands).
        file_format: File format to download; either 'shp' or 'pbf'.

    Returns:
        The full Geofabrik download URL.

    Raises:
        KeyError: If the ISO3 code is not found in DICT_GEOFABRIK.
        NotImplementedError: If an unsupported file format is specified.

    See Also:
        DICT_GEOFABRIK: Mapping of ISO3 to continent and country names.
    """
    # Retrieve Geofabrik definitions
    try:
        continent, country = DICT_GEOFABRIK[iso3]

    except KeyError as err:
        if iso3 == "RUS":
            raise KeyError(
                "Russia comes in two files. Please specify either RUS-A for the Asian "
                "or RUS-E for the European part."
            ) from err
        raise KeyError(
            "The provided iso3 seems not to be available on Geofabrik.de. You can clip "
            "it from the planet file or an adequate regional file, instead. See the OSM-flex "
            "methods in the clip module for this."
        ) from err

    # Set file extension
    if file_format == "shp":
        ext = "-free.shp.zip"
    elif file_format == "pbf":
        ext = ".osm.pbf"
    else:
        raise NotImplementedError(
            f"Invalid file format '{file_format}'. Please choose one of [shp, pbf]"
        )

    # Join to URL
    return urljoin(GEOFABRIK_URL, f"{continent}/{country}-latest{ext}")


def _download_file(download_url: str, filepath: Path, overwrite: bool = True) -> None:
    """Downloads a file from a URL to a local path.

    Args:
        download_url: The URL of the file to download.
        filepath: Destination path for the downloaded file.
        overwrite: If True, overwrite existing file.
    """
    if not Path(filepath).is_file() or overwrite:
        LOGGER.info(f"Download file: {filepath}")
        urllib.request.urlretrieve(download_url, filepath)
    else:
        LOGGER.info(f"Skip existing file: {filepath}")


def get_country_geofabrik(
    iso3: str,
    file_format: str = "pbf",
    save_path: str | Path = OSM_DATA_DIR,
    overwrite: bool = False,
) -> Path:
    """Downloads an OSM dataset for a country from Geofabrik.

    Args:
        iso3: ISO3 country code (e.g., 'NLD').
        file_format: Either 'pbf' or 'shp'.
        save_path: Directory to save the file.
        overwrite: If True, overwrite existing file.

    Returns:
        The local path to the downloaded file.

    See Also:
        DICT_GEOFABRIK: ISO3-to-Geofabrik region mapping.

    """
    download_url = _create_gf_download_url(iso3, file_format)
    filepath = Path(save_path, Path(download_url).name)
    filepath.parent.mkdir(exist_ok=True, parents=True)
    _download_file(download_url, filepath, overwrite)

    return filepath


# TODO: allow for several spelling options like "Central America", "Australia", ...
def get_region_geofabrik(
    region: str, save_path: str | Path = OSM_DATA_DIR, overwrite: bool = False
) -> Path:
    """Downloads an OSM dataset for an entire region from Geofabrik.

    Args:
        region: Name of the region (e.g., 'Europe', 'Africa', 'Asia').
        save_path: Directory to save the file.
        overwrite: If True, overwrite existing file.

    Returns:
        The local path to the downloaded file.
    """
    download_url = f"{GEOFABRIK_URL}{region.lower()}-latest.osm.pbf"
    filepath = Path(save_path, Path(download_url).name)
    filepath.parent.mkdir(exist_ok=True, parents=True)
    _download_file(download_url, filepath, overwrite)

    return filepath


def get_planet_file(
    save_path: str | Path = Path(OSM_DATA_DIR, "planet-latest.osm.pbf"),
    overwrite: bool = False,
) -> Path:
    """Downloads the full OSM planet file (~60 GB) from planet.openstreetmap.org.

    Args:
        save_path: Destination path for the file.
        overwrite: If True, overwrite existing file.

    Returns:
        The local path to the downloaded planet file.
    """
    _download_file(PLANET_URL, save_path, overwrite)
    return Path(save_path)
