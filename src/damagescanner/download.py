"""
This file is part of OSM-flex.
Copyright (C) 2023 OSM-flex contributors listed in AUTHORS.
OSM-flex is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, version 3.
OSM-flex is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
-----
downloading functions
"""

import logging
from pathlib import Path
import urllib.request
from urllib.parse import urljoin
from damagescanner.config import DICT_GEOFABRIK, GEOFABRIK_URL, PLANET_URL, OSM_DATA_DIR

LOGGER = logging.getLogger(__name__)

# =============================================================================
#  DOWNLOAD METHODS
# =============================================================================

def _create_gf_download_url(iso3, file_format):
    """
    create string with download-url from geofabrik

    Parameters
    ----------
    iso3 : str
        ISO3 code of country to download
    file_format : str
        Format in which file should be downloaded; ESRI Shapefiles ('shp')
        or osm-Protocolbuffer Binary Format ('pbf')

    Returns
    -------
    url: str
        Geofabrik download URL for the requested country.

    See also
    --------
    DICT_GEOFABRIK for exceptions / special regions.
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

def _download_file(download_url: str, filepath: Path, overwrite: bool = True):
    """Download a file located at an URL to a local file path

    Parameters
    ----------
    download_url : str
        URL of the file to download
    filepath : str or Path
        Local file path to store the file
    overwrite : bool, optional
        Overwrite existing files. If ``False``, the download will be skipped for
        existing files. Defaults to ``True``.
    """
    if not Path(filepath).is_file() or overwrite:
        LOGGER.info(f"Download file: {filepath}")
        urllib.request.urlretrieve(download_url, filepath)
    else:
        LOGGER.info(f"Skip existing file: {filepath}")

# TODO: decide whether to issue warnings for multi-country files

def get_country_geofabrik(iso3, file_format='pbf', save_path=OSM_DATA_DIR,
                          overwrite=False):
    """
    Download country files with all OSM map info from the provider
    Geofabrik.de.

    Parameters
    ----------
    iso3 : str
        ISO3 code of country to download
        Exceptions: Russia is divided into European and Asian part
        ('RUS-E', 'RUS-A'), Canary Islands are 'IC'.
    file_format : str
        Format in which file should be downloaded; options are
        ESRI Shapefiles (shp), which can easily be loaded into gdfs,
        or osm-Protocolbuffer Binary Format (pbf), which is smaller in
        size, but has a more complicated query syntax to load (functions
        are provided in the OSMFileQuery class).
    save_path : str or pathlib.Path
        Folder in which to save the file

    Returns
    -------
    filepath : Path
        The path to the downloaded file (``save_path`` + the Geofabrik filename)

    See also
    --------
    DICT_GEOFABRIK for exceptions / special regions.
    """

    download_url = _create_gf_download_url(iso3, file_format)
    filepath = Path(save_path, Path(download_url).name)
    filepath.parent.mkdir(exist_ok=True, parents=True)
    _download_file(download_url, filepath, overwrite)

    return filepath

# TODO: allow for several spelling options like "Central America", "Australia", ...
def get_region_geofabrik(region, save_path=OSM_DATA_DIR, overwrite=False):
    """
    Download regions files with all OSM map info from the provider
    Geofabrik.de
    
    Parameters
    ----------
    region: str
        one of Africa, Antarctica, Asia, Australia-and-Oceania,
        Central-America, Europe, North-America, South-America
    save_path : str or pathlib.Path
        Folder in which to save the file

    Returns
    -------
    filepath : Path
        The path to the downloaded file
    """

    download_url =  f'{GEOFABRIK_URL}{region.lower()}-latest.osm.pbf'
    filepath = Path(save_path, Path(download_url).name)
    filepath.parent.mkdir(exist_ok=True, parents=True)
    _download_file(download_url, filepath, overwrite)

    return filepath


def get_planet_file(save_path=Path(OSM_DATA_DIR,'planet-latest.osm.pbf'),
                    overwrite=False):
    """
    Download the entire planet file from the OSM server (ca. 60 GB).

    Parameters
    ----------
    save_path : str or pathlib.Path
        The path to store the file.

    Returns
    -------
    save_path : Path
        The path to the downloaded file. Returned for consistency with other download
        functions.
    """
    _download_file(PLANET_URL, save_path, overwrite)
    return Path(save_path)
