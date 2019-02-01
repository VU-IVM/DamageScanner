"""DamageScanner - a directe damage assessment toolkit

Copyright (C) 2019 Elco Koks. All versions released under the MIT license.
"""
import pkg_resources

__author__ = "Elco Koks"
__copyright__ = "Elco Koks"
__license__ = "MIT"

try:
    __version__ = pkg_resources.get_distribution(__name__).version
except Exception:
    __version__ = 'unknown'
    
__all__ = ['core','plot','vector','raster','sensitivity','utils']
