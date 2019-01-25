"""DamageScanner - a directe damage assessment toolkit
"""
import pkg_resources

# Define what is accessible directly on snkit, when a client writes::
#   from snkit import Network
from damagescanner.main import DamageScanner

__author__ = "Elco Koks"
__copyright__ = "Elco Koks"
__license__ = "GNU"


try:
    __version__ = pkg_resources.get_distribution(__name__).version
except pkg_resources.DistributionNotFound:
    __version__ = 'unknown'


# Define what should be imported as * when a client writes::
#   from snkit import *
__all__ = ['DamageScanner']