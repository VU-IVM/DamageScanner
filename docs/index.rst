.. DamageScanner documentation master file, created by
   sphinx-quickstart on Sat Feb 11 13:31:14 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==============================================================
DamageScanner: direct damage assessments for natural disasters
==============================================================

A python toolkit for direct damage assessments for natural disasters. Even though the method is initially developed for flood damage assessments, it can calculate damages for any hazard for which you just require a vulnerability curve (i.e. a one-dimensional relation). 

**Please note:** This package is still in development phase. In case of any problems, or if you have any suggestions for improvements, please raise an *issue*. 

Background
##########

This package is (loosely) based on the original DamageScanner, which calculated potential flood damages based on inundation depth and land use using depth-damage curves in the Netherlands. The DamageScanner was originally developed for the 'Netherlands Later' project `(Klijn et al., 2007) <https://www.rivm.nl/bibliotheek/digitaaldepot/WL_rapport_Overstromingsrisicos_Nederland.pdf>`_.  The original land-use classes were based on the Land-Use Scanner in order to evaluate the effect of future land-use change on flood damages. 

Installation
############

**Requirements:** `NumPy <http://www.numpy.org/>`_, `pandas <https://pandas.pydata.org/>`_, `geopandas <http://geopandas.org/>`_, `matplotlib <https://matplotlib.org/>`_, `rasterio <https://github.com/mapbox/rasterio>`_, `tqdm <https://github.com/tqdm/tqdm>`_, `xarray <https://docs.xarray.dev/en/stable/>`_, `pyproj <https://pyproj4.github.io/pyproj/stable/>`_ 


1. Open the python environment in your command prompt or bash in which you want to install this package.
2. Type ``pip install damagescanner`` and it should install itself into your python environment.
3. Now you can import the package like any other package!

OR:

1. Clone the repository or download the package on your computer and extract the folder.
2. Go to the DamageScanner folder in your command prompt or bash.
3. Type ``python setup.py install`` and it should install itself into your python environment.
4. Now you can import the package like any other package!

How to cite
###########

If you use the **DamageScanner** in your work, please cite the package directly:

* Koks. E.E. (2022). DamageScanner: Python tool for disaster damage assessments. Zenodo. http://doi.org/10.5281/zenodo.2551015

License
#######

Copyright (C) 2022 Elco Koks. All versions released under the **MIT license**.




.. toctree::
   :maxdepth: 1
   :caption: Contents:

   index
   getstarted

