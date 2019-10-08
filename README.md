[![Build Status](https://travis-ci.com/ElcoK/DamageScanner.svg?branch=master)](https://travis-ci.com/ElcoK/DamageScanner) [![Documentation Status](https://readthedocs.org/projects/damagescanner/badge/?version=latest)](https://damagescanner.readthedocs.io/en/latest/?badge=latest) [![PyPI version](https://badge.fury.io/py/damagescanner.svg)](https://badge.fury.io/py/damagescanner) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2551015.svg)](https://doi.org/10.5281/zenodo.2551015) 
[![PyPI - Downloads](https://img.shields.io/pypi/dm/damagescanner?color=yellow&label=Downloads)](https://pypistats.org/packages/damagescanner)

# DamageScanner 

Python toolkit for direct damage assessments for natural disasters.

Please refer to the [ReadTheDocs](http://damagescanner.readthedocs.io/) of this project for the full documentation of all functions. 

**Requirements:** [NumPy](http://www.numpy.org/), [pandas](https://pandas.pydata.org/), [geopandas](http://geopandas.org/), [matplotlib](https://matplotlib.org/), [rasterio](https://github.com/mapbox/rasterio), [tqdm](https://github.com/tqdm/tqdm) 

**Please note:** This package is still in development phase. In case of any problems, or if you have any suggestions for improvements, please raise an *issue*. 

## Background
This package is (loosely) based on the original DamageScanner, which calculated potential flood damages based on inundation depth and land use using depth-damage curves in the Netherlands. The DamageScanner was originally developed for the 'Netherlands Later' project [(Klijn et al., 2007)](https://www.rivm.nl/bibliotheek/digitaaldepot/WL_rapport_Overstromingsrisicos_Nederland.pdf).  The original land-use classes were based on the Land-Use Scanner in order to evaluate the effect of future land-use change on flood damages. 

This package aims to make this method widely available and for everyone to use. Next to a (generalized) function for estimating damages based on raster data, it also includes a damage assessment function using vector land-use data. 

Even though the method is initially developed for flood damage assessments, it can calculate damages for any hazard for which you just require a fragility curve (i.e. a one-dimensional relation). 

## Installation

1. Open the python environment in your command prompt or bash in which you want to install this package.
2. Type ``pip install damagescanner`` and it should install itself into your python environment.
3. Now you can import the package like any other package!

OR:

1. Clone the repository or download the package on your computer and extract the folder.
2. Go to the DamageScanner folder in your command prompt or bash.
3. Type ``python setup.py install`` and it should install itself into your python environment.
4. Now you can import the package like any other package!

## To-do:
* Make inputs for both the RasterScanner and VectorScanner more flexible. Catch common errors.
* Make plotting functions more flexible.
* Develop automated damage assessments using OpenStreetMap data.

## How to cite:
If you use the **DamageScanner** in your work, please cite the package directly:

* Koks. E.E. (2019). DamageScanner: Python tool for disaster damage assessments. Zenodo. http://doi.org/10.5281/zenodo.2551015

Here's an example BibTeX entry:

        @misc{damagescannerPython,
              author       = {Koks, E.E.},
              title        = {DamageScanner: Python tool for disaster damage assessments},
              year         = 2019,
              doi          = {10.5281/zenodo.2551015},
              url          = {http://doi.org/10.5281/zenodo.2551015}
        }

### License
Copyright (C) 2019 Elco Koks. All versions released under the [MIT license](LICENSE).


![IVM](http://ivm.vu.nl/en/Images/IVM_logo_rgb2_tcm234-851594.svg)
