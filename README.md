# DamageScanner: direct damage assessments for natural hazards

<img align="right" width="200" alt="Logo" src="https://raw.githubusercontent.com/ElcoK/DamageScanner/master/doc/ds_logo.png">


[![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/ElcoK/DamageScanner)
[![github license badge](https://img.shields.io/github/license/ElcoK/DamageScanner)](https://github.com/ElcoK/DamageScanner)
[![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)
[![Documentation Status](https://readthedocs.org/projects/damagescanner/badge/?version=latest)](https://damagescanner.readthedocs.io/en/latest/?badge=latest) 
[![PyPI version](https://badge.fury.io/py/damagescanner.svg)](https://badge.fury.io/py/damagescanner) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2551015.svg)](https://doi.org/10.5281/zenodo.2551015) 
[![PyPI - Downloads](https://img.shields.io/pypi/dm/damagescanner?color=yellow&label=Downloads)](https://pypistats.org/packages/damagescanner)


A python toolkit for direct damage assessments for natural hazards. Even though the method is initially developed for flood damage assessments, it can calculate damages for any hazard for which you just require a vulnerability curve (i.e. a one-dimensional relation). 

**Please note:** This package is still in development phase. In case of any problems, or if you have any suggestions for improvements, please raise an *issue*. 

## Background
This package is (loosely) based on the original DamageScanner, which calculated potential flood damages based on inundation depth and land use using depth-damage curves in the Netherlands. The DamageScanner was originally developed for the 'Netherlands Later' project [(Klijn et al., 2007)](https://www.rivm.nl/bibliotheek/digitaaldepot/WL_rapport_Overstromingsrisicos_Nederland.pdf).  The original land-use classes were based on the Land-Use Scanner in order to evaluate the effect of future land-use change on flood damages. 

## Installation
[![workflow pypi badge](https://img.shields.io/pypi/v/damagescanner.svg?colorB=blue)](https://pypi.python.org/project/damagescanner/)

**Requirements:** [NumPy](http://www.numpy.org/), [pandas](https://pandas.pydata.org/), [geopandas](http://geopandas.org/), [matplotlib](https://matplotlib.org/), [rasterio](https://github.com/mapbox/rasterio), [tqdm](https://github.com/tqdm/tqdm), 
[xarray](https://docs.xarray.dev/en/stable/), [pyproj](https://pyproj4.github.io/pyproj/stable/) 


1. Open the python environment in your command prompt or bash in which you want to install this package.
2. Type ``pip install damagescanner`` and it should install itself into your python environment.
3. Now you can import the package like any other package!

OR:

1. Clone the repository or download the package on your computer and extract the folder.
2. Go to the DamageScanner folder in your command prompt or bash.
3. Type ``python setup.py install`` and it should install itself into your python environment.
4. Now you can import the package like any other package!

## Create testing environment
Recommended option is to use a [miniconda](https://conda.io/miniconda.html)
environment to work in for this project, relying on conda to handle some of the
trickier library dependencies.

```bash

# Add conda-forge channel for extra packages
conda config --add channels conda-forge

# Create a conda environment for the project and install packages
conda env create -f environment.yml
activate ds_env

```

## Documentation
[![Documentation Status](https://readthedocs.org/projects/damagescanner/badge/?version=latest)](https://damagescanner.readthedocs.io/en/latest/?badge=latest) 

Please refer to the [ReadTheDocs](http://damagescanner.readthedocs.io/) of this project for the full documentation of all functions. 

## How to cite:
If you use the **DamageScanner** in your work, please cite the package directly:

* Koks. E.E. (2022). DamageScanner: Python tool for natural hazard damage assessments. Zenodo. http://doi.org/10.5281/zenodo.2551015

Here's an example BibTeX entry:

        @misc{damagescannerPython,
              author       = {Koks, E.E.},
              title        = {DamageScanner: Python tool for natural hazard damage assessments},
              year         = 2022,
              doi          = {10.5281/zenodo.2551015},
              url          = {http://doi.org/10.5281/zenodo.2551015}
        }

### License
Copyright (C) 2022 Elco Koks. All versions released under the [MIT license](LICENSE).
