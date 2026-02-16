# DamageScanner: direct damage assessments for natural hazards

<img align="right" width="200" alt="Logo" src="https://raw.githubusercontent.com/ElcoK/DamageScanner/main/docs/images/logo-dark.png">


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
To use `DamageScanner` in your project:

### Using `uv` (Recommended)
```bash
uv add damagescanner
```

### Using `pip`
```bash
pip install damagescanner
```

## Development & Testing
To set up a local environment for development or to run tests:

### Using `uv` (Recommended)
[uv](https://github.com/astral-sh/uv) is an extremely fast Python package manager and is the preferred way to set up the development environment.

```bash
# Clone the repository
git clone https://github.com/VU-IVM/DamageScanner.git
cd DamageScanner

# Create a virtual environment and install all optional dependencies
uv sync --all-groups
```

### Using Miniconda
If you prefer [Miniconda](https://docs.conda.io/en/latest/miniconda.html), use the provided `environment.yml` file:

```bash
# Add conda-forge channel for extra packages
conda config --add channels conda-forge

# Create environment and activate
conda env create -f environment.yml
conda activate ds-test
```

## Documentation
Please refer to the [ReadTheDocs](https://vu-ivm.github.io/DamageScanner/) of this project for the full documentation of all functions. 

## How to cite
If you use the **DamageScanner** in your work, please cite the package directly:

* Koks. E.E. (2022). DamageScanner: Python tool for natural hazard damage assessments. Zenodo. http://doi.org/10.5281/zenodo.2551015

Here's an example BibTeX entry:

```
@misc{damagescannerPython,
      author       = {Koks, E.E.},
      title        = {DamageScanner: Python tool for natural hazard damage assessments},
      year         = 2022,
      doi          = {10.5281/zenodo.2551015},
      url          = {http://doi.org/10.5281/zenodo.2551015}
}
```

## License
Copyright (C) 2022 Elco Koks. All versions released under the [MIT license](LICENSE).
