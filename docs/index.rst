.. damagescanner documentation master file, created by
   sphinx-quickstart on Sun Jan 27 15:46:05 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for the DamageScanner
=============================================

Python toolkit for direct damage assessments for natural disasters.

Background
----------
This package is (loosely) based on the original DamageScanner, which calculated 
potential flood damages based on inundation depth and land use using depth-damage curves 
in the Netherlands. The DamageScanner was originally developed for the 'Netherlands Later' 
`(Klijn et al., 2007) <https://www.rivm.nl/bibliotheek/digitaaldepot/WL_rapport_Overstromingsrisicos_Nederland.pdf/>`_
. The original land-use classes were based on the Land-Use Scanner in order to evaluate the effect of future land-use change on flood damages. 

This package aims to make this method widely available and for everyone to use. Next to a 
(generalized) function for estimating damages based on rasterdata, it also includes a damage 
assessment function using vector land-use data. 

Even though the method is initially developed for flood damage assessments, it can 
calculate damages for any hazard for which you just require a fragility curve (i.e. a one-dimensional relation). 

Installation
------------
1. Open the python environment in your command prompt or bash in which you want to install this package.
2. Type ``pip install damagescanner`` and it should install itself into your python environment.
3. Now you can import the package like any other package!

OR

1. Clone the repository or download the package on your computer and extract the folder.
2. Go to the DamageScanner folder in your command prompt or bash.
3. Type ``python setup.py install`` and it should install itself into your python environment.
4. Now you can import the package like any other package!


Future work:
------
* Make inputs even more flexible. Catch common errors.
* Make plotting more flexible.
* Develop automated damage assessments using OpenStreetMap data.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   getstarted
   main
   vector
   raster
   risk
   plot
   sensitivity

Indices and tables
==================

* :ref:`search`

