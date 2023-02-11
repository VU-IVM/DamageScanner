
=======================================
Getting started 
=======================================

After installation, the tool works fairly straightforward. 

If you want to do a raster-based calculation, use the **RasterScanner**. If you want 
to do a vector-based calculation, use the **VectorScanner**.

Please find below examples using the example data stored on the GitHub.


RasterScanner
-------------
.. code-block:: python
    :linenos:
 
    import os
    
    # import the RasterScanner
    from damagescanner.core import RasterScanner
    
    # set paths to the data
    inun_map = os.path.join(data_path,'data','inundation','inundation_map.tif')
    landuse_map = os.path.join(data_path,'data','landuse','landuse_map.tif')
    curve_path = os.path.join(data_path,'data','curves','curves.csv')
    maxdam_path = os.path.join(data_path,'data','curves','maxdam.csv')
        
    # run the RasterScanner and return a pandas DataFrame with loss per land-use class
    loss_df = RasterScanner(landuse_map,inun_map,curve_path,maxdam_path)[0]


VectorScanner
-------------
.. code-block:: python
    :linenos:
 
    # import necessary packages
    import os
    import numpy
    import pandas 
       
    # import the RasterScanner
    from damagescanner.core import VectorScanner
    
    # set paths to the data
    inun_map = os.path.join(data_path,'data','inundation','inundation_map.tif')
    landuse_map = os.path.join(data_path,'data','landuse','landuse.shp')

    # Create maximum damage dictionary
    maxdam = {"grass":5,
        "forest":10,
        "orchard":50,
        "residential":200,
        "industrial":300,
        "retail":300,
        "farmland":10,
        "cemetery":15,
        "construction":10,
        "meadow":5,
        "farmyard":5,
        "scrub":5,
        "allotments":10,
        "reservoir":5,
        "static_caravan":100,
        "commercial":300}
        
    # Create some dummy curves that will match the land-use classes
    curves = numpy.array(
            [[0,0],
            [50,0.2],
            [100,0.4],
            [150,0.6],
            [200,0.8],
            [250,1]])  
    
    curves = numpy.concatenate((curves,
                                numpy.transpose(numpy.array([curves[:,1]]*(len(maxdam)-1)))),
                               axis=1)
    
    curves = pandas.DataFrame(curves)
    curves.columns = ['depth']+list(maxdam.keys())
    curves.set_index('depth',inplace=True)    

    # run the VectorScanner and return the landuse map with damage values
    loss_df = VectorScanner(landuse,inun_map,curves,maxdam)
