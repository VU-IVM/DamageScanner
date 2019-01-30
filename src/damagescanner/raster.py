"""DamageScanner - a directe damage assessment toolkit

Raster specific functions
"""

from osgeo import gdal


def match(raster_in1,raster_in2):
    """
    In case of a mismatch between two rasters, return only the intersecting parts.
    
    Code adapted from http://sciience.tumblr.com/post/101722591382/finding-the-georeferenced-intersection-between-two
    
    Arguments:
        *raster_in1* : One of the two rasters to be clipped to the overlapping extent.
        
        *raster_in2* : One of the two rasters to be clipped to the overlapping extent.
        
    Returns:
        *array1* : Numpy Array of raster1
        
        *array2* : Numpy Array of raster2
        
        *intersection* : Bounding box of overlapping part
    
    """
    # Read rasterdata
    raster1 = gdal.Open(raster_in1)
    raster2 = gdal.Open(raster_in2)
    
    # load data
    band1 = raster1.GetRasterBand(1)
    band2 = raster2.GetRasterBand(1)
    gt1 = raster1.GetGeoTransform()
    gt2 = raster2.GetGeoTransform()
    
    # find each image's bounding box
    # r1 has left, top, right, bottom of dataset's bounds in geospatial coordinates.
    r1 = [gt1[0], gt1[3], gt1[0] + (gt1[1] * raster1.RasterXSize), gt1[3] + (gt1[5] * raster1.RasterYSize)]
    r2 = [gt2[0], gt2[3], gt2[0] + (gt2[1] * raster2.RasterXSize), gt2[3] + (gt2[5] * raster2.RasterYSize)]
    print('\t1 bounding box: %s' % str(r1))
    print('\t2 bounding box: %s' % str(r2))
    
    # find intersection between bounding boxes
    intersection = [max(r1[0], r2[0]), min(r1[1], r2[1]), min(r1[2], r2[2]), max(r1[3], r2[3])]
    if r1 != r2:
        print('\t** different bounding boxes **')
        # check for any overlap at all...
        if (intersection[2] < intersection[0]) or (intersection[1] < intersection[3]):
            intersection = None
            print('\t***no overlap***')
            return
        else:
            print('\tintersection:',intersection)
            left1 = int(round((intersection[0]-r1[0])/gt1[1])) # difference divided by pixel dimension
            top1 = int(round((intersection[1]-r1[1])/gt1[5]))
            col1 = int(round((intersection[2]-r1[0])/gt1[1])) - left1 # difference minus offset left
            row1 = int(round((intersection[3]-r1[1])/gt1[5])) - top1
            
            left2 = int(round((intersection[0]-r2[0])/gt2[1])) # difference divided by pixel dimension
            top2 = int(round((intersection[1]-r2[1])/gt2[5]))
            col2 = int(round((intersection[2]-r2[0])/gt2[1])) - left2 # difference minus new left offset
            row2 = int(round((intersection[3]-r2[1])/gt2[5])) - top2
            
            #print '\tcol1:',col1,'row1:',row1,'col2:',col2,'row2:',row2
            if col1 != col2 or row1 != row2:
                print("ERROR: Columns and rows still do not match! ***")
            # these arrays should now have the same spatial geometry though NaNs may differ
            array1 = band1.ReadAsArray(left1,top1,col1,row1)
            array2 = band2.ReadAsArray(left2,top2,col2,row2)

    else: # same dimensions from the get go
        col1 = raster1.RasterXSize # = col2
        row1 = raster1.RasterYSize # = row2
        array1 = band1.ReadAsArray()
        array2 = band2.ReadAsArray()
        
    return array1, array2, intersection