"""DamageScanner - a directe damage assessment toolkit

Raster specific functions
"""

import rasterio
import rasterio.transform
from rasterio.windows import Window


def match_and_load_rasters(raster_in1, raster_in2):
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
    with rasterio.open(raster_in1) as src1, rasterio.open(raster_in2) as src2:
        if src1.crs != src2.crs:
            raise ValueError("Different CRS: CRS must be the same.")
        if src1.res != src2.res:
            raise ValueError("Different resolution: Cell sizes must be the same.")

        top_delta = round((src2.bounds.top - src1.bounds.top) / src1.transform.e)
        bottom_delta = round(
            (src2.bounds.bottom - src1.bounds.bottom) / src1.transform.e
        )
        left_delta = round((src2.bounds.left - src1.bounds.left) / src1.transform.a)
        right_delta = round((src2.bounds.right - src1.bounds.right) / src1.transform.a)

        data1 = src1.read(
            1,
            window=Window(
                col_off=left_delta,
                row_off=top_delta,
                width=src1.width - left_delta + right_delta,
                height=src1.height - top_delta + bottom_delta,
            ),
        )
        data2 = src2.read(
            1,
            window=Window(
                col_off=abs(min(left_delta, 0)),
                row_off=abs(min(top_delta, 0)),
                width=max(src1.width, src2.width) - abs(left_delta) - abs(right_delta),
                height=max(src1.height, src2.height)
                - abs(top_delta)
                - abs(bottom_delta),
            ),
        )
        transform = rasterio.Affine(
            src1.transform.a,
            src1.transform.b,
            src1.transform.c + src1.transform.a * max(left_delta, 0),
            src1.transform.d,
            src1.transform.e,
            src1.transform.f + src1.transform.e * max(top_delta, 0),
        )

    return data1, data2, transform
