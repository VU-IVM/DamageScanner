"""Test core objects/concepts"""

# pylint: disable=C0103
import numpy
from damagescanner import RasterScanner

data_path = ".."


def test_RasterScanner():
    # dummy inundation array
    inundation = numpy.array(
        [
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0.5, 1, 1],
            [0, 0, 0, 1, 2, 1],
            [1, 0.25, 0, 0.25, 2, 0.5],
            [0, 0.5, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
        ]
    )

    # dummy landuse array
    landuse = numpy.array(
        [
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1, 4, 4, 4],
            [2, 2, 2, 2, 2, 2],
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 4, 4],
            [3, 3, 3, 1, 1, 1],
        ]
    )

    # dummy maxdam array
    maxdam = numpy.array([[1, 50], [2, 300], [3, 600], [4, 200]])

    # dummy curves array
    curves = numpy.array(
        [
            [0, 0, 0, 0, 0],
            [0.5, 0.2, 0.2, 0.2, 0.2],
            [1.5, 0.4, 0.4, 0.4, 0.4],
            [2, 0.6, 0.6, 0.6, 0.6],
            [2.5, 0.8, 0.8, 0.8, 0.8],
            [3, 1, 1, 1, 1],
        ]
    )

    assert (
        RasterScanner(landuse, inundation, curves, maxdam, cellsize=100)[0]
        .sum()
        .values[0]
        == 59500
    )
