from shapely.geometry import Polygon
from damagescanner.compute import download_osm


def test_download_osm():
    polygon = Polygon(
        [
            (5.7920340896, 50.6832293511),
            (6.2064643502, 50.6832293511),
            (6.2064643502, 50.84305974),
            (5.7920340896, 50.84305974),
            (5.7920340896, 50.6832293511),
        ]
    )
    file_paths = download_osm(polygon, overwrite=True)
    assert len(file_paths) == 3
    for file_path in file_paths:
        assert file_path.exists()
