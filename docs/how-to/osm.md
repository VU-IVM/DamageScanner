# Coupling with OpenStreetMap (OSM)

One of the powerful features of `DamageScanner` is the ability to work with **OpenStreetMap (OSM)** data as exposure input. This enables fast and global-scale exposure mapping, using freely available, crowdsourced infrastructure data.

This page outlines how to prepare and use OSM data for damage assessment workflows, using `.osm.pbf` files directly.

---

## 🌍 Why Use OSM Data?

- Freely available, global coverage
- Includes many infrastructure types (roads, buildings, energy, healthcare, etc.)
- Constantly updated by the global community
- Works well with vector-based workflows in `DamageScanner`

> ⚠️ Note: OSM data quality and completeness can vary by location. Always validate your input.

---

## ⚙️ How DamageScanner Works with OSM

DamageScanner includes built-in support for working directly with `.osm.pbf` files. The tool provides utilities to:

- Download `.osm.pbf` files from Geofabrik
- Extract specific infrastructure classes
- Clean and preprocess the data into usable `GeoDataFrame` or `.gpkg`

These functions are implemented in:
- [`osm.py`](https://github.com/VU-IVM/DamageScanner/blob/installation/src/damagescanner/osm.py)
- [`download.py`](https://github.com/VU-IVM/DamageScanner/blob/installation/src/damagescanner/download.py)

> ✅ This approach allows for highly automated and reproducible integration of OSM infrastructure.

You can extract features like this:

```python
from damagescanner.osm import read_osm_data

# Read critical infrastructure data from an OSM .pbf file
osm_path = "data/netherlands-latest.osm.pbf"
features = read_osm_data(osm_path, asset_type="road")
```

The `asset_type` must be one of the predefined keys such as `road`, `rail`, `power`, `healthcare`, `education`, etc. The output is a cleaned GeoDataFrame with valid geometries and an `object_type` column.

---

## 🧼 Manual Tagging (Optional)

If you're preparing custom OSM files manually (e.g. from `.shp` or `.gpkg`):

```python
# Reproject and tag manually loaded GeoDataFrame
gdf = gdf.to_crs("EPSG:32633")
gdf = gdf[gdf.geom_type.isin(["LineString", "MultiLineString"])]
gdf["object_type"] = "road"
```

> ⚠️ The required column name is always `object_type` — this is how DamageScanner matches features with vulnerability curves and max damage values.

Save to file if needed:

```python
gdf.to_file("data/cleaned_osm_roads.gpkg", driver="GPKG")
```

---

## 💥 Example Usage with DamageScanner

```python
from damagescanner import DamageScanner

hazard = "data/flood_depth_100yr.tif"
osm_exposure = "data/cleaned_osm_roads.gpkg"
curves = "data/vulnerability_curves.csv"
maxdam = "data/maxdam.csv"

scanner = DamageScanner(hazard, osm_exposure, curves, maxdam)
damage = scanner.calculate()
```

This calculates direct damages using cleaned vector infrastructure data from OSM.

---

## 🔁 Tips for Reproducibility

- Document the OSM extract date and source (e.g. Geofabrik region)
- Use consistent object naming (e.g. `road`, `bridge`, `hospital`) to match your vulnerability inputs
- Save your processed data in `.gpkg` or `.shp` format for reuse

---

## 📚 See Also

- [Overview](./overview.md)
- [Raster-based approach](./raster.md)
- [Vector-based approach](./vector.md)
- [GlobalInfraRisk OSM Guide](https://vu-ivm.github.io/GlobalInfraRisk/howto/using_osm.html)