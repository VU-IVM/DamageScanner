# Vector-based Approach

This page explains how to use `DamageScanner` with **vector-based exposure data** (e.g. shapefiles, GeoPackages, or OSM data). This approach is ideal for object-level damage estimation, such as buildings, roads, power plants, and other individual infrastructure assets.

---

## 🧠 When to Use Vector-Based Workflows

| Scenario | Raster-Based | Vector-Based |
|----------|--------------|---------------|
| Infrastructure objects (e.g. bridges, roads) | ❌ | ✅ |
| OpenStreetMap data | ❌ | ✅ |
| Detailed local studies | ❌ | ✅ |
| Exposure linked to specific geometries | ❌ | ✅ |

---

## 🔧 Required Inputs

As described in the [Overview](./overview.md), you need:

1. **Hazard raster** (e.g. flood depth, wind speed)
2. **Exposure vector** (e.g. buildings, roads — as `.shp`, `.gpkg`, `.pbf`, or `GeoDataFrame`)
3. **Vulnerability curves** (CSV or DataFrame)
4. **Maximum damage values** (CSV, dict, or DataFrame)

> ⚠️ Your exposure data must include a **column specifying asset type**, matching the keys used in the vulnerability curves and max damage data.

---

## 🧪 Minimal Working Example

```python
from damagescanner import DamageScanner

hazard = "data/hazard_flood_depth.tif"
feature_data = "data/infrastructure_osm.gpkg"
curves = "data/vulnerability_curves.csv"
maxdam = "data/maxdam.csv"

scanner = DamageScanner(hazard, feature_data, curves, maxdam)
damage = scanner.calculate()
```

The result is a `GeoDataFrame` of features with estimated direct damage values per asset.

---

## 🔍 Key Behavior

- The exposure file must contain geometry columns (Point, LineString, or Polygon)
- `DamageScanner` overlays each geometry with the hazard raster
- It then samples the hazard value and looks up the damage fraction from the corresponding vulnerability curve
- The damage is calculated as:

  \
  `exposure_area * damage_fraction * max_damage`

---

## ⚠️ Geometry Handling Tips

> ⚠️ **Mixed Geometry Types** — Avoid mixing Points and Polygons for the same asset type.
>
> ⚠️ **CRS Alignment** — Both hazard and vector data must be in the same coordinate reference system.
>
> ⚠️ **Object Column** — You must define a column that links each object to a vulnerability curve and max damage.

---

## 📚 See Also

- [Overview](./overview.md)
- [Raster-based approach](./raster.md)
- [Coupling with OSM](./osm.md)
