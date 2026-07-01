# Vector-based Approach

This page explains how to use `DamageScanner` with **vector-based exposure data** (e.g. shapefiles, GeoPackages, or OSM data). This approach is ideal for object-level damage estimation, such as buildings, roads, power plants, and other individual infrastructure assets.

---

## When to Use Vector-Based Workflows

| Scenario | Raster-Based | Vector-Based |
|----------|--------------|---------------|
| Infrastructure objects (e.g. bridges, roads) | ❌ | ✅ |
| OpenStreetMap data | ❌ | ✅ |
| Detailed local studies | ❌ | ✅ |
| Exposure linked to specific geometries | ❌ | ✅ |

---

## Required Inputs

As described in the [Overview](./overview.md), you need:

1. **Hazard raster** (e.g. flood depth, wind speed)
2. **Exposure vector** (e.g. buildings, roads — as `.shp`, `.gpkg`, `.pbf`, or `GeoDataFrame`)
3. **Vulnerability curves** (CSV or DataFrame)
4. **Maximum damage values** (CSV, dict, or DataFrame)

> For a detailed walkthrough with real data and visualizations, see our **[Vector-based example notebook](../examples/vector-based.ipynb)**.

---

## Key Behavior

- The exposure file must contain geometry columns (Point, LineString, or Polygon)
- `DamageScanner` overlays each geometry with the hazard raster
- It then samples the hazard value and looks up the damage fraction from the corresponding vulnerability curve
- The damage is calculated as:

  `exposure_area * damage_fraction * max_damage`

---

## 📚 See Also

- [Overview](./overview.md)
- [Raster-based approach](./raster.md)
- [Coupling with OSM](./osm.md)

## API Documentation

::: damagescanner.vector