# Raster-based Approach

This page explains how to use `DamageScanner` with **raster-based exposure and hazard data**. This approach is especially useful when working with gridded exposure datasets, such as population grids, land use rasters, or economic value rasters.

---

## When to Use Raster-Based Workflows

| Scenario | Raster-Based | Vector-Based |
|----------|--------------|---------------|
| Gridded exposure data (e.g. land cover) | ✅ | ❌ |
| Building footprints or linear infrastructure | ❌ | ✅ |
| Country- or global-scale risk analysis | ✅ | ✅ |
| Exposure not tied to individual assets | ✅ | ❌ |

---

## Required Inputs

You still need the same four key inputs as explained in the [Overview](./overview.md):

1. **Hazard raster** (e.g. flood depth, wind speed)
2. **Exposure raster** (e.g. land use, population density)
3. **Vulnerability curves** (CSV or DataFrame)
4. **Maximum damage values** (CSV, dict, or DataFrame)

> For a detailed walkthrough with real data and visualizations, see our **[Raster-based example notebook](../examples/raster-based.ipynb)**.

---

## Notes on Raster Behavior

- Exposure values are interpreted as **damageable units per cell** (e.g., people, € value, or m²)
- Output damage is computed as:

  \
  `exposure * damage_fraction * max_damage`

- The damage fraction is determined from the vulnerability curve, based on the hazard intensity in the same raster cell.

---

## Common Pitfalls

> ⚠️ **Misaligned Rasters** — Make sure both rasters have the same extent, resolution, and CRS. You can use `rasterio.warp` or `gdalwarp` to resample.

> ⚠️ **Missing Data** — NoData values can propagate through your analysis. Clean them or mask them before calculation.

> ⚠️ **Unit mismatch** — Hazard units and curve x-axis must match exactly (e.g., meters, m/s).

---

## 📚 See Also

- [Overview](./overview.md)
- [Vector-based approach](./vector.md)
- [Coupling with OSM](./osm.md)

## API Documentation

::: damagescanner.raster