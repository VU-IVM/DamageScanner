# Overview

`DamageScanner()` is a Python toolkit for direct damage assessments of natural hazards. While originally designed for flood risk analysis, it can be used for **any hazard** where vulnerability can be expressed as a **one-dimensional curve** (e.g., flood depth, wind speed, ground shaking).

The tool is optimized for both **raster-based** (e.g. land use) and **vector-based** (e.g., roads, power plants, buildings) damage assessments. This page walks you through how it works, and what is required for running a successful analysis.

> 📚 For a very extensive overview of real-world examples, please refer to the [GlobalInfraRisk documentation](https://vu-ivm.github.io/GlobalInfraRisk/).

---

## Core Workflow

The DamageScanner logic consists of three key steps:

1. **Exposure Analysis** — identifies what assets intersect the hazard.
2. **Damage Calculation** — estimates damage using vulnerability curves.
3. **Risk Assessment** — aggregates damage across hazard return periods.

---

## Inputs Required

The `DamageScanner` class requires four key inputs:

### 1. **Hazard Data**
- Raster (GeoTIFF, NetCDF) or path to raster files
- Represents hazard intensity (e.g., flood depth, wind speed)

### 2. **Exposure Data** (aka `feature_data`)
- Vector formats: `.shp`, `.gpkg`, `.pbf`, `.geoparquet`, or GeoDataFrames
- Raster exposure layers (`.tif`, `.nc`) also supported
- Automatically detects type from file extension

### 3. **Vulnerability Curves**
- CSV or `pandas.DataFrame`
- Relates hazard intensity to damage (as fraction of max damage)

| Intensity | residential | industrial | farmland |
|-----------|-------------|------------|----------|
| 0.0       | 0.00        | 0.00       | 0.00     |
| 0.5       | 0.25        | 0.15       | 0.10     |
| 1.0       | 0.50        | 0.35       | 0.25     |

> ⚠️ **Important:** The unit of the first column/index must match the unit of the hazard layer (e.g., meters for flood depth).

### 4. **Maximum Damage Values**
- Specifies the max value per asset type (e.g., €/m² or €/asset)
- Provided as `dict`, CSV, or DataFrame

| landuse     | damage |
|-------------|--------|
| residential | 1000   |
| industrial  | 5000   |
| farmland    | 50     |

---

## Interactive Examples

To see `DamageScanner` in action with real-world data, explore our interactive Jupyter Notebook examples:

- 📊 **[Vector-based assessment](../examples/vector-based.ipynb)** — Demonstrates how to use OSM polygons and points (e.g., land use, buildings) with raster hazard maps.
- 🗺️ **[Raster-based assessment](../examples/raster-based.ipynb)** — Shows how to perform rapid large-scale analysis using land-use rasters.

---

## 📚 Next Steps

- 📦 [Raster-based approach](raster.md)
- 📦 [Vector-based approach](vector.md)
- 🧭 [Coupling with OSM](osm.md)


