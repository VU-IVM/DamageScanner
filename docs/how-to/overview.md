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

> ⚠️ **Important:** The unit of the first column/index must match the unit of the hazard layer (e.g., meters for flood depth).

### 4. **Maximum Damage Values**
- Specifies the max value per asset type (e.g., €/m² or €/asset)
- Provided as `dict`, CSV, or DataFrame

---

## Example: Running DamageScanner

```python
from damagescanner import DamageScanner
import pandas as pd

# Paths to input files
hazard = "path/to/hazard_data.tif"
feature_data = "path/to/exposure_data.shp"
curves = "path/to/vulnerability_curves.csv"
maxdam = "path/to/maxdam.csv"

# Initialize
scanner = DamageScanner(hazard, feature_data, curves, maxdam)
```

---

## Step-by-Step Usage

### 1. `exposure()`
Identify which features overlap with the hazard.

```python
# Curves and maxdam must be empty DataFrames for exposure-only
scanner = DamageScanner(hazard, feature_data, pd.DataFrame(), pd.DataFrame())
exposed = scanner.exposure()
```

- Works with both vector and raster feature data
- Results in a GeoDataFrame with assets intersecting the hazard

---

### 2. `calculate()`
Applies vulnerability curves to estimate damage.

```python
results = scanner.calculate()
print(results.head())
```

- Requires valid curves and `maxdam`
- Damage is computed per asset, based on overlap and intensity

---

### 3. `risk()`
Calculates risk over multiple return periods.

```python
hazard_dict = {
    10: "hazard_10.tif",
    50: "hazard_50.tif",
    100: "hazard_100.tif"
}

risk_results = scanner.risk(hazard_dict)
```

- Computes **Expected Annual Damages (EAD)** using return periods
- Supports asset-specific curves and max damages

---

## Tips for Working with Geometry

> ⚠️ **Important:** Make sure each asset type uses one geometry type (Point or Polygon). Options:
>
> 1. Convert all assets to the same geometry
> 2. Split the damage calc by geometry type
> 3. Use custom object names like `substation_point` vs `substation_polygon`

---

## 📚 Next Steps

- 📦 [Raster-based approach](raster.md)
- 📦 [Vector-based approach](vector.md)
- 🧭 [Coupling with OSM](osm.md)


