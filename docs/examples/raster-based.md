# Raster-based Example (Kampen)

This example demonstrates how to run the raster-based flood damage assessment using a land-use map for the Kampen region in the Netherlands.

## 1. Import Dependencies

```python
import pandas as pd
from damagescanner.raster import RasterScanner
```

---

## 2. Load Input Data from GitHub

We use publicly available data from the [DamageScanner GitHub repository](https://github.com/VU-IVM/DamageScanner/tree/installation/data/kampen):

```python
hazard_file = "https://raw.githubusercontent.com/VU-IVM/DamageScanner/installation/data/kampen/hazard/1in1000_inundation_map.tif"
landuse_file = "https://raw.githubusercontent.com/VU-IVM/DamageScanner/installation/data/kampen/exposure/landuse_map.tif"
curve_file = "https://raw.githubusercontent.com/VU-IVM/DamageScanner/installation/data/kampen/vulnerability/curves.csv"
maxdam_file = "https://raw.githubusercontent.com/VU-IVM/DamageScanner/installation/data/kampen/vulnerability/maxdam.csv"
```

---

## 3. Load Vulnerability Data

```python
curves = pd.read_csv(curve_file)
maxdam = pd.read_csv(maxdam_file)
```

---

## 4. Run the Raster-Based Damage Assessment

```python
damage_df, damagemap, landuse_data, hazard_data = RasterScanner(
    landuse_file,
    hazard_file,
    curves,
    maxdam,
)
```

---

## 5. Display Output

```python
damage_df.head()
```

