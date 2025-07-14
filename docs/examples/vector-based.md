# Vector-based Example (Kampen)

This example demonstrates how to run the vector-based flood damage assessment using OSM exposure data for the Kampen region in the Netherlands.

## 1. Import Dependencies

```python
import pandas as pd
import damagescanner
from damagescanner.core import DamageScanner
```

---

## 2. Load Input Data from GitHub

We use publicly available data from the [DamageScanner GitHub repository](https://github.com/VU-IVM/DamageScanner/tree/installation/data/kampen):

```python
hazard_file = "https://raw.githubusercontent.com/VU-IVM/DamageScanner/installation/data/kampen/hazard/1in1000_inundation_map.tif"
exposure_file = "https://raw.githubusercontent.com/VU-IVM/DamageScanner/installation/data/kampen/exposure/kampen.osm.pbf"
curve_file = "https://raw.githubusercontent.com/VU-IVM/DamageScanner/installation/data/kampen/vulnerability/curves_osm.csv"
maxdam_file = "https://raw.githubusercontent.com/VU-IVM/DamageScanner/installation/data/kampen/vulnerability/maxdam_osm.csv"
```

---

## 3. Load Vulnerability Data

```python
curves = pd.read_csv(curve_file)
maxdam = pd.read_csv(maxdam_file)
```

---

## 4. Run the Vector-Based Damage Assessment

```python
scanner = DamageScanner(hazard_file, exposure_file, curves, maxdam)
features = scanner.calculate(asset_type='main_roads')
```

---

## 5. Display Output

```python
features[["object_type", "geometry", "damage"]].head()
```

---

