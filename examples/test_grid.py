import pandas as pd
import geopandas as gpd
import numpy as np
import sys
from pathlib import Path

# Add the src directory to Python path
repo_root = Path(__file__).parent.parent  # Adjust based on script location
sys.path.insert(0, str(repo_root / "src"))

# Now import
from damagescanner.vector import VectorScanner

data_path = Path(r"C:\Users\eks510\OneDrive - Vrije Universiteit Amsterdam\12_repositories\DamageScanner\data\kampen")

hazard = data_path / "hazard" / "1in100_inundation_map.tif"
exposure = data_path / "exposure" / "landuse.gpkg"
curves = data_path / "vulnerability" / "curves_landuse.csv"
maxdam = data_path / "vulnerability" / "maxdam_landuse.csv"

# Run both methods
result_gridded = VectorScanner(
    hazard_file=hazard,
    feature_file=exposure,
    curve_path=curves,
    maxdam_path=maxdam,
    object_col="landuse",
    gridded=True,
)

result_non_gridded = VectorScanner(
    hazard_file=hazard,
    feature_file=exposure,
    curve_path=curves,
    maxdam_path=maxdam,
    object_col="landuse",
    gridded=False,
)

print(f"Gridded features: {len(result_gridded)}")
print(f"Non-gridded features: {len(result_non_gridded)}")
print(f"\nGridded total damage: {result_gridded['damage'].sum():,.2f}")
print(f"Non-gridded total damage: {result_non_gridded['damage'].sum():,.2f}")
print(f"Difference: {abs(result_gridded['damage'].sum() - result_non_gridded['damage'].sum()):,.2f}")
print(f"Difference %: {abs(result_gridded['damage'].sum() - result_non_gridded['damage'].sum()) / result_non_gridded['damage'].sum() * 100:.1f}%")

# Check if same features are present
print(f"\n--- Feature count comparison ---")
print(f"Features only in gridded: {len(set(result_gridded.index) - set(result_non_gridded.index))}")
print(f"Features only in non-gridded: {len(set(result_non_gridded.index) - set(result_gridded.index))}")
print(f"Features in both: {len(set(result_gridded.index) & set(result_non_gridded.index))}")

# Compare coverage values for common features (coverage is a list)
common_idx = list(set(result_gridded.index) & set(result_non_gridded.index))
if common_idx:
    print(f"\n--- Coverage comparison for common features ---")
    
    # Sum coverage lists to get total coverage per feature
    gridded_coverage_sum = result_gridded.loc[common_idx, 'coverage'].apply(
        lambda x: np.sum(x) if hasattr(x, '__len__') else x
    )
    non_gridded_coverage_sum = result_non_gridded.loc[common_idx, 'coverage'].apply(
        lambda x: np.sum(x) if hasattr(x, '__len__') else x
    )
        
    # Count coverage items per feature
    gridded_coverage_len = result_gridded.loc[common_idx, 'coverage'].apply(
        lambda x: len(x) if isinstance(x, list) else 1
    )
    non_gridded_coverage_len = result_non_gridded.loc[common_idx, 'coverage'].apply(
        lambda x: len(x) if isinstance(x, list) else 1
    )
    
    print(f"Gridded total coverage (sum of all): {gridded_coverage_sum.sum():,.2f}")
    print(f"Non-gridded total coverage (sum of all): {non_gridded_coverage_sum.sum():,.2f}")
    print(f"Coverage sum difference: {abs(gridded_coverage_sum.sum() - non_gridded_coverage_sum.sum()):,.2f}")
    
    print(f"\nGridded total coverage items: {gridded_coverage_len.sum():,}")
    print(f"Non-gridded total coverage items: {non_gridded_coverage_len.sum():,}")

# Compare damage values for common features
if common_idx:
    print(f"\n--- Damage comparison for common features ---")
    
    gridded_damage = result_gridded.loc[common_idx, 'damage']
    non_gridded_damage = result_non_gridded.loc[common_idx, 'damage']
    
    damage_diff = (gridded_damage - non_gridded_damage).abs()
    
    print(f"Mean damage difference: {damage_diff.mean():,.2f}")
    print(f"Max damage difference: {damage_diff.max():,.2f}")
    print(f"Sum of damage differences: {damage_diff.sum():,.2f}")
    
    # Features where gridded > non-gridded
    gridded_higher = (gridded_damage > non_gridded_damage).sum()
    non_gridded_higher = (non_gridded_damage > gridded_damage).sum()
    equal = (gridded_damage == non_gridded_damage).sum()
    
    print(f"\nFeatures where gridded damage > non-gridded: {gridded_higher}")
    print(f"Features where non-gridded damage > gridded: {non_gridded_higher}")
    print(f"Features with equal damage: {equal}")

    # Find features with largest differences
    print(f"\n--- Top 10 features with largest absolute damage difference ---")
    top_diff_idx = damage_diff.nlargest(10).index
    
    comparison_df = pd.DataFrame({
        'gridded_damage': gridded_damage.loc[top_diff_idx],
        'non_gridded_damage': non_gridded_damage.loc[top_diff_idx],
        'difference': damage_diff.loc[top_diff_idx],
        'gridded_coverage_sum': gridded_coverage_sum.loc[top_diff_idx],
        'non_gridded_coverage_sum': non_gridded_coverage_sum.loc[top_diff_idx],
        'gridded_coverage_items': gridded_coverage_len.loc[top_diff_idx],
        'non_gridded_coverage_items': non_gridded_coverage_len.loc[top_diff_idx],
    })
    print(comparison_df.to_string())

    # Check if coverage difference explains damage difference
    print(f"\n--- Correlation analysis ---")
    coverage_diff = (gridded_coverage_sum - non_gridded_coverage_sum)
    correlation = coverage_diff.corr(damage_diff)
    print(f"Correlation between coverage diff and damage diff: {correlation:.4f}")