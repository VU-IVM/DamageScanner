from setuptools import setup

def long_description():
    with open('README.md', encoding='utf-8') as f:
        return f.read()

setup(
    name='damagescanner',    
    version='0.6.2',
    long_description=long_description(),
    long_description_content_type='text/markdown',    
    install_requires=[
    "geopandas",
    "matplotlib",	
    "numpy",
    "pyproj",
    "rasterio",
    "xarray",
    "rioxarray",
    "packaging",
    "pandas",
    "shapely",
    ],
    project_urls={
        'GitHub': 'https://github.com/ElcoK/DamageScanner'
    },

)