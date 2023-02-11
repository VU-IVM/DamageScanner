#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup DamageScanner package
"""
from glob import glob
from os.path import basename, splitext

from setuptools import find_packages
from setuptools import setup


def readme():
    """Read README contents"""
    with open("README.md", encoding="utf8") as f:
        return f.read()

setup(
    name="damagescanner",
    version="0.6.0",
    license="MIT License",
    description='Damage assessment tool for natural disasters',
    long_description=readme(),
    long_description_content_type="text/markdown",
    author='Elco Koks',
    author_email='elcokoks@gmail.com',
    url='https://github.com/ElcoK/DamageScanner',
    package_dir={"": "src"},  # Optional
    packages=find_packages(where="src"),
    python_requires=">=3.4, <4",
    include_package_data=True,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Utilities",
    ],

    install_requires=[
        "pandas",
        "geopandas",
        "xarray",
        "rasterio",
        "rioxarray",
        "numpy",
        "matplotlib",
        "tqdm"
    ]
)
