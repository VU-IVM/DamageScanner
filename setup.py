#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup DamageScanner package
"""
from glob import glob
from os.path import basename, splitext

from setuptools import find_packages
from setuptools import setup


def readme():
    """Read README contents
    """
    with open('README.md') as f:
        return f.read()


setup(
    name='damagescanner',
    version='0.4.0',
    license='MIT License',
    description='Damage assessment tool for natural disasters',
    long_description=readme(),
    long_description_content_type="text/markdown",
    author='Elco Koks',
    author_email='elcokoks@gmail.com',
    url='https://github.com/ElcoK/DamageScanner',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 1 - Planning',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: GIS',
        'Topic :: Utilities',
    ],
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    install_requires=[
        # eg: 'aspectlib==1.1.1', 'six>=1.7',
        'shapely>=1.6',
        'geopandas>=0.4.0',
	'pandas>=0.23.4',
	'rasterio>=1.0.8',
	'numpy>=1.15.2',
	'matplotlib>=3.0.0',
    'scipy>=1.0.0',
    'tqdm>=4.30.0'
	
    ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    entry_points={
        'console_scripts': [
            # eg: 'snkit = snkit.cli:main',
        ]
    },
)
