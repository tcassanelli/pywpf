#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
from setuptools import setup

setup(
    name='pypcaf',
    version='0.1',
    description='pypcaf - PCA-Folding',
    author='Tomas Cassanelli',
    author_email='tcassanelli@gmail.com',
    install_requires=[
        'setuptools',
        'numpy>=1.8',
        'scipy>=0.15',
        'astropy>=1.1',
        'matplotlib>=2.0',
        'fast_histogram',
        # 'cython>=0.26.1'
        ],
    setup_requires=['pytest-runner'],  # Enables tests
    tests_require=['pytest'],
    # ext_modules=extensions,
    packages=['pypcaf'],
    package_dir={
        'pypcaf': 'pypcaf',
        # 'pyoof.submodule': 'pypcaf/submodel',
        },
    include_package_data=True,
    long_description='''
    PyPCAF - Principal Componen Analysis (PCA) Folding.
    '''
    )
