#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
from setuptools import setup

setup(
    name='pywpf',
    version='0.1',
    description='PyWPF - Waterfall-PCA folding',
    author='Tomas Cassanelli',
    author_email='tcassanelli@protonmail.com',
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
    packages=['pywpf'],
    package_dir={
        'pywpf': 'pywpf',
        # 'pyoof.submodule': 'pypcaf/submodel',
        },
    include_package_data=True,
    long_description='''
    PyWPF - Waterfall-PCA folding
    '''
    )
