#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
from setuptools import setup
# from Cython.Build import cythonize

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
        'cython'
    ],
    # ext_modules=cythonize("pypcaf/pcaf_functions.pyx"),
    packages=['pypcaf'],
    package_dir={
        'pypcaf': 'pypcaf',
        # 'pyoof.submodule': 'pyoof/submodel',
        },
    package_data={
        '': ['pypcaf_sty.mplstyle'],
        # '': ['config_params.yaml']
        },
    long_description='''
    PyPCAF - Principal Componen Analysis (PCA) Folding.
    '''
    )
