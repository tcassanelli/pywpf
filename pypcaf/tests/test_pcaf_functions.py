#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import pytest
import os
import numpy as np
from astropy.utils.misc import NumpyRNGContext
from numpy.testing import assert_allclose, assert_equal
import pypcaf

# Plot style added from relative path
current_dir = os.path.dirname(__file__)
n0_path = os.path.join(current_dir, 'data', 'n0.npy')
folding_path = os.path.join(current_dir, 'data', 'folding_true.npz')


def test_nextpow2():

    with NumpyRNGContext(0):
        n = np.random.randint(1, 100)

    m = pypcaf.nextpow2(n)
    m_true = 6

    assert_allclose(m, m_true)


def test_flat_region_finder():

    with NumpyRNGContext(0):
        X = np.random.normal(0, 30, 1000)

    idx_max = pypcaf.flat_region_finder(X)
    idx_max_true = 108
    assert_allclose(idx_max, idx_max_true)


# Testing new function
def test_folding():

    n0 = np.load(n0_path)
    folding_true = np.load(folding_path)

    remainder, waterfall = pypcaf.folding(
        time=n0, dt=0.002793, T=0.089367, num_div=20
        )

    assert_allclose(remainder, folding_true['remainder'])
    assert_allclose(waterfall, folding_true['waterfall'])



