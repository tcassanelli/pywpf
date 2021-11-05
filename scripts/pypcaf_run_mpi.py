#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
os.environ['MPLCONFIGDIR'] = "/project/v/vanderli/cassane/"
import numpy as np
import pypcaf


# PyPCAF execution

files = [
    'n0',
    # 'n25', 'n50', 'n75', 'n100', 'n200', 'n250',
    # 'n275', 'n300', 'n325'
    ]

files = ["data_B0531+21.npy", "data_B0540-69.npy", "data_B0833-45.npy"]

M_list = np.linspace(10, 200, 20, dtype=int)

# crab scinet
pypcaf.pca_folding(
    times_path=f'/scratch/v/vanderli/cassane/{files[0]}',
    dt=0.0001,
    T_init=0.03362167003,
    iteration=1000,         # 10 us scanning
    delta=1e-8,             # 10 ns step
    num_div=M,
    merit_func=pypcaf.merit1,
    region_order=3,
    work_dir='/scratch/v/vanderli/cassane/data_crab/'
    )

pypcaf.plot_all(
    pypcaf_path=(
        f'/scratch/v/vanderli/cassane/data_crab/pypcaf_out/{_f}-{_s}'
        ),
    T_ref=0.03362167003
    )