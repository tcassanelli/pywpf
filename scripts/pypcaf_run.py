#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
os.environ['MPLCONFIGDIR'] = "/project/v/vanderli/cassane/"
import pypcaf

# B0531+21
# M = [200]
# base_dir = "/scratch/v/vanderli/cassane"
# pypcaf.pca_folding(
#     times_path=os.path.join(base_dir, "data_B0531+21.npy"),
#     dt=0.0001,
#     T_init=0.0336372543236884,
#     iteration=1000,
#     delta=5e-8,
#     num_div=M,
#     merit_func=pypcaf.merit1,
#     region_order=3,
#     work_dir=base_dir
#     )

# B0540-69
M = [20]
base_dir = "/scratch/v/vanderli/cassane"
pypcaf.pca_folding(
    times_path=os.path.join(base_dir, "data_B0540-69.npy"),
    dt=0.0001,
    T_init=0.05064997229,
    iteration=1000,
    delta=0.1e-9,
    num_div=M,
    merit_func=pypcaf.merit1,
    region_order=3,
    work_dir=base_dir
    )

# B0833-45
M = [20]
base_dir = "/scratch/v/vanderli/cassane"
pypcaf.pca_folding(
    times_path=os.path.join(base_dir, "data_B0833-45.npy"),
    dt=0.0001,
    T_init=0.08936644835,
    iteration=1000,
    delta=1e-9,
    num_div=M,
    merit_func=pypcaf.merit1,
    region_order=3,
    work_dir=base_dir
    )