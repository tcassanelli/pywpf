#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
os.environ['MPLCONFIGDIR'] = "/project/v/vanderli/cassane/"
import pypcaf

files = ["data_B0531+21.npy", "data_B0540-69.npy", "data_B0833-45.npy"]

# crab scinet
M = [2000]
base_dir = "/scratch/v/vanderli/cassane"
pypcaf.pca_folding(
    times_path=os.path.join(base_dir, "data_B0531+21.npy"),
    dt=0.0001,
    T_init=0.03362167003,
    iteration=1000,         # 10 us scanning
    delta=1e-8,             # 10 ns step
    num_div=M,
    merit_func=pypcaf.merit1,
    region_order=3,
    work_dir=base_dir
    )

pypcaf.plot_all(
    pypcaf_path=os.path.join(base_dir, "pypcaf_out/data_B0531+21-000"),
    T_ref=0.03362167003
    )


