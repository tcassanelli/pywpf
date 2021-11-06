#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
os.environ['MPLCONFIGDIR'] = "/project/v/vanderli/cassane/"
import pypcaf

M = [200]
base_dir = "/scratch/v/vanderli/cassane"
pypcaf.pca_folding(
    times_path=os.path.join(base_dir, "data_B0531+21.npy"),
    # dt=0.0001,
    dt=2.5e-7,
    T_init=0.0336372543236884,
    iteration=1000,         # 10 us scanning
    delta=1e-8,             # 10 ns step
    num_div=M,
    merit_func=pypcaf.merit1,
    region_order=3,
    work_dir=base_dir
    )
