#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import pypcaf

# pypcaf.pcaf_double(
#     time_path='../../data_pulsar/n500.npy',
#     dt=0.002793,
#     T_init=0.089367,
#     iterations=[100, 400],
#     deltas=[1e-7, 1e-8],
#     num_div=range(6, 11),
#     merit_func=pypcaf.merit1,
#     region_order=3,
#     use_previous=False,
#     work_dir=None
#     )

# pypcaf.pcaf_single(
#     time_path='../../data_pulsar/n3000.npy',
#     dt=0.002793,
#     T_init=0.089367,
#     iteration=1000,
#     delta=1e-8,
#     num_div=range(30, 31),
#     merit_func=pypcaf.merit1,
#     region_order=3,
#     work_dir=None
#     )

pypcaf.plot_all_periods(
    pypcaf_path='pypcaf_out/n3000-001',
    T_ref=0.08936715
    )
