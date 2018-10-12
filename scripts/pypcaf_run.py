#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import pypcaf

name = ['n85', 'n90', 'n95', 'n100']

# pypcaf.pcaf_double(
#     time_path='../../data_pulsar2/{}.npy'.format(name),
#     dt=0.002793,
#     T_init=0.089367,
#     iterations=[100, 400],
#     deltas=[1e-7, 1e-9],
#     num_div=range(3, 10),
#     merit_func=pypcaf.merit2,
#     region_order=3,
#     use_previous=False,
#     work_dir='../../data_pulsar2/'
#     )

for n in name:

    pypcaf.pcaf_single(
        time_path='../../data_pulsar2/{}.npy'.format(n),
        dt=0.002793,
        T_init=0.089367,
        iteration=1000,         # 10 us scanning
        delta=1e-8,             # 10 ns step
        num_div=range(3, 21),
        merit_func=pypcaf.merit2,
        region_order=3,
        work_dir='../../data_pulsar2/'
        )

    pypcaf.plot_all_periods(
        pypcaf_path='../../data_pulsar2/pypcaf_out/{}-000'.format(n),
        T_ref=0.08936715
        )
