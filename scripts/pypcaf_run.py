#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import pypcaf

# pypcaf.pcaf_double(
#     time_path='../../data_pulsar/n500.npy',
#     dt=0.002793,
#     T_init=0.089367,
#     iterations=[100, 400],
#     deltas=[1e-7, 1e-9],
#     num_div=range(6, 11),
#     merit_func=pypcaf.merit1,
#     region_order=3,
#     use_previous=False,
#     work_dir='../../data_pulsar/'
#     )

name = 'n0'

pypcaf.pcaf_single(
    time_path='../../data_pulsar/{}.npy'.format(name),
    dt=0.002793,
    T_init=0.089367,
    iteration=1000,         # 10 um scanning
    delta=1e-8,             # 10 ns step
    num_div=range(3, 23),
    merit_func=pypcaf.merit2,
    region_order=3,
    work_dir='../../data_pulsar/'
    )

pypcaf.plot_all_periods(
    pypcaf_path='../../data_pulsar/pypcaf_out/{}-001'.format(name),
    T_ref=0.08936715
    )
