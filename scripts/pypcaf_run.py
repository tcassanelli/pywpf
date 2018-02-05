#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import pypcaf

time_path = '../../data_pulsar/n500.npy'

# Initial values to start the execution
dt = 0.002793            # 4 ms, 0.002793 s
T_init = 0.089367        # Initial period, usualy well known
T_ref = 0.08936715

iteration1 = 100
delta1 = 1e-7

iteration2 = 400
delta2 = 1e-8

num_div = range(3, 21)  # 2 or more, not 0, or 1


pypcaf.pcaf(
    time_path=time_path,
    dt=dt,
    T_init=T_init,
    iterations=[iteration1, iteration2],
    deltas=[delta1, delta2],
    num_div=num_div,
    merit_func=pypcaf.merit1,
    region_finder=False,
    use_previous=False
    )

pypcaf.plot_all_periods(
    pypcaf_path='pypcaf_out/n500-000',
    T_ref=T_ref
    )
