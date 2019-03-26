#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import pypcaf

# PyPCAF execution

files = [
    'n0', 'n25', 'n50', 'n75', 'n100', 'n200', 'n250', 'n275', 'n300', 'n325'
    ]
series = ['000', '001']
merit_func = [pypcaf.merit2, pypcaf.merit3]

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


for f in files:

    # for m in merit_func:
    #     pypcaf.pcaf_single(
    #         time_path='../../data_pulsar2/{}.npy'.format(f),
    #         dt=0.002793,
    #         T_init=0.089367,
    #         iteration=1000,         # 10 us scanning
    #         delta=1e-8,             # 10 ns step
    #         num_div=range(3, 21),
    #         merit_func=m,
    #         region_order=3,
    #         work_dir='../../data_pulsar2/'
    #         )

    for s in series:
        pypcaf.plot_all(
            pypcaf_path='../../data_pulsar2/pypcaf_out/{}-{}'.format(f, s),
            T_ref=0.08936715
            )
