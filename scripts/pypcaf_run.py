#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import pypcaf

# PyPCAF execution

files = [
    'n0',
    # 'n25', 'n50', 'n75', 'n100', 'n200', 'n250',
    # 'n275', 'n300', 'n325'
    ]
series = ['000']
merit_func = [pypcaf.merit1, pypcaf.merit2, pypcaf.merit3]

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



# vela scinet
# for _s, _m in zip(series, merit_func):
#     for _f in files:
#         pypcaf.pcaf_single(
#             time_path=f'/scratch/v/vanderli/cassane/data_pulsar/{_f}.npy',
#             dt=0.002793,
#             T_init=0.089367,
#             iteration=1000,         # 10 us scanning
#             delta=1e-8,             # 10 ns step
#             num_div=range(3, 21),
#             merit_func=_m,
#             region_order=3,
#             work_dir='/scratch/v/vanderli/cassane/data_pulsar/'
#             )

#         pypcaf.plot_all(
#             pypcaf_path=(
#                 f'/scratch/v/vanderli/cassane/data_pulsar/pypcaf_out/{_f}-{_s}'
#                 ),
#             T_ref=0.08936715
#             )

# crab scinet
for _s, _m in zip(series, merit_func):
    for _f in files:
        pypcaf.pcaf_single(
            time_path=f'/scratch/v/vanderli/cassane/data_crab/{_f}.npy',
            dt=0.0001,
            T_init=0.033392,
            iteration=1000,         # 10 us scanning
            delta=1e-8,             # 10 ns step
            num_div=range(3, 21),
            merit_func=_m,
            region_order=3,
            work_dir='/scratch/v/vanderli/cassane/data_crab/'
            )

        pypcaf.plot_all(
            pypcaf_path=(
                f'/scratch/v/vanderli/cassane/data_crab/pypcaf_out/{_f}-{_s}'
                ),
            T_ref=0.03339241
            )


# vela
# for _s, _m in zip(series, merit_func):
#     for _f in files:
#         pypcaf.pcaf_single(
#             time_path=f'/Users/tomascassanelli/PCAF/data_vela/{_f}.npy',
#             dt=0.002793,
#             T_init=0.089367,
#             iteration=1000,         # 10 us scanning
#             delta=1e-8,             # 10 ns step
#             num_div=range(3, 21),
#             merit_func=_m,
#             region_order=3,
#             work_dir='/Users/tomascassanelli/PCAF/data_vela/'
#             )

#         pypcaf.plot_all(
#             pypcaf_path=(
#                 f'/Users/tomascassanelli/PCAF/data_vela/pypcaf_out/{_f}-{_s}'
#                 ),
#             T_ref=0.08936715
#             )

# import numpy as np
# pypcaf.folding(
#     time=np.load("/Users/tomascassanelli/PCAF/data_crab/n0.npy"),
#     dt=0.00001, T=0.03339241, num_div=5
#     )

# crab
# for _s, _m in zip(series, merit_func):
#     for _f in files:
#         pypcaf.pcaf_single(
#             time_path=f'/Users/tomascassanelli/PCAF/data_crab/{_f}.npy',
#             dt=0.0001,
#             T_init=0.03339241,
#             iteration=1000,
#             delta=1e-8,
#             num_div=range(3, 21),
#             merit_func=_m,
#             region_order=3,
#             work_dir='/Users/tomascassanelli/PCAF/data_crab/'
#             )

#         pypcaf.plot_all(
#             pypcaf_path=(
#                 f'/Users/tomascassanelli/PCAF/data_crab/pypcaf_out/{_f}-{_s}'
#                 ),
#             T_ref=0.0333924123
#             )
