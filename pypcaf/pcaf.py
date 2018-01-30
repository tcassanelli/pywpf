#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
import time
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from .pcaf_functions import find_period

__all__ = ['pcaf']


def pcaf(
    path_time, T_init, dt, iteration1, iteration2, delta1, delta2, num_div,
    merit_func
        ):

    start_time = time.time()

    print('\n ******* PyPCAF: finding pulsar period ******* \n')
    print('... Total number of num_div loops: {} ... \n'.format(len(num_div)))
    print('... Reading data ... \n')
    name = os.path.split(path_time)[1][:-4]  # input must be .npy file
    print('... Extracting period from {}.npy ... \n'.format(name))

    # calling the time array
    time_pulsar_data = np.load(path_time)

    EVALW1, EVALW2 = [], []
    SW1, SW2 = [], []
    MERIT1, MERIT2 = [], []

    pcaf_info = Table(
        names=[
            'num_div', 'dt', 'iter1', 'iter2', 'delta1', 'delta2',
            'T_init', 'T_est1', 'T_est2', 'idx1_max', 'idx2_max'
            ]
        )

    # M is the number of divisions
    i_loops = 1  # loop/iteration counter
    for M in num_div:

        print('... Computing loop {} ...\n'.format(i_loops))
        i_loops += 1

        (
            T, [idx1_max, idx2_max], [EValw1, EValw2], [Sw1, Sw2], [M1, M2]
            ) = find_period(
            time=time_pulsar_data,
            T_init=T_init,
            dt=dt,
            num_div=M,
            iteration1=iteration1,
            delta1=delta1,
            iteration2=iteration2,
            delta2=delta2,
            merit_func=merit_func
            )

        [T_init, T_est1, T_est2] = T

        EVALW1.append(EValw1)
        EVALW2.append(EValw2)

        MERIT1.append(M1)
        MERIT2.append(M2)

        SW1.append(Sw1)
        SW2.append(Sw2)

        pcaf_info.add_row([
            M, dt, iteration1, iteration2, delta1, delta2,
            T_init, T_est1, T_est2, idx1_max, idx2_max
            ])

    print(pcaf_info)
    print('\n... Storing data ... \n')

    # Making sub-directory to store data
    for j in ["%03d" % i for i in range(101)]:
        dir_name = os.path.join('pcaf_out', name + '-' + str(j))
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            break

    ascii.write(
        output=os.path.join(dir_name, 'pcaf_info.dat'),
        table=pcaf_info
        )

    for M, idx in zip(num_div, range(len(num_div))):
        np.savez(
            os.path.join(dir_name, 'pcaf_out_M{}'.format(M)),
            EVALW1=EVALW1[idx],
            EVALW2=EVALW2[idx],
            MERIT1=MERIT1[idx],
            MERIT2=MERIT2[idx],
            SW1=SW1[idx],
            SW2=SW2[idx]
            )

    final_time = np.round((time.time() - start_time) / 60, 2)
    print(
        '\n **** PyPCAF Completed at {} mins **** \n'.format(final_time)
        )
