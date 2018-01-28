import os
import time
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from pca_analysis import find_period
from merit_functions import merit1, merit2

start_time = time.time()

print('\n ******* PCA Folding: finding pulsar period ******* \n')

path_to_file = '../../data_pulsar'

# Initial values to start the execution
dt = 0.002793            # 4 ms, 0.002793 s
T_init = 0.089367  # Initial period, usualy well known

iter1 = 100
delta1 = 1e-7

iter2 = 400  # 1000, 400
delta2 = 1e-8

num_div = range(1, 9)

print('... total num_div: {} ... \n'.format(len(num_div)))
print('... Reading data ... \n')

# calling the time array
time_pulsar_data = np.load(os.path.join(path_to_file, 'n0.npy'))

EVALW1, EVALW2 = [], []
SW1, SW2 = [], []
MERIT1, MERIT2 = [], []

pcaf_info = Table(
    names=[
    'num_div', 'dt', 'iter1', 'iter2', 'delta1', 'delta2', 'T_init', 'T_est1',
    'T_est2', 'idx1_max', 'idx2_max'
        ]
    )

# M is the number of divisions
for M in num_div:

    (
        T, [idx1_max, idx2_max], [EValw1, EValw2], [Sw1, Sw2], [M1, M2]
        ) = find_period(
        time=time_pulsar_data,
        T_init=T_init,
        dt=dt,
        num_div=M,
        iter1=iter1,
        delta1=delta1,
        iter2=iter2,
        delta2=delta2,
        merit_func=merit2
        )

    [T_init, T_est1, T_est2] = T

    EVALW1.append(EValw1)
    EVALW2.append(EValw2)

    MERIT1.append(M1)
    MERIT2.append(M2)

    SW1.append(Sw1)
    SW2.append(Sw2)

    pcaf_info.add_row([
        M, dt, iter1, iter2, delta1, delta2, T_init, T_est1, T_est2,
        idx1_max, idx2_max
        ])

    print('... {} iteration(s) completed ...\n'.format(M))

print(pcaf_info)
print('\n... Storing data ... \n')

ascii.write(output='pcaf_info.dat', table=pcaf_info)

for M in range(len(num_div)):
    np.savez(
        'pcaf_out_M{}'.format(M),
        EVALW1=EVALW1[M],
        EVALW2=EVALW2[M],
        MERIT1=MERIT1[M],
        MERIT2=MERIT2[M],
        SW1=SW1[M],
        SW2=SW2[M]
        )

final_time = np.round((time.time() - start_time) / 60, 2)
print('\n **** PCA Folding Completed at {} mins **** \n'.format(final_time))
