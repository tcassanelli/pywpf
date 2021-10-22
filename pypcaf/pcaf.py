#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
import time as TIME
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from .pcaf_functions import find_period, CP_value

__all__ = ['pca_folding']


def pca_folding(
    times_path, dt, T_init, iteration, delta, num_div, merit_func,
    region_order, work_dir=None
        ):
    """
    Core function for the `~pypcaf` package. It computes one estimated period
    given an inital period and a ``times`` array (one dimenional), with its
    corresponent time bin.

    Parameters
    ----------
    times_path : `str`
        String with the path to the file. The file must be a text-like file,
        e.g. .csv, .txt, ext. It can also be a python binary file .npy. In any
        case it must be always a one dimensional array.
    dt : `float`
        time bin.
    T_init : `float`
        Initial period to start looking for a best estimate, ``T_est``.
    iterations : `list`
        Corresponds to ``[iteration1, iteration2]``.
    delta : `float`
        Increase of the period in each iteration. The recommended order of it
        is between ``1e-7`` and ``1e-8``.
    num_div : `int`
        Number of divisions made to the ``time`` array, which corresponds to the
        number of elements in a row of the waterfall matrix.
    merit_func : `function`
        It computes the merit function from eigenvalues and scalar arrays.
        Both of them should be a one dimensional array.
    region_order : `int`
        It makes use of the `~pypcaf.flat_region_finder` to search for the
        maximum in the selected merit function. If ``region_order = 1``,
        it will compute the ordinary maximum of the merit array, i.e.
        ``np.max(merit)``. This function defines the estimated period in both
        ``iterations``.
    work_dir : `str`
        Default is `None`, it will store the ``pypcaf_out/`` folder in
        elsewhere. The current configuration stores files next to the
        `~pypcaf.pcaf` script.
    """

    start_time = TIME.time()

    print('\n ******* PyPCAF: finding pulsar period (single) ******* \n')

    if not all(num_div[i] <= num_div[i + 1] for i in range(len(num_div) - 1)):
        raise TypeError('num_div has to be sorted')

    if not num_div[0] >= 3:
        raise TypeError('num_div has to be equal greater than 3')

    num_div = np.array(num_div)
    print('... Total number of num_div loops: {} ... \n'.format(num_div.size))

    base = os.path.basename(times_path)
    name = os.path.splitext(base)[0]

    print('... Extracting period from {} ... \n'.format(base))

    # calling the times array
    if os.path.splitext(base)[1] == '.npy':
        times = np.load(times_path)
    else:
        times = np.genfromtxt(times_path)

    EVALW, SW, MERIT = [], [], []

    pypcaf_info = Table(
        names=[
            'num_div', 'dt', 'iter', 'delta', 'T_init', 'T_est', 'idx_max',
            'region_order', 'MAX', 'STD', 'MEAN', 'CP'
            ],
        dtype=['int32', 'float128', 'int32'] + ['float128'] * 3 +
        ['int32'] * 2 + ['float128'] * 4
        )

    # M is the number of divisions
    i_loops = 1  # loop/iteration counter
    for i in range(num_div.size):

        # print('... Computing loop {} ...\n'.format(i_loops))
        i_loops += 1

        T_est, EValw, Sw, merit, idx_max = find_period(
            times=times,
            dt=dt,
            T_init=T_init,
            num_div=num_div[i],
            iteration=iteration,
            delta=delta,
            merit_func=merit_func,
            region_order=region_order
            )

        EVALW.append(EValw)
        MERIT.append(merit)
        SW.append(Sw)

        pypcaf_info.add_row([
            num_div[i], dt, iteration, delta, T_init, T_est, idx_max,
            region_order, merit.max()
            ] + [CP_value(merit=merit, idx_max=idx_max)])

        # Printing in every iteration another row
        if i == 0:
            pypcaf_info.pprint(max_width=-1)
        else:
            print(pypcaf_info.pformat(max_width=-1)[i + 2])

    # Printing summary
    # pypcaf_info.pprint(max_lines=-1, max_width=-1)
    print('\n... Storing data ... \n')

    if work_dir is None:
        work_dir = ''

    # Making sub-directory to store data
    for j in ["%03d" % i for i in range(101)]:
        dir_name = os.path.join(work_dir, 'pypcaf_out', name + '-' + str(j))
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
            break

    ascii.write(
        output=os.path.join(dir_name, 'info.dat'),
        table=pypcaf_info
        )

    for M, idx in zip(num_div, range(num_div.size)):
        np.savez(
            os.path.join(dir_name, 'M{}'.format(M)),
            EVALW=EVALW[idx],
            MERIT=MERIT[idx],
            SW=SW[idx]
            )

    final_time = np.round((TIME.time() - start_time) / 60, 2)
    print(
        '\n **** PyPCAF Completed at {} mins **** \n'.format(final_time)
        )
