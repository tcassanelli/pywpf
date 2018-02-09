#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
import time as TIME
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from .pcaf_functions import find_period, find_period2

__all__ = ['pcaf_single', 'pcaf_double']


def pcaf_single(
    time_path, dt, T_init, iteration, delta, num_div, merit_func,
    region_order, work_dir=None
        ):

    start_time = TIME.time()

    print('\n ******* PyPCAF: finding pulsar period (single) ******* \n')

    if not all(num_div[i] <= num_div[i + 1] for i in range(len(num_div) - 1)):
        raise TypeError('num_div has to be sorted')

    if not num_div[0] >= 3:
        raise TypeError('num_div has to be equal greater than 3')

    num_div = np.array(num_div)
    print('... Total number of num_div loops: {} ... \n'.format(num_div.size))

    base = os.path.basename(time_path)
    name = os.path.splitext(base)[0]

    print('... Extracting period from {} ... \n'.format(base))

    # calling the time array
    if os.path.splitext(base)[1] == '.npy':
        time = np.load(time_path)
    else:
        time = np.genfromtxt(time_path)

    EVALW, SW, MERIT = [], [], []

    pypcaf_info = Table(
        names=[
            'num_div', 'dt', 'iter', 'delta', 'T_init', 'T_est', 'idx_max',
            'region_order'
            ],
        dtype=['int32', 'float64', 'int32'] + ['float64'] * 3 + ['int32'] * 2
        )

    # M is the number of divisions
    i_loops = 1  # loop/iteration counter
    for i in range(num_div.size):

        print('... Computing loop {} ...\n'.format(i_loops))
        i_loops += 1

        T_est, EValw, Sw, merit, idx_max = find_period(
            time=time,
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
            num_div[i], dt, iteration, delta, T_init, T_est,
            idx_max, region_order
            ])

    # Printing summary
    pypcaf_info.pprint(max_lines=-1, max_width=-1)
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


def pcaf_double(
    time_path, dt, T_init, iterations, deltas, num_div, merit_func,
    region_order, use_previous, work_dir=None
        ):
    """
    Core functionfor the `~pypcaf` package. It computes two estimated periods
    given an inital period and a time array (one dimenional), with its
    corresponent time bin.

    Parameters
    ----------
    time_path : `str`
        String with the path to the file. The file must be a text-like file,
        e.g. .csv, .txt, ext. It can also be a python binary file .npy. In any
        case it must be always a one dimensional array.
    dt : `float`
        time bin.
    T_init : `float`
        Initial period to start looking for a best estimate, ``T_est``.
    iterations : `list`
        Corresponds to ``[iteration1, iteration2]``.
    deltas : `float`
        Corresponds to ``[delta1, delta2]``
    num_div : `list` or `~numpy.ndarray` or `range`
        Number of divisions  made to the time array, which corresponds to the
        number of elements in a row of the waterfall matrix. The number of
        divisions has to be equal greater than 3, and increasing sorted.
    merit_func : `function`
        It computes the merit function from eigenvalues and scalar arrays.
        Both of them should be a one dimensional array.
    region_order : `int`
        It makes use of the `~pypcaf.flat_region_finder` to search for the
        maximum in the selected merit function. If ``region_order = 1``,
        it will compute the ordinary maximum of the merit array, i.e.
        ``np.max(merit)``. This function defines the estimated period in both
        ``iterations``.
    use_previous : `bool`
        If `True` will run the first loop of the ``num_div`` `list`, and use
        its results, i.e. ``T_est2``, from the next loop in ``num_div``
        initial period, ``T_init``. If set to `False`, ``T_init`` will not
        change.
    work_dir : `str`
        Default is `None`, it will store the ``pypcaf_out/`` folder in
        elsewhere. The current configuration stores files next to the
        `~pypcaf.pcaf` script.
    """

    start_time = TIME.time()

    print('\n ******* PyPCAF: finding pulsar period (double) ******* \n')

    if not all(num_div[i] <= num_div[i + 1] for i in range(len(num_div) - 1)):
        raise TypeError('num_div has to be sorted')

    if not num_div[0] >= 3:
        raise TypeError('num_div has to be equal greater than 3')

    num_div = np.array(num_div)
    print('... Total number of num_div loops: {} ... \n'.format(num_div.size))

    base = os.path.basename(time_path)
    name = os.path.splitext(base)[0]

    print('... Extracting period from {} ... \n'.format(base))

    # calling the time array
    if os.path.splitext(base)[1] == '.npy':
        time = np.load(time_path)
    else:
        time = np.genfromtxt(time_path)

    EVALW1, EVALW2 = [], []
    SW1, SW2 = [], []
    MERIT1, MERIT2 = [], []

    pypcaf_info = Table(
        names=[
            'num_div', 'dt', 'iter1', 'iter2', 'delta1', 'delta2',
            'T_init', 'T_est1', 'T_est2', 'idx1_max', 'idx2_max',
            'region_order'
            ],
        dtype=['int32', 'float64'] + ['int32'] * 2 + ['float64'] * 5 +
        ['int32'] * 3
        )

    # M is the number of divisions
    i_loops = 1  # loop/iteration counter
    for i in range(num_div.size):

        print('... Computing loop {} ...\n'.format(i_loops))
        i_loops += 1

        if use_previous and i != 0:
            T_init = pypcaf_info['T_est2'][i - 1]
        else:
            pass

        T, EValw, Sw, merit, idx_max = find_period2(
            time=time,
            dt=dt,
            T_init=T_init,
            num_div=num_div[i],
            iterations=iterations,
            deltas=deltas,
            merit_func=merit_func,
            region_order=region_order
            )

        EVALW1.append(EValw[0])
        EVALW2.append(EValw[1])

        MERIT1.append(merit[0])
        MERIT2.append(merit[1])

        SW1.append(Sw[0])
        SW2.append(Sw[1])

        pypcaf_info.add_row(
            [num_div[i], dt] + iterations + deltas + T + idx_max +
            [region_order]
            )

    # Printing summary
    pypcaf_info.pprint(max_lines=-1, max_width=-1)
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
            EVALW1=EVALW1[idx],
            EVALW2=EVALW2[idx],
            MERIT1=MERIT1[idx],
            MERIT2=MERIT2[idx],
            SW1=SW1[idx],
            SW2=SW2[idx]
            )

    final_time = np.round((TIME.time() - start_time) / 60, 2)
    print(
        '\n **** PyPCAF Completed at {} mins **** \n'.format(final_time)
        )
