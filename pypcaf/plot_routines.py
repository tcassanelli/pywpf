#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import ascii
from scipy.constants import golden
from .pcaf_functions import flat_region_finder

__all__ = [
    'plot_waterfall', 'plot_lc', 'plot_period_double', 'plot_all_periods'
    ]

# Plot style added from relative path
plotstyle_dir = os.path.dirname(__file__)
plt.style.use(os.path.join(plotstyle_dir, 'pypcaf_sty.mplstyle'))


def plot_waterfall(waterfall, T):

    dt = T[1] - T[0]  # time bin
    M, N = waterfall.shape

    fig, ax = plt.subplots()
    im = ax.imshow(waterfall, cmap='viridis', interpolation='nearest')
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Total counts')
    ax.set_title('Waterfall rows: {}, dt = {} s'.format(M, dt))
    ax.set_xlabel('Bin s')
    ax.set_ylabel('Light curves')
    ax.grid('off')

    return fig


def plot_lc(lc, N, T):

    dt = T[1] - T[0]
    period_time = np.linspace(0, T, N + 1)[:-1]

    fig, ax = plt.subplots()
    ax.plot(
        period_time, lc, 'ro-',
        label='Period {} s'.format(T), linewidth=1.5
        )
    ax.set_title('Light curve dt = {} s'.format(dt))
    ax.set_xlabel('Time s')
    ax.set_ylabel('Total counts')
    ax.legend(loc='best')
    ax.grid('off')

    return fig


def plot_period_single(pypcaf_path, num_div, T_ref=None):

    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    pypcaf_out = os.path.join(
        pypcaf_path, 'M{}.npz'.format(num_div)
        )
    info = ascii.read(pypcaf_info)

    # index relative to the num_div
    idx = np.where(info['num_div'] == num_div)[0][0]

    # Estimated periods
    T_est = info['T_est'][idx]

    data = np.load(pypcaf_out)

    time_x = np.linspace(
        info['T_init'][idx] - info['delta'][idx] * info['iter'][idx] / 2,
        info['T_init'][idx] + info['delta'][idx] * info['iter'][idx] / 2,
        info['iter'][idx],
        endpoint=False
        )

    fig, ax = plt.subplots(figsize=(13, 13 / golden / 2))

    for i, alpha in zip(range(2), [1, .5]):

        # Eigenvalues 1 and 2
        ax.plot(
            time_x,
            data['EVALW'][:, i],
            label='$u_' + str(i + 1) + '$',
            color='b',
            alpha=alpha,
            marker='+'
            )

        # Scalar 1 and 2
        ax.plot(
            time_x,
            data['SW'][:, i],
            label='$s_' + str(i + 1) + '$',
            color='g',
            alpha=alpha,
            marker='^'
            )

    # Merit function
    ax.plot(
        time_x,
        data['MERIT'],
        label='$M$',
        color='r',
        linewidth=.5,
        marker='o'
        )

    # Virtical lines: T_est
    ax.axvline(
        x=info['T_init'][idx],
        color='darkorange',
        label='$T_\\mathrm{init}$',
        linewidth=.5
        )

    # Vertical lines
    ax.axvline(
        x=T_est,
        color='y',
        label='$T_\\mathrm{est1}$',
        linewidth=.5
        )

    # Adding reference period, optional
    if T_ref is not None:
        ax.axvline(
            x=T_ref,
            color='dimgray',
            label='$T_\\mathrm{ref}$',
            linewidth=.5
            )

    # Title & legends
    name = os.path.split(pypcaf_path)[1]
    ax.set_title(
        name + '. $T_\\mathrm{est}=' +
        str(T_est) + '$ s, $dt=' + str(info['dt'][idx]) +
        '$ s, $M=' + str(info['num_div'][idx]) + '$'
        )

    ax.set_xlabel('Time s')
    ax.set_ylabel('Scalar, eigenvalues and merit amplitude')
    ax.set_xlim(time_x[0], time_x[-1])
    ax.legend(loc='upper right')

    fig.tight_layout()

    return fig


def plot_period_double(pypcaf_path, num_div, T_ref=None):

    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    pypcaf_out = os.path.join(
        pypcaf_path, 'M{}.npz'.format(num_div)
        )
    info = ascii.read(pypcaf_info)

    # index relative to the num_div
    idx = np.where(info['num_div'] == num_div)[0][0]

    # Estimated periods
    T_est1 = info['T_est1'][idx]
    T_est2 = info['T_est2'][idx]

    data = np.load(pypcaf_out)

    time_x1 = np.linspace(
        info['T_init'][idx] - info['delta1'][idx] * info['iter1'][idx] / 2,
        info['T_init'][idx] + info['delta1'][idx] * info['iter1'][idx] / 2,
        info['iter1'][idx],
        endpoint=False
        )

    time_x2 = np.linspace(
        info['T_est1'][idx] - info['delta2'][idx] * info['iter2'][idx] / 2,
        info['T_est1'][idx] + info['delta2'][idx] * info['iter2'][idx] / 2,
        info['iter2'][idx],
        endpoint=False
        )

    time_x = [time_x1, time_x2]

    fig, ax = plt.subplots(nrows=2, figsize=(13, 13 / golden))

    for i, alpha in zip(range(2), [1, .5]):
        for j in range(2):

            # Eigenvalues 1 and 2
            ax[j].plot(
                time_x[j],
                data['EVALW{}'.format(j + 1)][:, i],
                label='$u_' + str(i + 1) + '$',
                color='b',
                alpha=alpha,
                marker='+'
                )

            # Scalar 1 and 2
            ax[j].plot(
                time_x[j],
                data['SW{}'.format(j + 1)][:, i],
                label='$s_' + str(i + 1) + '$',
                color='g',
                alpha=alpha,
                marker='^'
                )

        # Merit functions 1 and 2
        ax[i].plot(
            time_x[i],
            data['MERIT{}'.format(i + 1)],
            label='$M_' + str(i + 1) + '$',
            color='r',
            linewidth=.5,
            marker='o'
            )

    # Merit function 1 in second plot
    ax[1].plot(
        time_x[0],
        data['MERIT1'],
        label='$M_1$',
        color='r',
        linewidth=.3
        )

    # Virtical lines: T_est
    ax[0].axvline(
        x=info['T_init'][idx],
        color='darkorange',
        label='$T_\\mathrm{init}$',
        linewidth=.5
        )

    # Vertical lines
    ax[0].axvline(
        x=T_est1,
        color='y',
        label='$T_\\mathrm{est1}$',
        linewidth=.5
        )
    ax[1].axvline(
        x=T_est1,
        color='y',
        label='$T_\\mathrm{est1}$',
        linewidth=.5
        )
    ax[1].axvline(
        x=T_est2,
        color='m',
        label='$T_\\mathrm{est2}$',
        linewidth=.5
        )

    # Adding reference period, optional
    if T_ref is not None:
        for _ax in ax:
            _ax.axvline(
                x=T_ref,
                color='dimgray',
                label='$T_\\mathrm{ref}$',
                linewidth=.5
                )

    # Titles
    name = os.path.split(pypcaf_path)[1]
    ax[0].set_title(
        name + ' first iteration. $T_\\mathrm{est1}=' +
        str(T_est1) + '$ s, $dt=' + str(info['dt'][idx]) +
        '$ s, $M=' + str(info['num_div'][idx]) + '$'
        )
    ax[1].set_title(
        name + ' second iteration. $T_\\mathrm{est2}=' +
        str(T_est2) + '$ s, $dt=' + str(info['dt'][idx]) +
        '$ s, $M=' + str(info['num_div'][idx]) + '$'
        )

    for i in range(2):
        ax[i].set_xlabel('Time s')
        ax[i].set_ylabel('Scalar, eigenvalues and merit amplitude')
        ax[i].set_xlim(time_x[i][0], time_x[i][-1])
        ax[i].legend(loc='upper right')

    fig.tight_layout()

    return fig


def plot_all_periods(pypcaf_path, T_ref=None):

    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    info = ascii.read(pypcaf_info)
    num_div = info['num_div']

    # Printing summary
    print('... Making plots ... \n')
    info.pprint(max_lines=-1, max_width=-1)
    print('\n')

    path_plot = os.path.join(pypcaf_path, 'plots')
    if not os.path.exists(path_plot):
        os.makedirs(path_plot)

    # Iteration over all number of divisions in waterfall (M)
    if 'iter2' in info.keys():
        plot_func = plot_period_double
    else:
        plot_func = plot_period_single

    for M in num_div:
        fig = plot_func(
            pypcaf_path=pypcaf_path,
            num_div=M,
            T_ref=T_ref
            )
        fig.savefig(os.path.join(path_plot, 'M{}.pdf'.format(M)))
