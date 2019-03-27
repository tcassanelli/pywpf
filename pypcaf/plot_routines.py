#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from astropy.io import ascii
from scipy.constants import golden
from .pcaf_functions import flat_region_finder


__all__ = [
    'plot_waterfall', 'plot_lc', 'plot_period_single',
    'plot_period_double', 'plot_all', 'plot_all_merit', 'plot_all_scalar',
    'plot_average_merit'
    ]

# Plot style added from relative path
plotstyle_dir = os.path.dirname(__file__)
plt.style.use(os.path.join(plotstyle_dir, 'pypcaf_sty.mplstyle'))

# color definitions
color_code = {
    'T_est1': 'brown',
    'T_est2': 'purple',
    'T_ref': 'dimgray',
    'T_init': 'gold',
    's': 'green',
    'u': 'darkblue',
    'merit': 'red'
    }


# This plot may not be useful anymore
def plot_waterfall(waterfall, T):

    dt = T[1] - T[0]  # time bin
    M, N = waterfall.shape

    fig, ax = plt.subplots()
    im = ax.imshow(waterfall, cmap='viridis', interpolation='nearest')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.03)
    cb = fig.colorbar(im, cax=cax)

    cb.ax.set_ylabel('Total counts')
    ax.set_title('Waterfall rows: {}, dt = {} [s]'.format(M, dt))
    ax.set_xlabel('Bin [s]')
    ax.set_ylabel('Light curves')
    ax.grid(False)

    return fig


# This plot may not be useful anymore
def plot_lc(lc, N, T):

    dt = T[1] - T[0]
    period_time = np.linspace(0, T, N + 1)[:-1]

    fig, ax = plt.subplots()
    ax.plot(
        period_time, lc, 'ro-',
        label='Period {} s'.format(T), linewidth=1.5
        )
    ax.set_title('Light curve dt = {} [s]'.format(dt))
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Total counts')
    ax.legend(loc='best')
    ax.grid(False)

    return fig


def plot_average_merit(pypcaf_path, T_ref=None):
    """Plot the average of a list of num_div elements"""

    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    info = ascii.read(pypcaf_info)

    # Change in the future, all of the merits have same time sample
    time_x = np.linspace(
        info['T_init'][0] - info['delta'][0] * info['iter'][0] / 2,
        info['T_init'][0] + info['delta'][0] * info['iter'][0] / 2,
        info['iter'][0],
        endpoint=False
        )

    # Average over these num_div (M)
    averages = [
        range(3, 7), range(3, 10), range(3, 13), range(3, 16), range(3, 21)
        ]

    # this needs to be fixed!
    merit_sum = np.zeros((len(averages), info['iter'][0]))
    for i, num_div in enumerate(averages):
        for M in num_div:
            pypcaf_out = os.path.join(pypcaf_path, 'M{}.npz'.format(M))
            data = np.load(pypcaf_out)
            merit_sum[i] += data['MERIT']

    fig = plt.figure(figsize=(10, 10 / golden))
    fig.subplots_adjust(wspace=0, hspace=0)
    gs = gridspec.GridSpec(len(averages), 1)

    for i, num_div in enumerate(averages):
        ax = fig.add_subplot(gs[i])

        merit_avg = merit_sum[i] / len(num_div)

        ax.plot(
            time_x,
            merit_avg / merit_avg.max(),  # average normalized
            color=color_code['merit'],
            linewidth=1,
            # marker='o'
            )

        ax.fill_between(
            time_x, 0,
            merit_avg / merit_avg.max(),
            color=color_code['merit'],
            alpha=.5
            )

        if T_ref is not None:
            ax.axvline(
                x=T_ref,
                color=color_code['T_ref'],
                label='$T_{\\mathrm{ref}}$',
                linewidth=0.5
                )

        ax.set_xlim(time_x[0], time_x[-1])
        ax.set_ylim(0, 1)
        ax.set_facecolor('whitesmoke')
        ax.legend(loc='upper right')

        ax.set_ylabel(
            '$\\left<M_{3\\mathrm{-}' + str(num_div[-1]) + '}(u, s)\\right>$'
            )
        if i == 0:
            name = os.path.basename(os.path.dirname(pypcaf_path + '/'))
            ax.set_title(name + '. Average merit function (normalized)')
        if i != len(averages) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Time [s]')

    fig.tight_layout()

    return fig


def plot_all_merit(pypcaf_path, T_ref=None):
    """
    Plot for all the merit functions computed up to num_div iterations of
    pypcaf.
    """
    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    info = ascii.read(pypcaf_info)
    Ms = info['num_div']

    time_x = np.linspace(
        info['T_init'][0] - info['delta'][0] * info['iter'][0] / 2,
        info['T_init'][0] + info['delta'][0] * info['iter'][0] / 2,
        info['iter'][0],
        endpoint=False
        )

    fig, ax = plt.subplots(figsize=(12, 12 / golden))

    ax.axvline(
        x=info['T_init'][0],
        color=color_code['T_init'],
        label='$T_{\\mathrm{init}}$',
        linewidth=.5
        )

    if T_ref is not None:
        ax.axvline(
            x=T_ref,
            color=color_code['T_ref'],
            label='$T_{\\mathrm{ref}}$',
            linewidth=.5
            )

    label_T_est = '$T_{\\mathrm{est}}$'
    label_merit = '$M_{3\\mathrm{-}20}(u, s)$'

    for i, M in enumerate(Ms):

        pypcaf_out = os.path.join(
            pypcaf_path, 'M{}.npz'.format(M)
            )
        data = np.load(pypcaf_out)

        merit_norm = data['MERIT'] / data['MERIT'].max()

        ax.plot(
            time_x,
            merit_norm + i,
            color=color_code['merit'],
            label=label_merit,
            linewidth=1
            )
        label_merit = '_nolegend_'

        ax.fill_between(
            time_x, np.ones(time_x.size) * i,
            merit_norm + i,
            color=color_code['merit'],
            alpha=.5
            )

        # Estimated periods
        T_est = info['T_est'][i]

        ax.axvline(
            x=T_est,
            ymin=i / (Ms[-1] - 2),
            ymax=(i + 1) / (Ms[-1] - 2),
            color=color_code['T_est1'],
            label=label_T_est,
            linewidth=.5
            )
        label_T_est = '_nolegend_'

        # black horizontal lines
        ax.axhline(
            y=i,
            color='k',
            linewidth=0.2
            )

    ax.set_yticks(np.arange(Ms[0] - 4, Ms[-1] - 2, 1) + .5)
    ax.set_yticklabels(np.arange(Ms[0] - 1, Ms[-1] + 1, 1))

    ax.set_ylabel('Number of divisions ($M$)')
    ax.set_xlabel('Time [s]')

    # careful with this, it is only a work around!
    name = os.path.basename(os.path.dirname(pypcaf_path + '/'))
    ax.set_title(
        name + '. Merit functions (normalized) $' +
        str(Ms[0]) + '$-$' + str(Ms[-1]) + '$'
        )

    ax.grid(True)
    ax.set_xlim(time_x[0], time_x[-1])
    ax.set_ylim(0, Ms[-1] - 2)
    ax.set_facecolor('whitesmoke')
    ax.legend(loc='upper right')

    fig.tight_layout()

    return fig


def plot_all_eigenvalue(pypcaf_path, num_div, T_ref=None):

    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    pypcaf_out = os.path.join(pypcaf_path, 'M{}.npz'.format(num_div))
    info = ascii.read(pypcaf_info)

    # index relative to the num_div
    idx = np.where(info['num_div'] == num_div)[0][0]

    # Estimated periods
    T_est = info['T_est'][idx]

    data = np.load(pypcaf_out)

    time_x = np.linspace(
        info['T_init'][0] - info['delta'][0] * info['iter'][0] / 2,
        info['T_init'][0] + info['delta'][0] * info['iter'][0] / 2,
        info['iter'][0],
        endpoint=False
        )

    if 'iter2' in info.keys():
        EVALW = data['EVALW2']
    else:
        EVALW = data['EVALW']

    fig, ax = plt.subplots(figsize=(12, 12 / golden))

    label_u_i = '$u_i$'
    for i in range(num_div):
        EVAL_norm = EVALW[:, i] / EVALW[:, 0].max()
        shift = i

        ax.plot(
            time_x,
            EVAL_norm + shift,
            color=color_code['u'],
            linewidth=1,
            label=label_u_i
            )
        label_u_i = '_nolegend_'

        ax.fill_between(
            x=time_x,
            y1=np.ones(time_x.size) * shift,
            y2=EVAL_norm + shift,
            color=color_code['u'],
            alpha=.5
            )

    # Estimated and reference period
    ax.axvline(
        x=T_est,
        color=color_code['T_est1'],
        label='$T_\\mathrm{est1}$',
        linewidth=.5
        )

    if T_ref is not None:
        ax.axvline(
            x=T_ref,
            color=color_code['T_ref'],
            label='$T_\\mathrm{ref}$',
            linewidth=.5
            )

    ax.axvline(
        x=info['T_init'][0],
        color=color_code['T_init'],
        label='$T_\\mathrm{init}$',
        linewidth=.5
        )

    # black horizontal lines
    for i in range(num_div):
        ax.axhline(
            y=i,
            color='k',
            linewidth=0.2
            )

    ax.grid(True)
    ax.set_xlim(time_x[0], time_x[-1])
    ax.set_ylim(0, num_div)

    ax.set_yticks(np.arange(0, num_div, 1) + .5)
    ax.set_yticklabels(np.arange(1, num_div + 1, 1))

    ax.set_ylabel('Eigenvalue number $(u_i)$')
    ax.set_xlabel('Time [s]')

    name = os.path.split(pypcaf_path)[1]
    ax.set_title(
        name +
        '. Normalized eigenvalues $u_i/\\mathrm{max}(u_1)$ $T_\\mathrm{est}=' +
        str(T_est) + '$ s, $dt=' + str(info['dt'][idx]) +
        '$ s, $M=' + str(info['num_div'][idx]) + '$'
        )
    ax.set_facecolor('whitesmoke')
    ax.legend(loc='upper right')

    fig.tight_layout()

    return fig


def plot_all_scalar(pypcaf_path, num_div, T_ref=None):
    """
    for the moment only for a single iteration NO pcaf_double output
    """

    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    pypcaf_out = os.path.join(pypcaf_path, 'M{}.npz'.format(num_div))
    info = ascii.read(pypcaf_info)

    # index relative to the num_div
    idx = np.where(info['num_div'] == num_div)[0][0]

    # Estimated periods
    T_est = info['T_est'][idx]

    data = np.load(pypcaf_out)

    time_x = np.linspace(
        info['T_init'][0] - info['delta'][0] * info['iter'][0] / 2,
        info['T_init'][0] + info['delta'][0] * info['iter'][0] / 2,
        info['iter'][0],
        endpoint=False
        )

    if 'iter2' in info.keys():
        SW = data['SW2']
    else:
        SW = data['SW']

    fig, ax = plt.subplots(figsize=(12, 12 / golden))

    label_s_i = '$s_i$'
    for i in range(num_div):
        ax.plot(
            time_x,
            SW[:, i] + i,
            color=color_code['s'],
            linewidth=1,
            label=label_s_i
            )
        label_s_i = '_nolegend_'

        ax.fill_between(
            x=time_x,
            y1=np.ones(time_x.size) * i,
            y2=SW[:, i] + i,
            color=color_code['s'],
            alpha=.5
            )

    # Estimated and reference period
    ax.axvline(
        x=T_est,
        color=color_code['T_est1'],
        label='$T_\\mathrm{est1}$',
        linewidth=.5
        )

    if T_ref is not None:
        ax.axvline(
            x=T_ref,
            color=color_code['T_ref'],
            label='$T_\\mathrm{ref}$',
            linewidth=.5
            )

    ax.axvline(
        x=info['T_init'][0],
        color=color_code['T_init'],
        label='$T_\\mathrm{init}$',
        linewidth=.5
        )

    # black horizontal lines
    for i in range(num_div):
        ax.axhline(
            y=i,
            color='k',
            linewidth=0.2
            )

    ax.grid(True)
    ax.set_xlim(time_x[0], time_x[-1])
    ax.set_ylim(0, num_div)

    ax.set_yticks(np.arange(0, num_div, 1) + .5)
    ax.set_yticklabels(np.arange(1, num_div + 1, 1))

    ax.set_ylabel('Scalar number $(s_i)$')
    ax.set_xlabel('Time [s]')

    name = os.path.split(pypcaf_path)[1]
    ax.set_title(
        name + '. $T_\\mathrm{est}=' +
        str(T_est) + '$ s, $dt=' + str(info['dt'][idx]) +
        '$ s, $M=' + str(info['num_div'][idx]) + '$'
        )
    ax.set_facecolor('whitesmoke')
    ax.legend(loc='upper right')

    fig.tight_layout()

    return fig


def plot_period_single(pypcaf_path, num_div, T_ref=None, merit_func=None):

    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    pypcaf_out = os.path.join(
        pypcaf_path, 'M{}.npz'.format(num_div)
        )
    info = ascii.read(pypcaf_info)

    # index relative to the num_div
    idx = np.where(info['num_div'] == num_div)[0][0]

    data = np.load(pypcaf_out)

    time_x = np.linspace(
        info['T_init'][0] - info['delta'][0] * info['iter'][0] / 2,
        info['T_init'][0] + info['delta'][0] * info['iter'][0] / 2,
        info['iter'][0],
        endpoint=False
        )

    fig, ax = plt.subplots(figsize=(13, 13 / golden / 2))

    # Estimated periods
    if merit_func is not None:
        merit = merit_func(EValw=data['EVALW'], Sw=data['SW'])
        region_order = info['region_order'][idx]
        iteration = info['iter'][idx]
        delta = info['delta'][idx]
        T_init = info['T_init'][idx]

        if region_order > 1:
            idx_max = flat_region_finder(X=merit, n=region_order)
        else:
            idx_max = np.argmax(merit)

        T_est = T_init - iteration / 2 * delta + idx_max * delta
    else:
        T_est = info['T_est'][idx]
        merit = data['MERIT']

    # for i, alpha in enumerate([1, .5]):

    #     # Scalar 1 and 2
    #     ax.plot(
    #         time_x,
    #         data['SW'][:, i],
    #         label='$s_' + str(i + 1) + '$',
    #         color='g',
    #         alpha=alpha,
    #         marker='^'
    #         )

    # Scalar average Sw[:1]
    ax.plot(
        time_x,
        data['SW'][:, 1:].mean(axis=1),
        label='$\\left<s_{2\\mathrm{-}' + str(num_div) + '}\\right>$',
        color=color_code['s'],
        alpha=.5,
        marker='^'
        )

    # Scalar 1
    ax.plot(
        time_x,
        data['SW'][:, 0],
        label='$s_1$',
        color=color_code['s'],
        alpha=1,
        marker='^'
        )

    # Eigenvalue 1
    ax.plot(
        time_x,
        data['EVALW'][:, 0],
        label='$u_1$',
        color=color_code['u'],
        alpha=1,
        marker='+'
        )

    # Merit function
    ax.plot(
        time_x,
        merit,
        label='$M(u, s)$',
        color=color_code['merit'],
        linewidth=1,
        marker='o'
        )

    # Virtical lines: T_est
    ax.axvline(
        x=info['T_init'][0],
        color=color_code['T_init'],
        label='$T_\\mathrm{init}$',
        linewidth=.5
        )

    # Vertical lines
    ax.axvline(
        x=T_est,
        color=color_code['T_est1'],
        label='$T_\\mathrm{est1}$',
        linewidth=.5
        )

    # Adding reference period, optional
    if T_ref is not None:
        ax.axvline(
            x=T_ref,
            color=color_code['T_ref'],
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

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Scalar, eigenvalues and merit amplitude')
    ax.set_xlim(time_x[0], time_x[-1])
    ax.legend(loc='upper right')
    ax.set_facecolor('whitesmoke')

    fig.tight_layout()

    return fig


# Check all modifications added to other plots and update this one!
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

    for i, alpha in enumerate([1, .5]):
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
        ax[i].set_xlabel('Time [s]')
        ax[i].set_ylabel('Scalar, eigenvalues and merit amplitude')
        ax[i].set_xlim(time_x[i][0], time_x[i][-1])
        ax[i].legend(loc='upper right')

    fig.tight_layout()

    return fig


def plot_all(pypcaf_path, T_ref=None):

    pypcaf_info = os.path.join(pypcaf_path, 'info.dat')
    info = ascii.read(pypcaf_info)
    num_div = info['num_div']

    name = os.path.split(pypcaf_path)[1]

    # Printing summary
    print('... Plotting PyPCAF results {} ... \n'.format(name))
    # info.pprint(max_lines=-1, max_width=-1)
    # print('\n')

    path_plot = os.path.join(pypcaf_path, 'plots')
    if not os.path.exists(path_plot):
        os.makedirs(path_plot)

    # Iteration over all number of divisions in waterfall (M)
    if 'iter2' in info.keys():
        plot_func = plot_period_double
    else:
        plot_func = plot_period_single

    for M in num_div:
        fig1 = plot_func(
            pypcaf_path=pypcaf_path,
            num_div=M,
            T_ref=T_ref
            )
        fig2 = plot_all_scalar(
            pypcaf_path=pypcaf_path,
            num_div=M,
            T_ref=T_ref
            )
        fig3 = plot_all_eigenvalue(
            pypcaf_path=pypcaf_path,
            num_div=M,
            T_ref=T_ref
            )

        fig1.savefig(os.path.join(path_plot, 'M{}.pdf'.format(M)))
        fig2.savefig(os.path.join(path_plot, 'SW{}.pdf'.format(M)))
        fig3.savefig(os.path.join(path_plot, 'EVALW{}.pdf'.format(M)))

        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)

    fig4 = plot_all_merit(pypcaf_path=pypcaf_path, T_ref=T_ref)
    fig5 = plot_average_merit(pypcaf_path=pypcaf_path, T_ref=T_ref)

    fig4.savefig(os.path.join(path_plot, 'M_all.pdf'))
    fig5.savefig(os.path.join(path_plot, 'M_average.pdf'))

    plt.close(fig4)
    plt.close(fig5)
