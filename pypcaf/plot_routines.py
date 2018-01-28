import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

plt.style.use('pypcaf_sty.mplstyle')


def plot_waterfall(waterfall, T):

    dt = T[1] - T[0]
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


def plot_period(pcaf_path, num_div):

    # M, dt, iter1, iter2, delta1, delta2, T_init, T_est1, T_est2, idx1_max, idx2_max
    # [EValw1, EValw2], [Sw1, Sw2], [M1, M2]

    pcaf_info = os.path.join(pcaf_path, 'pcaf_info.dat')
    pcaf_out = os.path.join(pcaf_path, 'pcaf_out_M{}.npz'.format(num_div - 1))

    fig, ax = plt.subplots(nrows=2)

    info = ascii.read(pcaf_info)
    print(info)

    M = np.where(info['num_div'] == num_div)[0][0]

    data = np.load(pcaf_out)

    time_x1 = np.linspace(
        info['T_init'][M] - info['delta1'][M] * info['iter1'][M] / 2,
        info['T_init'][M] + info['delta1'][M] * info['iter1'][M] / 2,
        info['iter1'][M],
        endpoint=False
        )

    time_x2 = np.linspace(
        info['T_est1'][M] - info['delta2'][M] * info['iter2'][M] / 2,
        info['T_est1'][M] + info['delta2'][M] * info['iter2'][M] / 2,
        info['iter2'][M],
        endpoint=False
        )

    time_x = [time_x1, time_x2]


    # Eignevalues 1 and 2
    ax[0].plot(
        time_x[0],
        data['EVALW1'][:, 0],
        label='$u_1$',
        color='b',
        linewidth=.8,
        alpha=1,
        linestyle='-',
        marker='+'
        )
    ax[0].plot(
        time_x[0],
        data['EVALW1'][:, 1],
        label='$u_2$',
        color='b',
        linewidth=.8,
        alpha=.5,
        linestyle='-',
        marker='+'
        )

    ax[1].plot(
        time_x[1],
        data['EVALW2'][:, 0],
        label='$u_1$',
        color='b',
        linewidth=.8,
        alpha=1,
        linestyle='-',
        marker='+'
        )
    ax[1].plot(
        time_x[1],
        data['EVALW2'][:, 1],
        label='$u_2$',
        color='b',
        linewidth=.8,
        alpha=.5,
        linestyle='-',
        marker='+'
        )

    # Scalars 1 and 2
    ax[0].plot(
        time_x[0],
        data['SW1'][:, 0],
        label='$u_1$',
        color='g',
        linewidth=.8,
        alpha=1,
        linestyle='-',
        marker='^'
        )
    ax[0].plot(
        time_x[0],
        data['SW1'][:, 1],
        label='$u_2$',
        color='g',
        linewidth=.8,
        alpha=.5,
        linestyle='-',
        marker='^'
        )

    ax[1].plot(
        time_x[1],
        data['SW2'][:, 0],
        label='$u_1$',
        color='g',
        linewidth=.8,
        alpha=1,
        linestyle='-',
        marker='^'
        )
    ax[1].plot(
        time_x[1],
        data['SW2'][:, 1],
        label='$u_2$',
        color='g',
        linewidth=.8,
        alpha=.5,
        linestyle='-',
        marker='^'
        )

    # Virtical lines: T_est
    ax[0].axvline(
        x=info['T_est1'][M],
        color='y',
        label='$T_\\mathrm{est}$',
        linewidth=.8
        )
    ax[1].axvline(
        x=info['T_est2'][M],
        color='y',
        label='$T_\\mathrm{est}$',
        linewidth=.8
        )

    # Merit functions
    ax[0].plot(
        time_x[0],
        data['MERIT1'],
        label='$M_1$',
        color='r',
        linewidth=1,
        linestyle='-',
        marker='o'
        )
    ax[1].plot(
        time_x[0],
        data['MERIT1'],
        label='$M_1$',
        color='r',
        linewidth=0.3,
        linestyle='-',
        )
    ax[1].plot(
        time_x[1],
        data['MERIT2'],
        label='$M_2$',
        color='r',
        linewidth=1,
        linestyle='-',
        marker='o'
        )

    # Titles
    ax[0].set_title(
        'First iteration. $T={}$'.format(np.round(info['T_est1'][M], 3))
        )
    ax[1].set_title(
        'Second iteration. $T={}$'.format(np.round(info['T_est2'][M], 3))
        )

    for i in range(2):
        ax[i].grid()
        ax[i].set_xlabel('Time s')
        ax[i].set_ylabel('Scalar, eigenvalues and merit amplitude')
        ax[i].set_xlim(time_x[i][0], time_x[i][-1])
        ax[i].legend(loc='upper right')

    return fig


if __name__ == '__main__':

    plot_period(pcaf_path='./', num_div=2)
    plt.show()
