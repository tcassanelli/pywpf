from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# PRINCIPAL FUNCTIONS FOR PCA FOLDING
# AUTHOR: TOMAS CASSANELLI

# Collection of all functions that are meant to be used in pca_run.


def nextpow2(n):
    """
    Use nextpow2 to pad the signal you pass to FFT.

    Parameters
    ----------
    n : `int`

    Returns
    -------
    m : `int`
        Exponent of next higher power of 2.
    """

    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    m = int(np.log2(2 ** m_i))

    return m


def flat_region_finder(X, n=3):
    """
    To find the best location for the period a fine search has to be done, not
    only where the maximum is located but also where there is a plateau of
    maximums. This becomes clear when the merit function is plotted.

    Parameters
    ----------
    X : `~numpy.ndarray`
        Array to find maximum region, it is in the particular case the mstev (
        maximum to scalar eigenvalue).
    n : `int`
        Number of division that the plateau has to have, for observed data 3
        is more or less correct.

    Returns
    -------
    idx_final : `int`
        Index of the position of the maximum point between a plateau of n
        points. Eventually this is the maximum
        period position.
    maximum : `float`
        The maximum number chosen between the three local maximums.
    """

    flat_region = [sum(X[i:i + n]) for i in range(X.size - n + 1)]

    idx = np.argmax(flat_region)  # finding the maximum index

    # find the maximum between the n points
    idx_max = idx + np.argmax([X[idx], X[idx + 1], X[idx + 2]])

    return idx_max


def pre_analysis(time, dt, period, plot_check=False):
    """
    Given an initial frequncy, finds a better one using FFT.

    Parameters
    ----------
    time : `~numpy.ndarray` or list
        Observed periodicity time with the telescope.
    dt : float
        Bintime.
    period : float
        Estimated period or staring period.
    plot_check: boolean
        To decide if it is necesary an eye inspection.

    Returns
    -------
    bin_data : `~numpy.ndarray`
        Return the bin edges (length(hist)+1) from the computed histogram of a
        set of data time.
    frquency: float
        Returns the approximated frquency after an FFT search.
    """

    freq_start = 1 / period

    # Setting the x axis for histogram, represents the time discrete position
    time_axis = np.arange(0, time[-1], dt)

    # counts the number of values in time that are within each specified bin
    # range, time_axis
    bin_data = np.histogram(time, bins=time_axis)[0]

    fs = 1 / dt  # frequency step

    # Fast Fourier Transform computation
    NFFT = 2 ** nextpow2(bin_data.size)  # Length of the transformed axis of the output
    y = np.fft.fft(bin_data, NFFT)  # Computed FFT
    N = NFFT // 2 + 1  # indices to erase the mirror effect from FFT
    Y = np.abs(y[:N])  # cleaned from all mirror effects
    freq_axis = fs / 2 * np.linspace(0, 1, N)

    # To give a zero value to the first components, due to FFT
    k = 100
    for i in np.arange(0, k, 1):
        Y[i] = 0

    start = int(len(freq_axis) * freq_start * 2 / fs) - 1
    stop = int(len(freq_axis) * freq_start * 2 / fs * 2)
    Y_selected = Y[np.arange(start, stop, dtype=int)]

    # Selection of the index in the freq_axis array
    index = np.argmax(Y_selected)

    frequency = freq_axis[index + int(len(freq_axis) * freq_start * 2 / fs)]

    if plot_check:

        fig1, ax1 = plt.subplots()
        ax1.hist(bin_data, histtype='stepfilled')
        ax1.set_title('Histogram dt = ' + str(dt))
        ax1.set_ylabel('Photon counts')
        ax1.set_xlabel('Time in ' + str(dt) + ' s units')
        ax1.grid()

        fig2, ax2 = plt.subplots()
        ax2.plot(freq_axis, Y)
        ax2.set_title('FFT binned data')
        ax2.set_ylabel('Amplitude')
        ax2.set_xlabel('Frequency Hz')
        ax2.grid()

        plt.show()

    return bin_data, frequency


def folding(time, dt, T_estimated, num_div):
    """
    Folding algorithm using the waterfall diagrams. Time is data in a the .csv
    file. It is a column vector ``num_div`` is the number of divisions made to
    the time array (a.k.a. data) or rows in waterfall The period will only be
    an approximation, needs to be iterated to correct it!

    Parameters
    ----------
    time : `~numpy.ndarray` or list
        Observed periodicity time with the telescope.
    dt : `float`
        Bintime.
    T_estimated : `float`
        Estimated period or staring period.
    num_div : `int`
        Number of divisions made to the time array or rows in waterfall
        diagram. Later defined as M. It is also the number of eigenvectors.
    plot_check : `bool`
        If `True` will plot.

    Returns
    -------
    lc : `~numpy.ndarray`
        Light curve of the input period. It is a one column array.
    waterfall : `~numpy.ndarray`
        Matrix that contains an amount of num_div rows. The number of columns
        is N.
    """

    if T_estimated < dt:
        print('WARNING: Period cannot be smaller than bin size (dt)')

    # Length light-curve. It needs to be a division with no modulus
    # N represents the columns in the waterfall
    # It has to be chosen the int value over the approximation
    M = num_div
    N = round(T_estimated / dt)

    # Period array with dt_int step and N + 1 elements
    T = np.linspace(0, T_estimated, N + 1)

    # number of samples that will be considered for each row of the waterfall
    ns = time.size // M

    # Modulus divions left. Return element-wise remainder of division
    remainder = time % T_estimated

    lc = np.histogram(remainder, T)[0]  # light-curve

    w = []
    for line in range(M):
        # selection of each M in time array
        indices = range(ns * line, ns * (line + 1))
        # matrix that contains info for waterfall diagram
        w.append(np.histogram(remainder[indices], bins=T)[0])
    waterfall = np.array(w)

    return lc, waterfall


def pca(waterfall):
    """
    Finds PCs, eigenvalues and signal matrix from waterfall M x N matrix.
    M: rows, number of segments in which the whole adquisition has been
    divided.
    N: columns, number of bins in folding period. Number of eigenvalues.

    Parameters
    ----------
    waterfall : `~numpy.ndarray`
        Matrix that contains an Amount of num_div rows (M). The number of columns is N.
    plot_check: boolean
        To decide if it is necesary an eye inspection.

    Returns
    -------
    V_sorted : `list`
        Eigenvalues from the covariance of the normalized waterfall matrix. Means mu = 0 and std = 1.
        Organized in decreasent order. Also known as variance.
    PC_sorted : `~numpy.ndarray`
        Unitary eigenvectors from the covariance of the normalized waterfall matrix that correspond to each
        eigenvalue. Sorted in the same way as the V_sorted list. One column, PC_sorted[:, i] represents one PC.
    cov : `~numpy.ndarray`
        Covariance of the normalized waterfall matrix.
    norm : `~numpy.ndarray`
        Normalization of the waterfall matrix. This is done for each row, (x - <x>)/std(x).
    signal : `~numpy.ndarray`
        Is the transposed PC_sorted times the normalized waterfall matrix.
    """

    M, N = waterfall.shape  # This should be the waterfall matrix
    mean = np.mean(waterfall, axis=1).reshape(M, 1)

    # carful, different in matlab!
    std = np.std(waterfall, axis=1, ddof=1).reshape(M, 1)

    # Normalization waterfall matrix to mean=0 and std=1
    norm = (waterfall - mean) / std  # (x - <x>)/std(x)
    # Covariance matrix
    cov = 1 / (N - 1) * np.dot(norm, norm.T)

    # Eigenvalue, Eigenvector
    # PC[:, i] is the eigenvector corresponding to V[i] eigenvalue
    EVal, EVec = np.linalg.eig(cov)

    # Sorting the eigenvalues and vectors
    EVal_sorted = np.sort(EVal.real)[::-1]             # eigenvalues
    EVec_sorted = EVec[:, np.argsort(EVal.real)[::-1]]  # eigenvectors or PCs
    K = np.dot(EVec_sorted.T, norm)                # information (signal) matrix

    # # Plot to visualize the PCs
    # if plot_check:
    #     width = 0.8
    #     ind = np.arange(0, len(V_sorted))

    #     fig1, ax1 = plt.subplots()
    #     ax1.bar(ind, V_sorted, width=width)
    #     ax1.set_xlabel('Component value')
    #     ax1.set_ylabel('Eigenvalue amplitude')
    #     ax1.set_title('PCA values')
    #     ax1.set_xticks(ind + width/2)
    #     ax1.set_xticklabels(np.arange(1, len(V) + 1, dtype=int))
    #     ax1.grid()
    #     ax1.set_ylim([-0.1, V_sorted[0] + 0.1])
    #     ax1.set_xlim([-0.1, len(V)])

    #     fig2, ax2 = plt.subplots()
    #     im2 = ax2.imshow(S, interpolation='nearest', aspect='auto')
    #     cb2 = fig2.colorbar(im2, ax=ax2)
    #     cb2.set_label('Norm(0, 1) counts')
    #     ax2.set_title('Signal = PC.T * normalized')
    #     ax2.set_xlabel('Bins s')
    #     ax2.set_ylabel('Light curves')
    #     ax2.grid()

    #     fig3, ax3 = plt.subplots()
    #     im3 = ax3.imshow(norm, interpolation='nearest', aspect='auto')
    #     cb3 = fig3.colorbar(im3, ax=ax3)
    #     cb3.set_label('Norm(0, 1) counts')
    #     ax3.set_title('Normalized waterfall')
    #     ax3.set_xlabel('Bins s')
    #     ax3.set_ylabel('Light curves')
    #     ax3.grid()

    #     plt.show()

    return EVal_sorted, EVec_sorted, K


def delta_finder(T_init, iterations, delta, time, dt, num_div, merit_func):
    """
    Finds the best period given an initial starting point, a number of iterations and a step to look for.
    It is the most inportant function which define the method of selection and the merit function of the script!

    Parameters
    ----------
    period : float
        Estimated period or staring period.
    iterations : int
        Interger number to iterate the main function loop.
    delta : float
        Increase of the period in each iteration. The orther of it is between 1e-7 - 4e-9.
    time : `~numpy.ndarray` or list
        Observed periodicity time with the telescope.
    dt : float
        Bintime.
    num_div : int
        Number of divisions made to the time array or rows in waterfall diagram. Later defined as M. It is
        also the number of eigenvectors.

    Returns
    -------
    period_final : float
        Optimum period of the iteration.
    V_array : `~numpy.ndarray`
        Values of all the eigenvalues, expressed as a np.array. i. e. V_array[:, 0] contains all the
        eigenvalues of the first position, or maximum eigenvalue. It has a length of the number of iterations.
    S_array : `~numpy.ndarray`
        Values of the first three scalars, expressed as a nu.array. i. e. S_array[:, 0] contains all the
        scalars of the first position. It has a length of the number of iterations. It is computed from the
        result of the hyperdimensional unitary vector times the eigenvalues (dot product), then the maximum
        absolute value per iteration is chosen.
    mstev : `~numpy.ndarray`
        Maximum scalar times the (selected) eigenvalue. It is the merit function selected to choose the right
        iteration. Represents the maximum scalar in a row (or in one iteration) minus the average of all the rest
        in the same interation. Then is multiplicated by the associated eigenvalue from the maximum scalar selected.
    max_idx : int
        For each iteration a step is added to the final period, this number of steps selected is the maximum index.
        Notice that the period search starts from (period - iterations / 2 * delta).
    """
    # makes an interval from central period, [period - i/2 * delta, period + i/2 * delta]

    M = num_div
    T_iterated = T_init - iterations / 2 * delta

    eigenvalues = []
    scalar = []  # Scalar matrix
    u = np.ones((M, 1)) / np.sqrt(M)  # hyperdiagonal unitary vector

    for i in range(iterations):
        waterfall = folding(
            time=time, dt=dt, T_estimated=T_iterated, num_div=M
            )[1]

        EVal, EVec, K = pca(waterfall=waterfall)

        # It is a vector with scalar_to_save = [s0, s1, s2, ...] for the num_div value
        scalar.append(np.sum(EVec * u, axis=0))
        # Both values are in decreasing order

        eigenvalues.append(EVal)

        T_iterated += delta

    # scalar and eigenvalues from waterfall N x M matrix
    Sw = np.abs(scalar)  # Sw[:, 0] represents all iteration for the first eigenvector

    EValw = np.array(eigenvalues)  # EValw[:, 0] represents all iteration for the first eigenvalue

    merit = merit_func(EValw=EValw, Sw=Sw)

    idx_max = flat_region_finder(merit)

    T_estimated = T_init - iterations / 2 * delta + idx_max * delta

    return T_estimated, EValw, Sw, merit, idx_max


def find_period(
    time, T_init, dt, num_div, iter1, delta1, iter2, delta2, merit_func, noisy_signal=True
        ):
    """
    Finds the optimal period using PCA. Encapsulates two iterations in one.

    Parameters
    ----------
    time : `~numpy.ndarray` or list
        Observed periodicity time with the telescope.
    period : float
        Estimated period or staring period.
    dt : float
        Bintime.
    num_div : int
        Number of divisions made to the time array or rows in waterfall diagram. Later defined as M. It is
        also the number of eigenvectors.
    iter1 : int
        Interger number to iterate the delta_finder function. Usually with a value of 100.
    delta1 : float
        Increase of the period in each iter1. The orther of it is between 1e-7.
    iter2 : int
        Interger number to iterate the delta_finder function. Usually with a value of 500.
    delta2 : float
        Increase of the period in each iter2. The orther of it is between 4e-9.
    noisy_signal : boolean
        If True the first iteration will be made using the function pre_analysis which looks for the
        best FFT frequency.

    Returns
    -------
    1/freq : float
        Initial given period of search. Check in literature of every object to find a good start.
    period_start1 : float
        Initial period of search in the case of a noisy signal. It first looks for the FFT.
    period_final1, 2 : float
        Best period from the first and second iteration. The starting period of the second iteration
        corresponds to period_final2.
    var_iter1, 2 : `~numpy.ndarray`
        See V_array in delta_finder function. 1 and 2 for first and second iterations.
    scalar_iter1, 2 : `~numpy.ndarray`
        See S_array in delta_finder function. 1 and 2 for first and second iterations.
    mstev_iter1, 2 : `~numpy.ndarray`
        See mstev in delta_finder function. 1 and 2 for first and second iterations.
    max_index1,, 2 : int
        See max_idx in delta_finder function. 1 and 2 for first and second iterations.
    """

    # if noisy_signal:
    #     period_start1 = period
    # else:
    #     freq_start = pre_analysis(time, dt, period)[1]
    #     period_start1 = 1 / freq_start

    (
        T_est1, EValw1, Sw1, M1, idx1_max
        ) = delta_finder(
        T_init=T_init,
        iterations=iter1,
        delta=delta1,
        time=time,
        dt=dt,
        num_div=num_div,
        merit_func=merit_func
        )

    (
        T_est2, EValw2, Sw2, M2, idx2_max
        ) = delta_finder(
        T_init=T_est1,
        iterations=iter2,
        delta=delta2,
        time=time,
        dt=dt,
        num_div=num_div,
        merit_func=merit_func
        )

    T = [T_init, T_est1, T_est2]

    return T, [idx1_max, idx2_max], [EValw1, EValw2], [Sw1, Sw2], [M1, M2]


if __name__ == '__main__':
    """
    Before running the full computation, pca_run, test the program with several iterations
    see how this behaves with the plots!
    """

    # file_name = 'FILE NAME' # Contains the time array

    # dt = 0.002793  # 4 ms, 0.002793 s
    # period_start = 0.089367  # Initial period, usualy well known

    # num_div = 20

    # time = np.genfromtxt('data_pulsar/' + file_name + '.csv')
    # lc, water = folding(time, dt, period_start, num_div, plot_check=True)

    # V_sorted, PC_sorted, cov, norm, signals = pca(water, True)

    # Testing the functions
    dt = 0.002793
    period = 0.089367
    num_div = 20

    time = np.genfromtxt('../../data_pulsar/rmr0.csv')

    lc, waterfall = folding(time, dt, period, num_div, plot_check=True)

    print('lc.shape: ', lc.shape)
    print('waterfall.shape: ', waterfall.shape)
