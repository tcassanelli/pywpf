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


def new_fold(time, dt, period, num_div, plot_check=False):
    """
    Folding algorithm using the waterfall diagrams. Time is data in a the .csv
    file. It is a column vector
    num_div is the number of divisions made to the time array (aka data) or rows in waterfall The period
    will only be an approximation, needs to be iterated to correct it!

    Parameters
    ----------
    time : `~numpy.ndarray` or list
        Observed periodicity time with the telescope.
    dt : float
        Bintime.
    period : float
        Estimated period or staring period.
    num_div : int
        Number of divisions made to the time array or rows in waterfall diagram. Later defined as M. It is
        also the number of eigenvectors.
    plot_check: boolean
        To decide if it is necesary an eye inspection.

    Returns
    -------
    lc : `~numpy.ndarray`
        Light curve of the input period. It is a one column array.
    waterfall: `~numpy.ndarray`
        Matrix that contains an amount of num_div rows. The number of columns
        is Nint.
    """

    if period < dt:
        print('WARNING: Period cannot be smaller than bin size (dt)')

    # Length light-curve. It needs to be a division with no modulus
    # N represents the columns in the waterfall
    # It has to be chosen the int value over the approximation
    Nint = round(period / dt)
    dt = period / Nint  # dt recalculated so it becomes an interger

    # Period division in Nint*dt
    period_div_dt = np.linspace(0, period, Nint + 1)

    # number of samples that will be considered for each row of the waterfall
    num_samples = np.floor(time.size / num_div)

    # Modulus divions left. Return element-wise remainder of division
    remainder = np.mod(time, period)

    # for each line in the waterfall diagram
    for line in range(num_div):
        # selection of each num_div in time array
        indices = np.arange(
            num_samples * line, num_samples * (line + 1), dtype=int
            )
        # matrix that contains info for waterfall diagram
        if line == 0:
            waterfall = np.histogram(remainder[indices], bins=period_div_dt)[0]
        else:
            waterfall = np.vstack(
                (waterfall, np.histogram(remainder[indices],
                bins=period_div_dt)[0])
                )

    # Light-Curve plot
    lc = np.histogram(remainder, period_div_dt)[0]
    period_time_one = np.arange(0, period, dt)

    # Stacking two periods together for visualization
    lc2 = np.hstack((lc, lc))
    period_time_two = np.arange(0, 2 * period, dt)

    if plot_check:
        fig1, ax1 = plt.subplots()
        ax1.plot(period_time_two, lc2, 'ro-', label='Period ' + str(period) + ' s', linewidth=1.5)
        ax1.set_title('Light curve dt = ' + str(dt) + ' s')
        ax1.set_xlabel('Time s')
        ax1.set_ylabel('Total counts')
        ax1.legend(loc='best')
        ax1.grid()

        fig2, ax2 = plt.subplots()
        im2 = ax2.imshow(waterfall, cmap=plt.cm.jet, interpolation='nearest', aspect='auto')
        cb = fig2.colorbar(im2, ax=ax2)
        cb.set_label('Total counts')
        ax2.set_title('Waterfall rows: ' + str(num_div) + ', dt = ' + str(dt) + ' s')
        ax2.set_xlabel('Bin s')
        ax2.set_ylabel('Light curves')
        ax2.grid()

        plt.show()

    return lc, waterfall


def fast_pca(waterfall, plot_check=False):
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
    V, PC = np.linalg.eig(cov)

    V_sorted = np.sort(V.real)[::-1].tolist() # Eigenvalue
    j_indices = np.argsort(V.real)[::-1]
    PC_sorted = PC[:, j_indices] # Eigenvector or PCs

    signals = np.dot(PC_sorted.T, norm) # Information matrix, not clear whar represents!

    # Plot to visualize the PCs
    if plot_check:
        width = 0.8
        ind = np.arange(0, len(V_sorted))

        fig1, ax1 = plt.subplots()
        ax1.bar(ind, V_sorted, width=width)
        ax1.set_xlabel('Component value')
        ax1.set_ylabel('Eigenvalue amplitude')
        ax1.set_title('PCA values')
        ax1.set_xticks(ind + width/2)
        ax1.set_xticklabels(np.arange(1, len(V) + 1, dtype=int))
        ax1.grid()
        ax1.set_ylim([-0.1, V_sorted[0] + 0.1])
        ax1.set_xlim([-0.1, len(V)])

        fig2, ax2 = plt.subplots()
        im2 = ax2.imshow(signals, interpolation='nearest', aspect='auto')
        cb2 = fig2.colorbar(im2, ax=ax2)
        cb2.set_label('Norm(0, 1) counts')
        ax2.set_title('Signal = PC.T * normalized')
        ax2.set_xlabel('Bins s')
        ax2.set_ylabel('Light curves')
        ax2.grid()

        fig3, ax3 = plt.subplots()
        im3 = ax3.imshow(norm, interpolation='nearest', aspect='auto')
        cb3 = fig3.colorbar(im3, ax=ax3)
        cb3.set_label('Norm(0, 1) counts')
        ax3.set_title('Normalized waterfall')
        ax3.set_xlabel('Bins s')
        ax3.set_ylabel('Light curves')
        ax3.grid()

        plt.show()

    return V_sorted, PC_sorted, cov, norm, signals


def delta_finder(period, iterations, delta, time, dt, num_div):
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
    period_iter = period - iterations / 2 * delta

    VARIANCE = []
    SCALAR = []  # Scalar matrix
    unit_vec = np.ones((num_div, 1)) / np.sqrt(num_div)  # unitary vector

    for i in range(iterations):
        waterfall = new_fold(time, dt, period_iter, num_div)[1]
        eigenvalues, eigenvectors, _, _, _ = fast_pca(waterfall)

        # It is a vector with scalar_to_save = [s0, s1, s2, ...] for the num_div value
        scalar_to_save = np.sum(eigenvectors * unit_vec, axis=0).tolist()

        SCALAR.append(scalar_to_save)  # Both values are in decreasing order
        VARIANCE.append(eigenvalues)

        period_iter += delta

    S_array = np.abs(np.array(SCALAR))  # S_array[:, 0] represents all iteration for the first eigenvector
    V_array = np.array(VARIANCE)  # V_array[:, 0] represents all iteration for the first eigenvalue

    # Correspondent eigenvalue to the maximum selected scalar
    # V_corr = np.choose(np.argmax(S_array, axis=1), V_array.T)   # has a 32 lim!
    S_array_argmax = np.argmax(S_array, axis=1)

    V_corr = V_array[range(S_array_argmax.size), S_array_argmax]

    S_avg = [] # max scalar minus its average
    M = len(S_array[0])
    N = len(S_array[:, 0])
    for i in range(0, N):
        noise = np.sum(S_array[i]) - np.max(S_array[i])
        S_avg.append(np.max(S_array[i]) - noise / (M - 1))
    S_avg_array = np.array(S_avg)

    # (maximum scalar minus average) times the associated eigenvalue
    mstev = S_avg_array * V_corr # mstev = Maximum Scalar Times EigenValue

    max_idx = flat_region_finder(mstev.tolist())[0]

    period_final = period - iterations / 2 * delta + max_idx * delta

    return period_final, V_array, S_array, mstev, max_idx


def find_period(
    time, period, dt, num_div, iter1, delta1, iter2, delta2, noisy_signal=True
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

    if noisy_signal:
        period_start1 = period
    else:
        freq_start = pre_analysis(time, dt, period)[1]
        period_start1 = 1 / freq_start

    period_final1, var_iter1, scalar_iter1, mstev_iter1, max_index1 = \
    delta_finder(period_start1, iter1, delta1, time, dt, num_div)

    period_start2 = period_final1
    period_final2, var_iter2, scalar_iter2, mstev_iter2, max_index2 = \
    delta_finder(period_start2, iter2, delta2, time, dt, num_div)

    return [period, period_start1, period_final1, period_final2], [var_iter1, var_iter2], \
    [scalar_iter1, scalar_iter2], [mstev_iter1, mstev_iter2], [max_index1, max_index2]


if __name__ == '__main__':
    """
    Before running the full computation, pca_run, test the program with several iterations
    see how this behaves with the plots!
    """

    file_name = 'FILE NAME' # Contains the time array

    dt = 0.002793 # 4 ms, 0.002793 s
    period_start = 0.089367 # Initial period, usualy well known

    num_div = 20

    time = np.genfromtxt('data_pulsar/' + file_name + '.csv')
    lc, water = new_fold(time, dt, period_start, num_div, plot_check=True)

    V_sorted, PC_sorted, cov, norm, signals = fast_pca(water, True)
