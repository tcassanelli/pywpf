#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import numpy as np
from fast_histogram import histogram1d


__all__ = [
    'nextpow2', 'flat_region_finder', 'pre_analysis', 'folding',
    'folding_fast', 'pca', 'find_period', 'cp_value', 'rms_value'
    ]


def nextpow2(n):
    """
    Use `pypcaf.nextpow2` to pad the signal you pass to FFT.

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
    Finds the highest (or maximum) plateau region in data, the plateau has
    length ``n``, then among the ``n`` points returns its maximum.

    Parameters
    ----------
    X : `~numpy.ndarray`
        One dimensional array where the maximum needs to be found.
    n : `int`
        Number of points in the plateau region. It is also named
        ``region_order``.

    Returns
    -------
    idx_max : `int`
        Index or position in the ``X`` array where the maximum of a plateau of
        length ``n`` is located.
    """

    flat_region = [sum(X[i:i + n]) for i in range(X.size - n + 1)]

    idx = np.argmax(flat_region)  # finding the maximum index

    # find the maximum between the n points
    idx_max = idx + np.argmax([X[idx], X[idx + 1], X[idx + 2]])

    return idx_max


def pre_analysis(times, dt, T_init):
    """
    Given a one dimensional times array, ``times``, the function computes its
    FFT to find a better period than the given one, ``T_init``. Be aware that
    depending on the signal-to-noise of the times array, may or may not find a
    correct value. For really noisy signals do not use this function.

    Parameters
    ----------
    times : `~numpy.ndarray`
        One dimensional times array with certain hidden periodicity, e.g.
        pulsar period.
    dt : `float`
        times bin.
    T_init : `float`
        Initial period.

    Returns
    -------
    T_est : `float`
        Estimated period from the ``times`` array, found by a FFT analysis.
    """

    if T_init < dt:
        raise TypeError(
            'Initial period (T_init) cannot be smaller than times bin (dt)'
            )

    f_start = 1 / T_init  # starting frequency

    # Setting the x axis for histogram
    times_n = int(round(times[-1] / dt) + 1)  # number of elements in times array
    times_x = np.linspace(0, times[-1], times_n)

    # counts the number of values in the times array that are within each
    # specified bin range, times_x
    hist = np.histogram(times, bins=times_x)[0]

    f_step = 1 / dt  # frequency step
    NFFT = 2 ** nextpow2(hist.size)  # length of the transform
    # power of two for ease of calculation

    freq_raw = np.fft.fft(a=hist, n=NFFT)     # Computed FFT
    freq_abs = np.abs(freq_raw[:NFFT // 2 + 1])   # Erasing mirror effect
    freq_axis = f_step / 2 * np.linspace(0, 1, NFFT // 2 + 1)

    # Erasing the first 100 elements from FFT usually they are noisy
    for i in range(100):
        freq_abs[i] = 0

    start = freq_axis.size * f_start * 2 // f_step - 1
    stop = freq_axis.size * f_start * 2 // f_step * 2

    # Selection of the index in the freq_axis array
    idx = np.argmax(freq_abs[range(start, stop)])
    f_est = freq_axis[idx + freq_axis.size * f_start * 2 // f_step]
    T_est = 1 / f_est

    return T_est


def rms_value(X):
    """
    Computes the root-mean-square value for an array.

    Parameters
    ----------
    X : `~numpy.ndarray`
        One or two dimensional array.
    """

    return np.sqrt(X.dot(X) / X.size)


def cp_value(merit, idx_max, frac=0.25):
    """
    Global value signal-to-noise, named Confidence Paramerter (PC). Quantifies
    the goodness fo a merit function after the PCA-Folding has been executed.
    It isolates the peak by a ``frac`` portion of the whole array. Then it
    computes its average and tandard deviation excluding the peak region.

    Parameters
    ----------
    merit : `~numpy.ndarray`
        One dimentional array which contains the merit function, computed
        trhough the eigenvector and eigenvalues, output from the PCA-Folding.
        The merit function determines the best period found in PCA-Folding.
    idx_max : `int`
        Position of the maximum value in the merit function (or maximum in a
        plateau region, depending on ``region_finder``). This index is
        then used to compute the estimated period from the selected
        ``iteration`` (see `pypcaf.find_period` function).
    frac : `float`
        Porcentage of the region to isolate the maximum peak.

    Returns
    -------
    sigma : `float`
    mean : `float`
    cp : `float`
    """

    region = merit.size * frac

    idx_left = int(idx_max - region / 2)
    idx_right = int(idx_max + region / 2)

    if idx_left < 0:
        idx_left = 0
    if idx_right > merit.size - 1:
        idx_right = -1

    X = np.hstack((merit[0:idx_left], merit[idx_right:-1]))

    sigma = np.std(X)
    mean = np.average(X)

    cp = (merit.max() - mean) / sigma

    return sigma, mean, cp


def folding(times, dt, T, num_div):
    """
    Classical folding algorithm with waterfall (diagram) implementation. The
    ``times`` array is reshaped respect to a specific period, ``T``, and placed
    as a waterfall diagram with a number of division of ``num_div`` or ``M``.
    The number of divisions represents the number of elements in a row (and
    later it will represent the number of eigenvalues from in a waterfall PCA).

    Parameters
    ----------
    times : `~numpy.ndarray`
        One dimensional times array with certain hidden periodicity, e.g.
        pulsar period.
    dt : `float`
        times bin.
    T : `float`
        Period used to compute the waterfall matrix, :math:`N \\timess M`. ``N``
        is strictly dependent on the period and the times bin used.
    num_div : `int`
        Number of divisions made to the times array, which corresponds to the
        number of elements in a row of the waterfall matrix.

    Returns
    -------
    remainder : `~numpy.ndarray`
        Remainder of the folding, ready to use in the `~pypcaf.light_curve`
        light-curve.
    waterfall : `~numpy.ndarray`
        :math:`N \\timess M`. ``N`` matrix (two dimensional array). The
        waterfall matrix depends on the four inputs from the `~pypfac.folding`
        function. i.e. ``times``, ``dt``, ``T``, ``num_div``.
    """

    if T < dt:
        raise TypeError('Period (T) cannot be smaller than times bin (dt)')

    # Light-curve needs to have a division with no modulus
    M = num_div        # M for ease of notation
    N = round(T / dt)  # It will only select the integer value

    # Recalculating period with a (N + 1) step
    bins = np.linspace(0, T, N + 1)

    # number of samples that will be considered for each row of the waterfall
    ns = times.size // M

    # Modulus from division, it returns an element-wise remainder
    remainder = times % T

    waterfall = np.zeros((M, N), dtype=times.dtype)
    for m in range(M):
        indices = range(ns * m, ns * (m + 1))
        waterfall[m, :] = np.histogram(remainder[indices], bins=bins)[0]

    return remainder, waterfall


def folding_fast(times, dt, T, num_div):
    """
    Newer and faster folding algorithm with waterfall (diagram)
    implementation. The ``times`` array is reshaped respect to a specific
    period, ``T``, and placed as a waterfall diagram with a number of division
    of ``num_div`` or ``M``. The number of divisions represents the number of
    elements in a row (and later it will represent the number of eigenvalues
    from in a waterfall PCA).

    Parameters
    ----------
    times : `~numpy.ndarray`
        One dimensional times array with certain hidden periodicity, e.g.
        pulsar period.
    dt : `float`
        times bin.
    T : `float`
        Period used to compute the waterfall matrix, :math:`N \\timess M`. ``N``
        is strictly dependent on the period and the times bin used.
    num_div : `int`
        Number of divisions made to the times array, which corresponds to the
        number of elements in a row of the waterfall matrix.

    Returns
    -------
    remainder : `~numpy.ndarray`
        Remainder of the folding, ready to use in the `~pypcaf.light_curve`
        light-curve.
    waterfall : `~numpy.ndarray`
        :math:`N \\timess M`. ``N`` matrix (two dimensional array). The
        waterfall matrix depends on the four inputs from the `~pypfac.folding`
        function. i.e. ``times``, ``dt``, ``T``, ``num_div``.
    """

    if T < dt:
        raise TypeError('Period (T) cannot be smaller than times bin (dt)')

    # Light-curve needs to have a division with no modulus
    M = num_div        # M for ease of notation
    N = round(T / dt)  # It will only select the integer value

    # number of samples that will be considered for each row of the waterfall
    ns = times.size // M

    # Modulus from division, it returns an element-wise remainder
    remainder = times % T

    def hist(a, bins):
        return np.histogram(a=a, bins=bins)[0]

    if times.dtype == np.float64 or times.dtype == np.dtype('<f8'):
        waterfall = np.apply_along_axis(
            func1d=histogram1d,
            axis=1,
            arr=remainder[:M * ns].reshape(M, ns),
            bins=N,
            range=[0, T]
            )
    else:
        waterfall = np.apply_along_axis(
            func1d=hist,
            axis=1,
            arr=remainder[:M * ns].reshape(M, ns),
            bins=np.linspace(0, T, N + 1)
            )

    return remainder, waterfall


def light_curve(remainder, T_folding):
    return np.histogram(remainder, T_folding)[0]


def pca_signal(EVec, waterfall):
    """
    Computes the signal or K matrix after the PCA analysis.

    Returns
    -------
    K : `~numpy.ndarray`
        Signal matrix, one of the outputs from the PCA analysis. It is
        currently not been used in the `pypcaf` main core. It represents the
        transposed ``EVec_sorted`` timess the normalized waterfall matrix.
    """

    M = waterfall.shape[0]

    mean = np.mean(waterfall, axis=1).reshape(M, 1)

    # sample standard deviation
    std = np.std(waterfall, axis=1, ddof=1).reshape(M, 1)

    # Normalization waterfall matrix to mean=0 and std=1
    norm = (waterfall - mean) / std

    K = np.dot(EVec.T, norm)

    return K


def pca(waterfall):
    """
    It returns the eigenvalues and eigenvectors (PCA) from the covariance
    matrix of a normalized waterfall array of length, math:`N\\times M`.
    math:`M` represents the number of eigenvalues and eigenvectors.

    Parameters
    ----------
    waterfall : `~numpy.ndarray`
        Array that contains an Amount of ``num_div`` rows (math:`M`). The
        number of elements in a column is N (related to the number of
        iterations in the function `~pypcaf.find_period`). It also is the
        same output as mentioned in function `~pypcaf.folding`.

    Returns
    -------
    EVal_sorted : `~numpy.ndarray`
        Eigenvalues from the covariance matrix, it is always one dimensional
        array. It has also been sorted, in such a way that their values are
        from maximum to minimum. Since the waterfall matrix has been
        normalized their values are of mean equals zero and std 1.
    EVec_sorted : `~numpy.ndarray`
        Eigenvectors from the covariance matrix, it is a multidimensional
        array, and it depends on the number of dimensions of the waterfall
        array. Their values have been sorted in the same way as the ones from
        the eigenvalues, ``EVal_sorted``, viz. their correspondent eigenvalue
        and eigenvector are in the same position in each array.
    signal : `bool`
        If `True` it will compute the signal matrix.
    """

    M, N = waterfall.shape

    mean = np.mean(waterfall, axis=1).reshape(M, 1)
    std = np.std(waterfall, axis=1, ddof=1).reshape(M, 1)

    # Normalization waterfall matrix to mean=0 and std=1
    norm = (waterfall - mean) / std           # (x - <x>) / std(x)
    cov = 1 / (N - 1) * np.dot(norm, norm.T)  # Covariance matrix

    # EVec[:, i] is the eigenvector corresponding to EVal[i] eigenvalue
    EVal, EVec = np.linalg.eig(cov)

    # Sorting the eigenvalues and vectors
    sorted_positions = np.argsort(EVal.real)[::-1]
    EVal_sorted = EVal[sorted_positions]           # eigenvalues
    EVec_sorted = EVec[:, sorted_positions]        # eigenvectors or PCs

    return EVal_sorted, EVec_sorted


def find_period(
    times, dt, T_init, iteration, delta, num_div, merit_func,
    region_order
        ):
    """
    Finds the best period given an initial starting point, ``T_init``, and
    ``iteration`` and an increase step in every iteration, ``delta``.It also
    computes the merit function, derived from the eigenvalues and
    eigenvectors in `~pypcaf.pca` function.

    Parameters
    ----------
    times : `~numpy.ndarray`
        One dimensional times array with certain hidden periodicity, e.g.
        pulsar period.
    dt : `float`
        times bin.
    T_init : `float`
        Initial period to start looking for a best estimate, ``T_est``.
    iteration : `int`
        Number of iterations in which a ``delta`` increase will be applied.
        The initial period, ``T_init``, is the center point.
    delta : `float`
        Increase of the period in each iteration. The recommended order of it
        is between ``1e-7`` and ``1e-8``.
    num_div : `int`
        Number of divisions made to the times array, which corresponds to the
        number of elements in a row of the waterfall matrix.
    merit_func : `function`
        It computes the merit function from eigenvalues and scalar arrays.
        Both of them should be a one dimensional array.
    region_order : `int`
        It makes use of the `~pypcaf.flat_region_finder` to search for the
        maximum in the selected merit function. If ``region_order = 1``,
        it will compute the ordinary maximum of the merit array, i.e.
        ``np.max(merit)``. This function defines the estimated period after
        one ``iteration``.

    Returns
    -------
    T_est : `float`
        Estimated period from the iteration. The period corresponds to
        ``T_est = T_init - iteration / 2 * delta + idx_max * delta``.
    EValw : `~numpy.ndarray`
        Set of eigenvalues of math:`N` length, i.e. the number of
        iterations, ``iteration``. EValw[:, 0] represents all iteration for the
        first eigenvalue.
    Sw : `~numpy.ndarray`
        Set of scalars of math:`N` length, i.e. the number of
        iterations, ``iteration``. The scalar value is the projection of each
        for the eigenvectors into the hyper-diagonal unitary vector. In other
        words, the dot product of the (absolute value) eigenvector timess the
        same dimension unitary vector. Sw[:, 0] represents all iteration for
        the first eigenvector.
    merit : `~numpy.ndarray`
        Output from the evaluation of ``merit_func``.
    idx_max : `int`
        Position of the maximum value in the merit function (or maximum in a
        plateau region, depending on ``region_finder``). This index is
        then used to compute the estimated period from the selected
        ``iteration``.
    """

    M = num_div
    # iters = np.linspace(0, iteration, iteration, endpoint=False, dtype=int)
    # T_iteration = (iters - np.mean(iters, dtype=int)) * delta
    T_iteration = np.linspace(
        start=T_init - delta * iteration / 2,
        stop=T_init + delta * iteration / 2,
        num=iteration,
        endpoint=False,
        dtype=times.dtype
        )

    EValw = np.zeros((iteration, M), dtype=T_iteration.dtype)
    Sw = np.zeros_like(EValw)
    u = np.ones((M, 1)) / np.sqrt(M)  # hyper-diagonal unitary vector

    for t, T_iterated in enumerate(T_iteration):
        # Computing a new waterfall for every iteration
        waterfall = folding_fast(times=times, dt=dt, T=T_iterated, num_div=M)[1]
        # waterfall = folding(times=times, dt=dt, T=T_iterated, num_div=M)[1]
        EValw[t, :], EVec = pca(waterfall=waterfall)
        Sw[t, :] = np.abs((EVec * u).sum(axis=0))

    # scalar and eigenvalues from waterfall N x M matrix
    # Sw[:, 0] all iterations for the first scalar
    # EValw[:, 0] all iterations for the first eigenvalue

    # Future development decouple the period calculation
    # using an independent function that calls the merit func
    merit = merit_func(EValw=EValw, Sw=Sw)

    if region_order > 1:
        idx_max = flat_region_finder(X=merit, n=region_order)
    else:
        idx_max = merit.argmax()
    T_est = T_iteration[idx_max]

    return T_est, EValw, Sw, merit, idx_max
