#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import numpy as np

__all__ = ['chi_squared_folding']


def chi_squared_folding(times, T_init, dt, delta, iteration):
    """
    Simple chi-squared epoch folding algorithm.
    """

    if T_init < dt:
        raise TypeError('Period cannot be smaller than time bin (dt)')

    N = round(T_init / dt)
    bins = np.linspace(0, 1, N + 1)

    Tlow = T_init - iteration / 2 * delta
    Thigh = T_init + iteration / 2 * delta

    if Tlow < 0.:
        raise TypeError('iteration // 2 * dt is larger than T_init')

    freqs = np.linspace(1 / Tlow, 1 / Thigh, iteration, endpoint=False)
    t0 = (times[-1] + times[0]) / 2

    chisq_min = 0.
    chisq = np.zeros(freqs.shape, dtype=times.dtype)
    for k, freq in enumerate(freqs):
        phases = freq * (times - t0) % 1.0
        fold, edges = np.histogram(phases, bins=bins)

        chisq_iter = np.sum((fold - fold.mean()) ** 2 / fold)
        if chisq_iter > chisq_min:
            chisq[k] = chisq_iter

    T_list = 1 / freqs
    T_opt = T_list[chisq.argmax()]

    return chisq, T_list, T_opt
