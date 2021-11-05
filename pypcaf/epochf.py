#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import numpy as np

__all__ = ['epoch_folding']


def epoch_folding(times, T_init, dt, delta, iteration):
    """
    Simple chi-squared epoch folding algorithm.
    """

    if T_init < dt:
        raise TypeError('Period cannot be smaller than time bin (dt)')

    N = round(T_init / dt)
    bins = np.linspace(0, 1, N + 1)

    T_iteration = np.linspace(
        start=T_init - delta * iteration / 2,
        stop=T_init + delta * iteration / 2,
        num=iteration,
        endpoint=False,
        dtype=times.dtype
        )

    freqs = 1 / T_iteration.copy()
    t0 = (times[-1] + times[0]) / 2

    chisq_min = 0.
    chisq = np.zeros(freqs.shape, dtype=times.dtype)
    for k, freq in enumerate(freqs):
        phases = freq * (times - t0) % 1.0
        fold, edges = np.histogram(phases, bins=bins)

        chisq_iter = np.sum((fold - fold.mean()) ** 2 / fold)
        if chisq_iter > chisq_min:
            chisq[k] = chisq_iter

    T_est = T_iteration[chisq.argmax()]

    return chisq, T_iteration, T_est
