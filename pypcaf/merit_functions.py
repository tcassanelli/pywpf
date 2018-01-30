#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import numpy as np

__all__ = ['merit1', 'merit2']


def merit1(EValw, Sw):

    M = Sw[0].size
    N = Sw[:, 0].size
    Sw_argmax = np.argmax(Sw, axis=1)

    EValw_corr = EValw[range(Sw_argmax.size), Sw_argmax]
    Sw_max = np.max(Sw, axis=1)

    return Sw_max * EValw_corr


def merit2(EValw, Sw):

    # Correspondent eigenvalue to the maximum selected scalar
    # V_corr = np.choose(np.argmax(Sw, axis=1), EValw.T)   # has a 32 lim!
    M = Sw[0].size
    N = Sw[:, 0].size
    Sw_argmax = np.argmax(Sw, axis=1)

    EValw_corr = EValw[range(Sw_argmax.size), Sw_argmax]

    S_avg = []  # max scalar minus its average
    for i in range(0, N):
        noise = np.sum(Sw[i]) - np.max(Sw[i])
        S_avg.append(np.max(Sw[i]) - noise / (M - 1))
    S_avg_array = np.array(S_avg)

    # (maximum scalar minus average) times the associated eigenvalue
    return S_avg_array * EValw_corr
