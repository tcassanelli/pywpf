#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import numpy as np

__all__ = ['merit1', 'merit2']


def merit1(EValw, Sw):
    """
    Merit function calculation. Computes the merit function out of the
    eigenvalues and scalar product after the waterfall diagram has been
    calculated.

    Parameters
    ----------
    EValw : `~numpy.ndarray`
        Set of eigenvalues of math:`N` length, i.e. the number of
        iterations, ``iteration``. EValw[:, 0] represents all iteration for the
        first eigenvalue.
    Sw : `~numpy.ndarray`
        Set of scalars of math:`N` length, i.e. the number of
        iterations, ``iteration``. The scalar value is the projection of each
        for the eigenvectors into the hyper-diagonal unitary vector. In other
        words, the dot product of the (absolute value) eigenvector times the
        same dimension unitary vector. Sw[:, 0] represents all iteration for
        the first eigenvector.

    Returns
    -------
    merit : `~numpy.ndarray`
        The merit function is calculated by selecting the maximum scalar over
        each iteration and then it is multiplied by its correspondent
        eigenvalue.
    """

    Sw_argmax = np.argmax(Sw, axis=1)
    EValw_corr = EValw[range(Sw_argmax.size), Sw_argmax]
    Sw_max = np.max(Sw, axis=1)
    merit = Sw_max * EValw_corr

    return merit


def merit2(EValw, Sw):
    """
    Merit function calculation. Computes the merit function out of the
    eigenvalues and scalar product after the waterfall diagram has been
    calculated.

    Parameters
    ----------
    EValw : `~numpy.ndarray`
        Set of eigenvalues of math:`N` length, i.e. the number of
        iterations, ``iteration``. EValw[:, 0] represents all iteration for the
        first eigenvalue.
    Sw : `~numpy.ndarray`
        Set of scalars of math:`N` length, i.e. the number of
        iterations, ``iteration``. The scalar value is the projection of each
        for the eigenvectors into the hyper-diagonal unitary vector. In other
        words, the dot product of the (absolute value) eigenvector times the
        same dimension unitary vector. Sw[:, 0] represents all iteration for
        the first eigenvector.

    Returns
    -------
    merit : `~numpy.ndarray`
        The merit function is calculated by selecting the maximum scalar over
        each iteration and then finding its correspondent eigenvalue. Later a
        modified scalar is found by adding a noise term to each iteration in
        the scalar. Finally this modified scalar is multiplied by the
        correspondent eigenvalues (same oreder as the maximum scalar).
    """

    M = Sw[0].size
    N = Sw[:, 0].size
    Sw_argmax = np.argmax(Sw, axis=1)

    EValw_corr = EValw[range(Sw_argmax.size), Sw_argmax]

    S_mod = []  # max scalar minus its average
    for i in range(0, N):
        noise = np.sum(Sw[i]) - np.max(Sw[i])
        S_mod.append(np.max(Sw[i]) - noise / (M - 1))

    merit = np.array(S_mod) * EValw_corr

    return merit
