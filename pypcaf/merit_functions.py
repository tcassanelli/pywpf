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

    sorted_positions = np.argmax(Sw, axis=1)
    Sw_max = Sw[range(sorted_positions.size), sorted_positions]
    EValw_corr = EValw[range(sorted_positions.size), sorted_positions]
    merit = Sw_max * EValw_corr

    return merit


def merit2(EValw, Sw):
    """
    Merit function calculation. Computes the merit function by selecting the
    first scalar, ``Sw``. The first scalar is constructed in
    `pypcaf.find_period` and `pypcaf.pca` functions. First the eigenvalues are
    sorted (maximum in the zero position) and then the eigenvector are ordered
    according to the same relocation done in the sorting of the eigenvalues.
    Later the hyper-diagonal is multiplied by the new sort of the
    eigenvectors. Finally the first scalar is multiplied by the first
    eigenvalues.

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
        The merit function is calculated by selecting only the first position
        or first scalar and then multiplied it by its correspondent eigenvalue
        which coincides with the first eigenvalue as well.
    """

    merit = Sw[:, 0] * EValw[:, 0]

    return merit


def merit3(EValw, Sw):

    weight = Sw[:, 0] - np.mean(Sw[:, 1:], axis=1)
    weight[weight < 0] = 0

    return EValw[:, 0] * weight
