#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
import numpy as np
from astropy.io import ascii
import glob


def npz_to_csv(path_pypcaf_out):

    info = ascii.read(os.path.join(path_pypcaf_out, 'info.dat'))
    num_div = info['num_div']

    for m in num_div:

        path_to_data = os.path.join(
            path_pypcaf_out, 'M{}.npz'.format(m)
            )

        pcaf_out_M = np.load(path_to_data)

        for key, value in pcaf_out_M.items():
            path_to_save = os.path.join(
                path_pypcaf_out, 'M{}_'.format(m) + key + '.csv'
                )
            np.savetxt(path_to_save, value)


def npy_to_csv(path_data):
    '''
    The default file extention for this function is .npy
    '''

    if path_data.endswith('.npy'):

        path, name = os.path.split(path_data)
        path_to_save = os.path.join(path, name[:-4] + '.csv')
        np.savetxt(path_to_save, np.load(path_data))

    else:
        'Wrong format, only .npy'


def double_check(path_data):
    '''
    Double check for the .npy data converted to .csv
    '''

    npy = np.sort(glob.glob(os.path.join(path_data, '*.npy')))
    csv = np.sort(glob.glob(os.path.join(path_data, '*.csv')))

    for i, j in zip(npy, csv):
        if not np.array_equal(np.load(i), np.genfromtxt(j)):
            print('the array {} and array {} are not equal')
    
    print('all arrays are fine!')


if __name__ == '__main__':

    pypcaf_out = glob.glob('../../data_pulsar2/pypcaf_out/*')

    for i in pypcaf_out:
        npz_to_csv(path_pypcaf_out=i)

    # data = glob.glob('../../data_pulsar2/*.npy')

    # for j in data:
    #     npy_to_csv(path_data=j)


    # double_check('../../data_pulsar2/')
