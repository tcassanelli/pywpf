#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Tomas Cassanelli
import os
import numpy as np
from astropy.io import ascii


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


if __name__ == '__main__':

    import glob

    folders = glob.glob('/Users/tomascassanelli/PCAF/data_pulsar/pypcaf_out/*')

    for i in folders:
        npz_to_csv(path_pypcaf_out=i)
