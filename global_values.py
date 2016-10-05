# Global Values Computation
from astropy.io import ascii
import csv
import sys
import os.path
import time
import numpy as np
from scipy import stats
from pca_analysis import flat_region_finder

def rms(x):
    return np.sqrt(x.dot(x)/x.size)

def significant(vector):
    if flat_region_finder(vector)[1] / (np.mean(vector) + rms(vector)) > 1:
        to_print = 'YES'
    else:
        to_print = 'NO'
    return to_print

path = 'rmr3000/method14/'
file_ = 'rmr3000_test68'

name_ = 'gv_' + file_

data = ascii.read(path + file_ + '.csv')
mstev = [np.load(path + file_ + '_mstev_iter1.npy'), np.load(path + file_ + '_mstev_iter2.npy')]

# PEAK_SIGNIFICANT = (peak/(avg+rms) > 1
titles = ['FILE', 'NUM_DIV','ITER','BINTIME', 'PERIOD_ITER', 'MEAN', 'RMS', 'PEAK_SELECTED', 'PEAK_SIGNIFICANT', 'PEAK_REAL']

N = len(data)

for iteration in range(0, 2):
    merit = mstev[iteration]
    for k in range(0, N):

        bintime = data['BINTIME'][k]
        p2 = data['PERIOD_ITER' + str(iteration+1)][k]
        num_div = data['NUM_DIV'][k]

        TO_SAVE = [file_, num_div, iteration + 1, bintime, p2, np.mean(merit[k]), rms(merit[k]), \
        flat_region_finder(merit[k])[1], significant(merit[k]), np.max(merit[k])]

        if os.path.isfile(name_ + '.csv'):
            with open(name_ + '.csv', 'a') as outcsv:
                writer2 = csv.writer(outcsv)
                writer2.writerow(TO_SAVE)
        else:
            with open(name_ + '.csv', 'a') as outcsv:	
                writer = csv.DictWriter(outcsv, fieldnames = titles)
                writer.writeheader()
                writer2 = csv.writer(outcsv)
                writer2.writerow(TO_SAVE)