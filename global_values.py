# Global Values Computation
from astropy.io import ascii
import csv
import os.path
import numpy as np
from pca_analysis import flat_region_finder, rms, significant

path = 'rmr3000/method14/'
file_ = 'rmr3000_test68'

name_ = 'gv_' + file_

# Check if the sub directory is already created
if not os.path.isdir('glob_var'):
    os.makedirs('glob_var')

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

        if os.path.isfile('glob_var/' + name_ + '.csv'):
            with open('glob_var/' + name_ + '.csv', 'a') as outcsv:
                writer2 = csv.writer(outcsv)
                writer2.writerow(TO_SAVE)
        else:
            with open('glob_var/' + name_ + '.csv', 'a') as outcsv:	
                writer = csv.DictWriter(outcsv, fieldnames = titles)
                writer.writeheader()
                writer2 = csv.writer(outcsv)
                writer2.writerow(TO_SAVE)