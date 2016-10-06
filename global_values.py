# Global Values Computation
from astropy.io import ascii
import csv
import os.path
import numpy as np
from statistics import pstdev
from pca_analysis import flat_region_finder
from scipy.stats import kurtosis, skew


def significant(vector):
    """
    Checks if a peak inside a one dimentional array is or not significant.

    Parameters
    ----------
    vector : numpy.ndarray
    
    Returns
    -------
    m : str 
        I will return 'YES' if what was stated is true, else 'NO'.
    value_ : float
        It is the ratio between the max peak and the rms + mean
    """
    value_ = (flat_region_finder(vector)[1] - np.mean(vector)) / pstdev(vector)

    if value_ > 1:
        to_print = 'YES'
    else:
        to_print = 'NO'
    return to_print, value_

rmr = ['rmr0', 'rmr100', 'rmr500', 'rmr1600', 'rmr2000', 'rmr2600', 'rmr3000']
test = ['62', '63', '64', '65', '66', '67', '68']

for i in range(len(rmr)):

    path = rmr[i] + '/method14/'
    file_ = rmr[i] + '_test' + test[i] 

    name_ = 'gv_' + file_

    # Check if the sub directory is already created
    if not os.path.isdir('glob_var'):
        os.makedirs('glob_var')

    data = ascii.read(path + file_ + '.csv')
    mstev = [np.load(path + file_ + '_mstev_iter1.npy'), np.load(path + file_ + '_mstev_iter2.npy')]

    # PEAK_SIGNIFICANT = (peak/(avg+rms) > 1
    titles = ['FILE', 'NUM_DIV','ITER','BINTIME', 'PERIOD_ITER', 'MEAN', 'RMS','PEAK_SELECTED', 'PEAK_SIGNIFICANT', \
    'PEAK_REAL', 'KURTOSIS', 'SKEWNESS']

    N = len(data)

    for iteration in range(0, 2):
        merit = mstev[iteration]
        for k in range(0, N):

            bintime = data['BINTIME'][k]
            p2 = data['PERIOD_ITER' + str(iteration+1)][k]
            num_div = data['NUM_DIV'][k]

            TO_SAVE = [file_, num_div, iteration + 1, bintime, p2, np.mean(merit[k]), pstdev(merit[k]),\
            flat_region_finder(merit[k])[1], significant(merit[k])[1], np.max(merit[k]), kurtosis(merit[k]), skew(merit[k])]

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