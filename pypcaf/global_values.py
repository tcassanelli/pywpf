 # Global Values Computation
from astropy.io import ascii
import csv
import os.path
import numpy as np
from statistics import pstdev
from pca_analysis import flat_region_finder
from scipy.stats import kurtosis, skew


def significant(vector, iteration=0):
    """
    Checks if a peak inside a one dimentional array is or not significant with respect to the noise.

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

    # Select carefully this values, choose according to plot_test.py
    if iteration == 0:
        pos_i = 40
        pos_f = 60
    else:
        pos_i = 40
        pos_f = 60

    new_vector = np.hstack((vector[:pos_i], vector[pos_f:]))

    value_ = (flat_region_finder(vector)[1] - np.mean(new_vector)) / pstdev(new_vector)

    return value_, np.mean(new_vector), pstdev(new_vector)


rmr = ['rmr0', 'rmr100', 'rmr500', 'rmr1600', 'rmr2000', 'rmr2600', 'rmr3000']
test = ['71', '72', '73', '74', '75', '76', '77']

for i in range(len(rmr)):

    path = rmr[i] + '/'
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
            p2 = data['PERIOD_ITER' + str(iteration+1)][k] # Final period selected
            num_div = data['NUM_DIV'][k]
            mean = significant(merit[k])[1] # Mean without the peak
            rms = significant(merit[k])[2] # RMS without the peak

            TO_SAVE = [file_, num_div, iteration + 1, bintime, p2, mean, rms, flat_region_finder(merit[k])[1], \
            significant(merit[k])[0], np.max(merit[k]), kurtosis(merit[k]), skew(merit[k])]

            if os.path.isfile('glob_var/' + name_ + '.csv'):
                with open('glob_var/' + name_ + '.csv', 'a') as outcsv:
                    writer2 = csv.writer(outcsv)
                    writer2.writerow(TO_SAVE)
            else:
                with open('glob_var/' + name_ + '.csv', 'a') as outcsv:
                    writer = csv.DictWriter(outcsv, fieldnames=titles)
                    writer.writeheader()
                    writer2 = csv.writer(outcsv)
                    writer2.writerow(TO_SAVE)
