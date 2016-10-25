from __future__ import division
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.io import ascii


# Function to desribe the noisy output as a Gaussian

def gauss(x, amp, mu, sigma):
    """
    Evaluate Gaussian.

    Parameters
    ----------
    A : float
        Amplitude.
    mu : float
        Mean.
    std : float
        Standard deviation.

    """
    return amp * np.exp(-(x - mu)**2 / (2 * sigma**2))


def gauss_fit(x, y):

    # Typical interval where the peak is located
    y_ = y[130:275]
    x_ = x[130:275]

    # Working with small numbers in the curve_fit package returns problems
    coeff, var_matrix = curve_fit(gauss, np.linspace(0, 100, len(y_)), y_)

    return x_, gauss(np.linspace(0, 100, len(y_)), *coeff)


if __name__ == "__main__":

    test_num = '62'
    rmr = 'rmr0'
    version = 'v3'

    path = rmr + '/method14/'
    mstev2 = np.load(path + rmr + '_test' + test_num + '_mstev_iter2.npy')[0]
    data = ascii.read(path + rmr + '_test' + test_num + '.csv')


    k = 0
    x_axis_1 = np.linspace(data['PERIOD_START'][k] - data['DELTA1'][k] * data['ITER1'][k] / 2, \
        data['PERIOD_START'][k] + data['DELTA1'][k] * (data['ITER1'][k] - 1) / 2, data['ITER1'][k])
    x_axis_2 = np.linspace(data['PERIOD_ITER1'][k] - data['DELTA2'][k] * data['ITER2'][k] / 2, \
        data['PERIOD_ITER1'][k] + data['DELTA2'][k] * (data['ITER2'][k] - 1) / 2, data['ITER2'][k])

    y = mstev2
    x = x_axis_2

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()
    ax.plot(x, y, 'bo')
    ax.plot(gauss_fit(x, y)[0], gauss_fit(x, y)[1], linewidth=2, c='r')
    plt.show()
