from astropy.io import ascii
import os
import numpy as np

# Translating .npy files to csv

# Inputs to change the file names, do not forget!
path = 'rmr3000/rmr3000_test78'
save = 'rmr3000_test78'

S0_i1 = np.load(path + '_S0_iter1.npy')
S1_i1 = np.load(path + '_S1_iter1.npy')

S0_i2 = np.load(path + '_S0_iter2.npy')
S1_i2 = np.load(path + '_S1_iter2.npy')

V0_i1 = np.load(path + '_V0_iter1.npy')
V1_i1 = np.load(path + '_V1_iter1.npy')

V0_i2 = np.load(path + '_V0_iter2.npy')
V1_i2 = np.load(path + '_V1_iter2.npy')

mstev1 = np.load(path + '_mstev_iter1.npy')
mstev2 = np.load(path + '_mstev_iter2.npy')

np.savetxt('csv_data/' + save + '_S0_iter1.csv', S0_i1)
np.savetxt('csv_data/' + save + '_S1_iter1.csv', S1_i1)

np.savetxt('csv_data/' + save + '_S0_iter2.csv', S0_i2)
np.savetxt('csv_data/' + save + '_S1_iter2.csv', S1_i2)

np.savetxt('csv_data/' + save + '_V0_iter1.csv', V0_i1)
np.savetxt('csv_data/' + save + '_V1_iter1.csv', V1_i1)

np.savetxt('csv_data/' + save + '_V0_iter2.csv', V0_i2)
np.savetxt('csv_data/' + save + '_V1_iter2.csv', V1_i2)

np.savetxt('csv_data/' + save + '_mstev_iter1.csv', mstev1)
np.savetxt('csv_data/' + save + '_mstev_iter2.csv', mstev2)


# Calculation the time axis
data = ascii.read(path + '.csv')
length = len(data)


time_axis_iter1 = []
time_axis_iter2 = []
for k in range(length):

    x_axis_1 = np.linspace(data['PERIOD_START'][k] - data['DELTA1'][k] * data['ITER1'][k] / 2, \
        data['PERIOD_START'][k] + data['DELTA1'][k] * (data['ITER1'][k]) / 2, data['ITER1'][k], endpoint=False)
    x_axis_2 = np.linspace(data['PERIOD_ITER1'][k] - data['DELTA2'][k] * data['ITER2'][k] / 2, \
        data['PERIOD_ITER1'][k] + data['DELTA2'][k] * (data['ITER2'][k]) / 2, data['ITER2'][k], endpoint=False)

    time_axis_iter1.append(x_axis_1)
    time_axis_iter2.append(x_axis_2)


np.savetxt('csv_data/' + save + '_time_iter1.csv', time_axis_iter1)
np.savetxt('csv_data/' + save + '_time_iter2.csv', time_axis_iter2)
