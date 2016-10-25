from __future__ import division
import numpy as np
from pca_analysis import *
import csv
import os.path
from progressbar import ProgressBar

# SECONDARY PROGRAM FOR PCA FOLDING
# AUTHOR: TOMAS CASSANELLI

# pca_run.py runs the functions from pca_analysis and stores them in .csv and .npy formats
# The purpose is to find the Maximum Scalar To EigenValue (mstev) or merit function to find the periods.
# Beware of the computing time!
# Run again the program will not overwrite files, it will append them!

#################### INITIAL PARAMETERS ####################

file_name = 'FILE_NAME'  # Contains the time array
output_name = file_name + '_output1'

# Number of rows in waterfall for second interation
interval_to_compute = [3, 4, 5]

dt = 0.002793  # 4 ms, 0.002793 s
period_start = 0.089367  # Initial period, usualy well known

iter1 = 100
delta1 = 1e-7

iter2 = 400
delta2 = 10e-9

noise_signal = True
# Change to Flase to search period starting with FFT
# Use it only when signal is very clean (i.e. Crab Pulsar)

############################################################

print('Starting pca_run.py')
print('Loading CSV files ...')

# Generation of the data file to use
time = np.genfromtxt('data_pulsar/' + file_name + '.csv')

titles = ('BINTIME', 'PERIOD_START', 'NUM_DIV', 'ITER1', 'DELTA1', 'PERIOD_ITER1', 'ITER2', 'DELTA2', 'PERIOD_ITER2')

# Per = GIVEN, FFT, ITER1, ITER2
print('Starting iteration of ' + str(len(interval_to_compute)) + ' loops')
progress = ProgressBar()
print(progress)

# Check if the directory is already created
if not os.path.isdir(file_name):
    os.makedirs(file_name)

# Generation of the dictionaries that will contain the output npy files
var_files_iter1 = []
var_files_iter2 = []
for i in range(0, interval_to_compute[0]):
    var_files_iter1.append('V' + str(i) + '_iter1')
    var_files_iter2.append('V' + str(i) + '_iter2')
var_files = var_files_iter1 + var_files_iter2

scalar_files_iter1 = []
scalar_files_iter2 = []
for i in range(0, interval_to_compute[0]):
    scalar_files_iter1.append('S' + str(i) + '_iter1')
    scalar_files_iter2.append('S' + str(i) + '_iter2')
scalar_files = scalar_files_iter1 + scalar_files_iter2

variance = {}
for i in var_files:
    if os.path.isfile(file_name + '/' + output_name + '_' + i + '.npy'):
        variance[i] = np.load(file_name + '/' + output_name + '_' + i + '.npy').tolist()
    else:
        variance[i] = []

scalar = {}
for i in scalar_files:
    if os.path.isfile(file_name + '/' + output_name + '_' + i + '.npy'):
        scalar[i] = np.load(file_name + '/' + output_name + '_' + i + '.npy').tolist()
    else:
        scalar[i] = []

mstev_files = ['mstev_iter1', 'mstev_iter2']  # max scalar to eigenvalues (merit function)
mstev = {}
for i in mstev_files:
    if os.path.isfile(file_name + '/' + output_name + '_' + i + '.npy'):
        mstev[i] = np.load(file_name + '/' + output_name + '_' + i + '.npy').tolist()
    else:
        mstev[i] = []

# Being iteration for number of rows in waterfall. Higher num_div higher noise!
for num_div in interval_to_compute:

    # Loop in the main function from pca_analysis.py.
    Per, [var_iter1, var_iter2], [scalar_iter1, scalar_iter2], [mstev_iter1, mstev_iter2], _ \
     = find_period(time, period_start, dt, num_div, iter1, delta1, iter2, delta2, noisy_signal=noise_signal)

    # Per = [period, period_start1, period_final1, period_final2]

    TO_SAVE = [dt, Per[1], num_div, iter1, delta1, Per[2], iter2, delta2, Per[3]]

    if os.path.isfile(file_name + '/' + output_name + '.csv'):
        with open(file_name + '/' + output_name + '.csv', 'a') as outcsv:
            writer2 = csv.writer(outcsv)
            writer2.writerow(TO_SAVE)
    else:
        with open(file_name + '/' + output_name + '.csv', 'a') as outcsv:
            writer = csv.DictWriter(outcsv, fieldnames = titles)
            writer.writeheader()
            writer2 = csv.writer(outcsv)
            writer2.writerow(TO_SAVE)

    # iteration in half of the scalar and variables array
    for i in range(0, interval_to_compute[0]):
        variance[var_files_iter1[i]].append(var_iter1[:, i])
        variance[var_files_iter2[i]].append(var_iter2[:, i])
        scalar[scalar_files_iter1[i]].append(scalar_iter1[:, i])
        scalar[scalar_files_iter2[i]].append(scalar_iter2[:, i])

    mstev['mstev_iter1'].append(mstev_iter1)
    mstev['mstev_iter2'].append(mstev_iter2)

    # Finish with progress bar update
    print(progress + 10 / len(interval_to_compute))

# Store files with a .npy extension, to open use np.load()
for j in var_files:
    np.save(file_name + '/' + output_name + '_' + j + '.npy', np.array(variance[j]))
for j in scalar_files:
    np.save(file_name + '/' + output_name + '_' + j + '.npy', np.array(scalar[j]))
for j in mstev_files:
    np.save(file_name + '/' + output_name + '_' + j + '.npy', np.array(mstev[j]))

print('It is done!')
