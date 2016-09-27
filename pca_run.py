from __future__ import division
import numpy as np
from pca_analysis import *
import csv
import sys
import os.path
import time
from progressbar import ProgressBar

# pca_run.py runs the functions from pca_analysis and stores them in .csv and .npy formats
# approx 15 min with num_div = 20

print('Starting pca_run.py')
print('Loading CSV files ...')

# Initial Parameter list
path = 'data_pulsar/'
file_name = 'rmr0'

# Generation of the data file to use
time = np.genfromtxt(path + file_name + '.csv')

method = 'method5'

# Name of the file created with CSV extension.
output_name = 'rmr0_test18'

dt = 0.004 # 4 ms
period_start = 0.089367 # s
# num_div = 20

iter1 = 100
delta1 = 1e-7 # 1e-7

# Give zero value to not compute them.
iter2 = 500 # before 500
delta2 = 4e-9 # before 4e-9 ns


titles = ('BINTIME', 'PERIOD_START', 'NUM_DIV', 'ITER1', 'DELTA1', 'ITER2', 'DELTA2', 'PERIOD_ITER1', \
	'PERIOD_ITER2', 'MAX_IND1', 'MAX_IND2', 'MAX_VAR1', 'MAX_VAR2')

# Interval to compute the num_div 
interval_to_compute = [4, 6, 8, 10, 12, 14, 16, 18, 20]

# Per = GIVEN, FFT, ITER1, ITER2
print('Starting iteration of ' + str(len(interval_to_compute)) + ' loops')
progress = ProgressBar()
print(progress)

# Check if the directory is already created
if not os.path.isdir(file_name):
	os.makedirs(file_name)
if not os.path.isdir(file_name + '/' + method):
	os.makedirs(file_name + '/' + method)

path_to_save = file_name + '/' + method + '/'

# Generation of the dictionaries that will contain the output npy files
var_files = ['V1_iter1', 'V2_iter1', 'V3_iter1', 'V1_iter2', 'V2_iter2', 'V3_iter2']
variance = {}
for i in var_files:
	if os.path.isfile(path_to_save + output_name + '_' + i + '.npy'):
		variance[i] = np.load(path_to_save + output_name + '_' + i + '.npy').tolist()
	else:
		variance[i] = []

scalar_files = ['S1_iter1', 'S2_iter1', 'S3_iter1', 'S1_iter2', 'S2_iter2', 'S3_iter2']
scalar = {}
for i in scalar_files:
	if os.path.isfile(path_to_save + output_name + '_' + i + '.npy'):
		scalar[i] = np.load(path_to_save + output_name + '_' + i + '.npy').tolist()
	else:
		scalar[i] = []

mstev_files = ['mstev_iter1', 'mstev_iter2'] # max scalar to eigenvalues
mstev = {}
for i in mstev_files:
	if os.path.isfile(path_to_save + output_name + '_' + i + '.npy'):
		mstev[i] = np.load(path_to_save + output_name + '_' + i + '.npy').tolist()
	else:
		mstev[i] = []

# Being iteration for number of rows in waterfall. Higher num_div higher noise!
for num_div in interval_to_compute:
		
	Per, [var_iter1, var_iter2], [scalar_iter1, scalar_iter2], [mstev_iter1, mstev_iter2], max_index = find_period(time, period_start, \
		dt, num_div, iter1, delta1, iter2, delta2, noisy_signal=True)

	INPUT = [dt, period_start, num_div, iter1, delta1, iter2, delta2]
	OUTPUT = [Per[2], Per[3], max_index[0], max_index[1], np.max(var_iter1[0]), np.max(var_iter2[0])]

	TO_SAVE = INPUT + OUTPUT

	if os.path.isfile(path_to_save + output_name + '.csv'):
		with open(path_to_save + output_name + '.csv', 'a') as outcsv:
			writer2 = csv.writer(outcsv)
			writer2.writerow(TO_SAVE)
	else:
		with open(path_to_save + output_name + '.csv', 'a') as outcsv:	
			writer = csv.DictWriter(outcsv, fieldnames = titles)
			writer.writeheader()
			writer2 = csv.writer(outcsv)
			writer2.writerow(TO_SAVE)

	# iteration in half of the scalar and variables array
	for i in range(0, 3): 
		variance[var_files[:3][i]].append(var_iter1[i])
		variance[var_files[-3:][i]].append(var_iter2[i])
		scalar[scalar_files[:3][i]].append(scalar_iter1[i])
		scalar[scalar_files[-3:][i]].append(scalar_iter2[i])

	mstev['mstev_iter1'].append(mstev_iter1)
	mstev['mstev_iter2'].append(mstev_iter2)

	# Finish with progress bar update
	print(progress + 10 / len(interval_to_compute))

for j in var_files:
	np.save(path_to_save + output_name + '_' + j + '.npy', np.array(variance[j]))
for j in scalar_files:
	np.save(path_to_save + output_name + '_' + j + '.npy', np.array(scalar[j]))
for j in mstev_files:
	np.save(path_to_save + output_name + '_' + j + '.npy', np.array(mstev[j]))

print('It is done!')