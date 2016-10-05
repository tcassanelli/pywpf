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

# Initial Parameter list
file_name = 'rmr0'
method = 'method14'
test = '62'

# file_names = ['rmr500', 'rmr1600', 'rmr2000', 'rmr2600', 'rmr3000']
# method = 'method14'
# tests = ['64', '65', '66', '67', '68']

# for x in range(0, 5):
#     file_name = file_names[x]
#     test = tests[x]

# Interval to compute the num_div 
# Number of eigenvalues to save is the minimum interval_to_compute[0]
interval_to_compute = [7, 8, 9, 10]

dt = 0.002793 # 4 ms, 0.002793 s
period_start = 0.089367 # s

iter1 = 100
delta1 = 1e-7 # 1e-7

# Give zero value to not compute them.
iter2 = 400 # before 500, 400
delta2 = 10e-9 # before 4e-9 ns, 10e-9 s


print('Starting pca_run.py')
print('Loading CSV files ...')

# Generation of the data file to use
path = 'data_pulsar/'
time = np.genfromtxt(path + file_name + '.csv')

# Name of the file created with CSV extension.
output_name = file_name + '_test' + test
path_to_save = file_name + '/' + method + '/'

titles = ('BINTIME', 'PERIOD_START', 'NUM_DIV', 'ITER1', 'DELTA1', 'PERIOD_ITER1', 'ITER2', 'DELTA2', 'PERIOD_ITER2')

# Per = GIVEN, FFT, ITER1, ITER2
print('Starting iteration of ' + str(len(interval_to_compute)) + ' loops')
progress = ProgressBar()
print(progress)

# Check if the directory is already created
if not os.path.isdir(file_name):
	os.makedirs(file_name)
if not os.path.isdir(file_name + '/' + method):
	os.makedirs(file_name + '/' + method)

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
	if os.path.isfile(path_to_save + output_name + '_' + i + '.npy'):
		variance[i] = np.load(path_to_save + output_name + '_' + i + '.npy').tolist()
	else:
		variance[i] = []

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
		
	Per, [var_iter1, var_iter2], [scalar_iter1, scalar_iter2], [mstev_iter1, mstev_iter2], max_index \
	 = find_period(time, period_start, dt, num_div, iter1, delta1, iter2, delta2, noisy_signal=True)

	# Per = [period, period_start1, period_final1, period_final2]

	TO_SAVE = [dt, Per[1], num_div, iter1, delta1, Per[2], iter2, delta2, Per[3]]

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
	for i in range(0, interval_to_compute[0]): 
		variance[var_files_iter1[i]].append(var_iter1[:, i])
		variance[var_files_iter2[i]].append(var_iter2[:, i])
		scalar[scalar_files_iter1[i]].append(scalar_iter1[:, i])
		scalar[scalar_files_iter2[i]].append(scalar_iter2[:, i])

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