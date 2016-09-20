from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from pca_analysis import *
import scipy.io as sio
import csv
import sys
import os.path
import time
from progressbar import ProgressBar

path = '/Users/tomascassanelli/Dropbox/PCA.A&A/Verroi/MATLAB/'

file_name = 'rmr0'
matlab_file = sio.loadmat(path + file_name + '.mat')
time = matlab_file['time']

method = 'method1'

# Name of the file created with CSV extension.
output_name = 'rmr0_2iter'

dt = 0.004 # 4 ms
period_start = 0.089367 # s
# num_div = 20

iter1 = 100
delta1 = 1e-7 # 1e-7

# Give zero value to not compute them.
iter2 = 1000 # before 1000
delta2 = 1e-9 # before 1e-9 ns

# Testing data
# bindata, frequency_fft = pre_analysis(time, dt, period_start, True)
# lc, waterfall = new_fold(time, dt, period_start, num_div, True)
# V, PC, cov, norm, signals = fast_pca(waterfall, True)

titles = ('BINTIME', 'PERIOD_START', 'NUM_DIV', 'ITER1', 'DELTA1', 'ITER2', 'DELTA2', 'PERIOD_ITER1', \
	'PERIOD_ITER2', 'MAX_VAR1', 'MAX_VAR2')

# Per = GIVEN, FFT, ITER1, ITER2

progress = ProgressBar()
print('Starting pca_run.py')
print(progress)

# Interval to compute the num_div 
interval_to_compute = np.arange(15, 20, dtype=int)

# Check if the directory is already created
if not os.path.isdir(file_name):
	os.makedirs(file_name)
if not os.path.isdir(file_name + '/' + method):
	os.makedirs(file_name + '/' + method)

path_to_save = file_name + '/' + method + '/'

# if os.path.isfile(path_to_save + output_name + '_V.npy'):
# 	Eigenvalues_all1 = np.load(path_to_save + output_name + '_V.npy').tolist()
# else: 
# 	Eigenvalues_all1 = []

for num_div in interval_to_compute:
		
	Per, Var, Eigenvalues, Eigenvectors = find_period(time, period_start, dt, num_div, \
		iter1, delta1, iter2, delta2, noisy_signal=True)

	INPUT = [dt, period_start, num_div, iter1, delta1, iter2, delta2]
	OUTPUT = [Per[2], Per[3], np.max(Var[0]), np.max(Var[1])]

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

	# Eigenvalues_all1.append(Eigenvalues[0])

	# Finish with progress bar update
	print(progress + 10 / len(interval_to_compute))

# Set it to True if you want to save the files
if False:
	np.save(path_to_save + output_name + '_V', np.array(Eigenvalues_all1))


print('It is done!')