from astropy.io import ascii
import datetime
import os
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from matplotlib.backends.backend_pdf import PdfPages

# Inputs to change the file names, do not forget!
path = 'rmr500/method14/'
test_num = '69'
method = 'Method 14'
rmr = 'rmr500'
version = 'v1'

# To plot only 3 scalars and eigevalues.
S0_i1 = np.load(path + rmr + '_test' + test_num + '_S0_iter1.npy')
S1_i1 = np.load(path + rmr + '_test' + test_num + '_S1_iter1.npy')

S0_i2 = np.load(path + rmr + '_test' + test_num + '_S0_iter2.npy')
S1_i2 = np.load(path + rmr + '_test' + test_num + '_S1_iter2.npy')

V0_i1 = np.load(path + rmr + '_test' + test_num + '_V0_iter1.npy')
V1_i1 = np.load(path + rmr + '_test' + test_num + '_V1_iter1.npy')

V0_i2 = np.load(path + rmr + '_test' + test_num + '_V0_iter2.npy')
V1_i2 = np.load(path + rmr + '_test' + test_num + '_V1_iter2.npy')

mstev1 = np.load(path + rmr + '_test' + test_num + '_mstev_iter1.npy')
mstev2 = np.load(path + rmr + '_test' + test_num + '_mstev_iter2.npy')

data = ascii.read(path + rmr + '_test' + test_num + '.csv')


with PdfPages('plots/test' + str(test_num) + '_' + str(version) + '.pdf') as pdf:

	for k in range(0, len(data['BINTIME'])):

		# Normalization for mstev1
		# mstev1_norm = (mstev1[k] - np.min(mstev1[k])) / (np.max(mstev1[k]) - np.min(mstev1[k]))
		# mstev2_norm = (mstev2[k] - np.min(mstev1[k])) / (np.max(mstev1[k]) - np.min(mstev1[k]))

		x_axis_1 = np.linspace(data['PERIOD_START'][k] - data['DELTA1'][k] * data['ITER1'][k] / 2, \
			data['PERIOD_START'][k] + data['DELTA1'][k] * (data['ITER1'][k]) / 2, data['ITER1'][k], endpoint=False)
		x_axis_2 = np.linspace(data['PERIOD_ITER1'][k] - data['DELTA2'][k] * data['ITER2'][k] / 2, \
			data['PERIOD_ITER1'][k] + data['DELTA2'][k] * (data['ITER2'][k]) / 2, data['ITER2'][k], endpoint=False)

		fig = plt.figure(figsize=(20,10))
		fig.subplots_adjust(hspace=.26, left=0.04, bottom=0.07, right=0.97, top=0.94, wspace=0.2)

		ax1 = fig.add_subplot(211)
		ax1.plot(x_axis_1, S0_i1[k], 'g+-', linewidth=1, label='First scalar')
		ax1.plot(x_axis_1, S1_i1[k], 'g+-', linewidth=0.8, label='Second scalar', alpha=0.5)
		# ax1.plot(x_axis_1, S3_0[k], 'g+-', linewidth=0.8, label='Third scalar', alpha=0.3)
		ax1.plot(x_axis_1, V0_i1[k], 'b+-', linewidth=1, label='First eigenvalue')
		ax1.plot(x_axis_1, V1_i1[k], 'b+-', linewidth=0.8, label='Second eigenvalue', alpha=0.5)
		# ax1.plot(x_axis_1, V3_0[k], 'b+-', linewidth=0.8, label='Third eigenvalue', alpha=0.3)
		ax1.plot(x_axis_1, mstev1[k], 'r.-', linewidth=2.5, label='Merit function 1')
		ax1.axvline(x=0.08936715, color='y', linewidth=1, label='Folding Period = 0.08936715')
		ax1.axvline(x=data['PERIOD_ITER1'][k], color='m', linewidth=1, label='Best Period = ' + \
			str(data['PERIOD_ITER1'][k]))
		ax1.grid()
		ax1.legend(loc='upper left', fontsize=10)
		ax1.set_xlabel('Time [s]. Iterations = ' + str(data['ITER1'][k]) + ', delta = ' + str(data['DELTA1'][k]))
		ax1.set_ylabel('Scalar, eigvalue and merit amplitude')
		ax1.set_title(str(method) + '. ' + str(rmr) + ' Iteration 1. Waterfall rows = ' + \
			str(data['NUM_DIV'][k]) + '. dt = ' + str(data['BINTIME'][k]) + ' s.')
		ax1.set_xlim([data['PERIOD_START'][k] - data['DELTA1'][k] * (data['ITER1'][k] + 1) / 2, \
			data['PERIOD_START'][k] + data['DELTA1'][k] * (data['ITER1'][k] - 1)/2])

		ax2 = fig.add_subplot(212)
		ax2.plot(x_axis_2, S0_i2[k], 'g+-', linewidth=1, label='First scalar')
		ax2.plot(x_axis_2, S1_i2[k], 'g+-', linewidth=0.8, label='Second scalar', alpha=0.5)
		# ax2.plot(x_axis_2, S3_1[k], 'g+-', linewidth=0.8, label='Third scalar', alpha=0.3)
		ax2.plot(x_axis_2, V0_i2[k], 'b+-', linewidth=1, label='First eigenvalue')
		ax2.plot(x_axis_2, V1_i2[k], 'b+-', linewidth=0.8, label='Second eigenvalue', alpha=0.5)
		# ax2.plot(x_axis_2, V3_1[k], 'b+-', linewidth=0.8, label='Third eigenvalue', alpha=0.3)
		ax2.plot(x_axis_2, mstev2[k], 'r.-', linewidth=2.5, label='Merit function 2')
		ax2.plot(x_axis_1, mstev1[k], 'r.-', linewidth=1, label='Merit function 1')
		ax2.axvline(x=0.08936715, color='y', linewidth=1, label='Folding Period = 0.08936715')
		ax2.axvline(x=data['PERIOD_ITER2'][k], color='m', linewidth=1, label='Best Period = ' + \
			str(data['PERIOD_ITER2'][k]))
		ax2.grid()
		ax2.legend(loc='upper left', fontsize=10)
		ax2.set_xlabel('Time [s]. Iterations = ' + str(data['ITER2'][k]) + ', delta = ' + str(data['DELTA2'][k]))
		ax2.set_ylabel('Scalar, eigvalue and merit amplitude')
		ax2.set_title(str(method) + '. ' + str(rmr) + ' Iteration 2. Waterfall rows = ' + \
			str(data['NUM_DIV'][k]) + '. dt = ' + str(data['BINTIME'][k]) + ' s.')
		ax2.set_xlim([data['PERIOD_ITER1'][k] - data['DELTA2'][k] * (data['ITER2'][k] + 1) / 2, \
			data['PERIOD_ITER1'][k] + data['DELTA2'][k] * (data['ITER2'][k] - 1)/2])

		pdf.savefig()  # saves the current figure into a pdf page

	plt.close()

	d = pdf.infodict()
	d['Title'] = 'PCA Analysis'
	d['Author'] = u'Tomas Cassanelli'
	d['Keywords'] = 'PCA'
	d['ModDate'] = datetime.datetime.today()

# To open the already stored file
os.system('open ' + 'plots/test' + str(test_num) + '_' + str(version) + '.pdf')
