from astropy.io import ascii
import datetime
import os
import numpy as np
from noise_to_gauss import gauss, gauss_fit
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from progressbar import ProgressBar

# Inputs to change the file names, do not forget!
rmr = '1600'
path = 'rmr1600/rmr1600_test65'
name_output = 'test65_v4'

# To plot only 2 scalars and 2 eigevalues.
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

data = ascii.read(path + '.csv')


print('Starting iteration of ' + str(len(data['BINTIME'])) + ' loops')
progress = ProgressBar()
print(progress)

with PdfPages('plots/' + name_output + '.pdf') as pdf:

	for k in range(0, len(data['BINTIME'])):

		x_axis_1 = np.linspace(data['PERIOD_START'][k] - data['DELTA1'][k] * data['ITER1'][k] / 2, \
			data['PERIOD_START'][k] + data['DELTA1'][k] * (data['ITER1'][k]) / 2, data['ITER1'][k], endpoint=False)
		x_axis_2 = np.linspace(data['PERIOD_ITER1'][k] - data['DELTA2'][k] * data['ITER2'][k] / 2, \
			data['PERIOD_ITER1'][k] + data['DELTA2'][k] * (data['ITER2'][k]) / 2, data['ITER2'][k], endpoint=False)

		fig = plt.figure(figsize=(20,10))
		fig.subplots_adjust(hspace=.26, left=0.04, bottom=0.07, right=0.97, top=0.94, wspace=0.2)

		ax1 = fig.add_subplot(211)
		ax1.plot(x_axis_1, S0_i1[k], 'g+-', linewidth=1, label='First scalar')
		ax1.plot(x_axis_1, S1_i1[k], 'g+-', linewidth=0.8, label='Second scalar', alpha=0.5)
		ax1.plot(x_axis_1, V0_i1[k], 'b+-', linewidth=1, label='First eigenvalue')
		ax1.plot(x_axis_1, V1_i1[k], 'b+-', linewidth=0.8, label='Second eigenvalue', alpha=0.5)
		ax1.plot(x_axis_1, mstev1[k], 'r.-', linewidth=2.5, label='Merit function 1')
		ax1.axvline(x=0.08936715, color='y', linewidth=1, label='Folding Period = 0.08936715')
		ax1.axvline(x=data['PERIOD_ITER1'][k], color='m', linewidth=1, label='Best Period = ' + \
			str(data['PERIOD_ITER1'][k]))
		ax1.grid()
		ax1.legend(loc='upper right', fontsize=10)
		ax1.set_xlabel('Time [s]. Iterations = ' + str(data['ITER1'][k]) + ', delta = ' + str(data['DELTA1'][k]))
		ax1.set_ylabel('Scalar, eigvalue and merit amplitude')
		ax1.set_title('Noise ' + rmr + '. Iteration 1. Waterfall rows = ' + \
			str(data['NUM_DIV'][k]) + '. dt = ' + str(data['BINTIME'][k]) + ' s.')
		ax1.set_xlim([data['PERIOD_START'][k] - data['DELTA1'][k] * (data['ITER1'][k] + 1) / 2, \
			data['PERIOD_START'][k] + data['DELTA1'][k] * (data['ITER1'][k] - 1)/2])

		ax2 = fig.add_subplot(212)
		ax2.plot(x_axis_2, S0_i2[k], 'g+-', linewidth=1, label='First scalar')
		ax2.plot(x_axis_2, S1_i2[k], 'g+-', linewidth=0.8, label='Second scalar', alpha=0.5)
		ax2.plot(x_axis_2, V0_i2[k], 'b+-', linewidth=1, label='First eigenvalue')
		ax2.plot(x_axis_2, V1_i2[k], 'b+-', linewidth=0.8, label='Second eigenvalue', alpha=0.5)
		ax2.plot(x_axis_2, mstev2[k], 'r.-', linewidth=2.5, label='Merit function 2')
		ax2.plot(x_axis_1, mstev1[k], 'r.-', linewidth=1, label='Merit function 1')

		ax2.plot(gauss_fit(x_axis_2, mstev2[k])[0], gauss_fit(x_axis_2, mstev2[k])[1], 'dimgrey', linewidth=3)
		ax2.axvline(x=gauss_fit(x_axis_2, mstev2[k])[0][np.argmax(gauss_fit(x_axis_2, mstev2[k])[1])], \
			color='coral', linewidth=1, label='Gauss Period = ' + \
			str(gauss_fit(x_axis_2, mstev2[k])[0][np.argmax(gauss_fit(x_axis_2, mstev2[k])[1])]))

		ax2.axvline(x=0.08936715, color='y', linewidth=1, label='Folding Period = 0.08936715')
		ax2.axvline(x=data['PERIOD_ITER2'][k], color='m', linewidth=1, label='Best Period = ' + \
			str(data['PERIOD_ITER2'][k]))
		ax2.grid()
		ax2.legend(loc='upper left', fontsize=10)
		ax2.set_xlabel('Time [s]. Iterations = ' + str(data['ITER2'][k]) + ', delta = ' + str(data['DELTA2'][k]))
		ax2.set_ylabel('Scalar, eigvalue and merit amplitude')
		ax2.set_title('Noise ' + rmr + '. Iteration 2. Waterfall rows = ' + \
			str(data['NUM_DIV'][k]) + '. dt = ' + str(data['BINTIME'][k]) + ' s.')
		ax2.set_xlim([data['PERIOD_ITER1'][k] - data['DELTA2'][k] * (data['ITER2'][k] + 1) / 2, \
			data['PERIOD_ITER1'][k] + data['DELTA2'][k] * (data['ITER2'][k] - 1)/2])

		pdf.savefig()  # saves the current figure into a pdf page

		print(progress + 10 / len(data['BINTIME']))

	plt.close()

	d = pdf.infodict()
	d['Title'] = 'PCA Analysis'
	d['Author'] = u'Tomas Cassanelli'
	d['Keywords'] = 'PCA'
	d['ModDate'] = datetime.datetime.today()

print('It is done!')

# To open the already stored file
os.system('open plots/' + name_output + '.pdf')
