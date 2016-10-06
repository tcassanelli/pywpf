from astropy.io import ascii
import datetime
import os
import numpy as np
from cycler import cycler
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from matplotlib.backends.backend_pdf import PdfPages

# Importing data from global values
path = 'glob_var/'
data62 = ascii.read(path + 'gv_rmr0_test62.csv')
data63 = ascii.read(path + 'gv_rmr100_test63.csv')
data64 = ascii.read(path + 'gv_rmr500_test64.csv')
data65 = ascii.read(path + 'gv_rmr1600_test65.csv')
data66 = ascii.read(path + 'gv_rmr2000_test66.csv')
data67 = ascii.read(path + 'gv_rmr2600_test67.csv')
data68 = ascii.read(path + 'gv_rmr3000_test68.csv')

DATA = [data62, data63, data64, data65, data66, data67, data68]
M = len(DATA)

# Half value, to separate between two iterations
N = int(len(DATA[0]) / 2)
cm = plt.get_cmap('hot')
list_color = [cm(2.3 * i / N) for i in range(N)]

with PdfPages('plots/global_vals.pdf') as pdf:

	fig1 = plt.figure(figsize=(20,10))
	fig1.subplots_adjust(hspace=.26, left=0.04, bottom=0.07, right=0.97, top=0.94, wspace=0.13)

	ax1_1 = fig1.add_subplot(221)
	ax1_1.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax1_1.plot(DATA[i]['NUM_DIV'][:N], DATA[i]['MEAN'][:N], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax1_1.grid()
	ax1_1.set_title('First iteration mean value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax1_1.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax1_1.set_ylabel('MSTEV units')
	ax1_1.legend(loc='upper right', fontsize=10)

	ax1_2 = fig1.add_subplot(222)
	ax1_2.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax1_2.plot(DATA[i]['NUM_DIV'][N:], DATA[i]['MEAN'][N:], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax1_2.grid()
	ax1_2.set_title('Second iteration mean value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax1_2.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax1_2.set_ylabel('MSTEV units')
	ax1_2.legend(loc='lower left', fontsize=10)

	ax1_3 = fig1.add_subplot(223)
	ax1_3.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax1_3.plot(DATA[i]['NUM_DIV'][:N], DATA[i]['RMS'][:N], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax1_3.grid()
	ax1_3.set_title('First iteration RMS value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax1_3.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax1_3.set_ylabel('MSTEV units')
	ax1_3.legend(loc='upper left', fontsize=10)

	ax1_4 = fig1.add_subplot(224)
	ax1_4.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax1_4.plot(DATA[i]['NUM_DIV'][N:], DATA[i]['RMS'][N:], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax1_4.grid()
	ax1_4.set_title('Second iteration RMS value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax1_4.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax1_4.set_ylabel('MSTEV units')
	ax1_4.legend(loc='upper left', fontsize=10)

	pdf.savefig()

	fig2 = plt.figure(figsize=(20,10))
	fig2.subplots_adjust(hspace=.26, left=0.04, bottom=0.07, right=0.97, top=0.94, wspace=0.13)

	ax2_1 = fig2.add_subplot(221)
	ax2_1.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax2_1.plot(DATA[i]['NUM_DIV'][:N], DATA[i]['PEAK_SELECTED'][:N], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax2_1.grid()
	ax2_1.set_title('First iteration peak selected value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax2_1.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax2_1.set_ylabel('MSTEV units')
	ax2_1.legend(loc='upper left', fontsize=10)

	ax2_2 = fig2.add_subplot(222)
	ax2_2.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax2_2.plot(DATA[i]['NUM_DIV'][N:], DATA[i]['PEAK_SELECTED'][N:], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax2_2.grid()
	ax2_2.set_title('Second iteration peak selected value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax2_2.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax2_2.set_ylabel('MSTEV units')
	ax2_2.legend(loc='upper left', fontsize=10)

	ax2_3 = fig2.add_subplot(223)
	ax2_3.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax2_3.plot(DATA[i]['NUM_DIV'][:N], DATA[i]['PEAK_SIGNIFICANT'][:N], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax2_3.grid()
	ax2_3.set_title('First iteration peak significant value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax2_3.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax2_3.set_ylabel('MSTEV units')
	ax2_3.legend(loc='upper left', fontsize=10)

	ax2_4 = fig2.add_subplot(224)
	ax2_4.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax2_4.plot(DATA[i]['NUM_DIV'][N:], DATA[i]['PEAK_SIGNIFICANT'][N:], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax2_4.grid()
	ax2_4.set_title('Second iteration peak significant value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax2_4.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax2_4.set_ylabel('MSTEV units')
	ax2_4.legend(loc='upper left', fontsize=10)

	pdf.savefig()

	fig3 = plt.figure(figsize=(20,10))
	fig3.subplots_adjust(hspace=.26, left=0.04, bottom=0.07, right=0.97, top=0.94, wspace=0.13)

	ax3_1 = fig3.add_subplot(221)
	ax3_1.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax3_1.plot(DATA[i]['NUM_DIV'][:N], DATA[i]['KURTOSIS'][:N], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax3_1.grid()
	ax3_1.set_title('First iteration kurtosis value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax3_1.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax3_1.set_ylabel('MSTEV units')
	ax3_1.legend(loc='upper left', fontsize=10)

	ax3_2 = fig3.add_subplot(222)
	ax3_2.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax3_2.plot(DATA[i]['NUM_DIV'][N:], DATA[i]['KURTOSIS'][N:], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax3_2.grid()
	ax3_2.set_title('Second iteration kurtosis value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax3_2.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax3_2.set_ylabel('MSTEV units')
	ax3_2.legend(loc='upper left', fontsize=10)

	ax3_3 = fig3.add_subplot(223)
	ax3_3.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax3_3.plot(DATA[i]['NUM_DIV'][:N], DATA[i]['SKEWNESS'][:N], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax3_3.grid()
	ax3_3.set_title('First iteration skewness value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax3_3.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax3_3.set_ylabel('MSTEV units')
	ax3_3.legend(loc='upper left', fontsize=10)

	ax3_4 = fig3.add_subplot(224)
	ax3_4.set_prop_cycle(cycler('color', list_color))
	for i in range(M):
		ax3_4.plot(DATA[i]['NUM_DIV'][N:], DATA[i]['SKEWNESS'][N:], '.-', linewidth=2, label=DATA[i]['FILE'][0])
	ax3_4.grid()
	ax3_4.set_title('Second iteration skewness value. dt = ' + str(DATA[0]['BINTIME'][0]))
	ax3_4.set_xlabel('Number of eigenvalues (rows waterfall)')
	ax3_4.set_ylabel('MSTEV units')
	ax3_4.legend(loc='upper left', fontsize=10)

	pdf.savefig()

	plt.close()

	d = pdf.infodict()
	d['Title'] = 'PCA Analysis'
	d['Author'] = u'Tomas Cassanelli'
	d['Keywords'] = 'Global Values'
	d['ModDate'] = datetime.datetime.today()


