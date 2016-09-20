PCA Folding Analysis.
Author: Tomas Cassanelli
Date: 16-Sep-2016

Summary:
The purpose of this program is to see the effectivity and efficiency that a PCA folding analysis has compared to the Epoch folding, which is the standard used for pulsar period computations.

Program: 
The program is subdivided in two scripts, the first one, pca_analysis.py, has all the functions and calculations and the second one, pca_run.py, runs the program and stores the generated matrices or plots in the working directory. Another extention is meant to be build which will have a GUI.

What is wanted to be tested is the following

INPUTS: This are the files or numbers introduce to the program to start the calculations, between then thera are two critical values that we are interested in.

1. Matlab .mat file, which contains the time array, that comes from the instrument, i.e. rmr0.mat.
2. hte bintime selected, in the program is also called dt.
3. The period_start, or the given period to start the iteration. Depending on how noisy or clean is the time array this can be computed or approximated since the beginning using a FFT. If the signal is not clear it is better to start with a close value to the real period.
4. Number of iterations, iter1, iter2.
5. delta1, delta2, which are the amount of time added or substracted depending on the iterations. They are of the order of the ns.
6. The num_div, which is in fact the number of how the time array is going to be divided. Also corresponds to the number of rows in the waterfall diagram.

METHODS: To obtain the maximum and best PCA selection or in other words obtain the best period, it is necessary a selection method. Between then there are several.

1. Chosee the maximum between the first V or eigenvalues.
2. Maximize the SNR, where SNR = V[0]/V[1], or may be different.

OUTPUTS: The principal outputs of the program will be the second and first periods from the iterations. It is possible to obtain other values from the residual data,

1. period_final1, period_final2
2. Function dependent on the eigenvalues, organized by the method selection (x num_div), norm and signal matrices for the optimum selected, again depends on the selection method. 
3. The chi-squared plot -> you need to find out how to make it!!

How data will be stored (after computations):
Folder: Method1
Files: numpy.save -> .npy file for egigenvalues (method selection), norm and signal.
ascii -> for INPUT data and final period.