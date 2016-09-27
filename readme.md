# PCA Folding Analysis.
### Author: Tomas Cassanelli
### Staring date: 16-Sep-2016

## SUMMARY
The purpose of this program is to see the effectivity and efficiency that a PCA folding analysis has compared to the **Epoch Folding**, which is the standard used for pulsar period computations.

## PROGRAM:
The program is subdivided in two scripts, the first one, `pca_analysis.py`, it has all the functions and calculations and the second one, `pca_run.py`, runs the program and stores the generated matrices or plots in the working directory. Another extention is meant to be build which will be a GUI (soon).

### INPUTS:
1. CSV file which contains the time array that comes from the instrument. It is the time array that comes straight from the astronomical instrument when a pulsar is observed.
2. The bintime selected, in the program is also called `dt`.
3. The `period_start`, or the given period to start the iteration. Depending on how noisy or clean is the time array this can be computed or approximated since the beginning using a FFT. If the signal is not clear it is better to start with a close value to the real period.
4. Number of iterations, `iter1`, `iter2`, there are two big iterations loops, the second one should be the one with more iteration steps.
5. `delta1`, `delta2`, are the amount of time added or substracted depending on the iterations. They are of the order of the ns.
6. The `num_div`, which is in fact the number of how the time array is going to be divided. Also corresponds to the number of rows in the waterfall diagram.

### METHODS: To obtain the maximum and best PCA selection or in other words obtain the best period, it is necessary a selection method. Between then there are several.

1. Chosee the maximum between the first V or eigenvalues.
2. Maximize the SNR, where SNR = V[0]/V[1], or may be different. Where V represent the variance or the principal conponents eigenvalues.
3. A variant from method 2, using a SNR = V[0]/(V[1] + V[2] + ...)
4. Scalar maximum hyperdiagonal, selecting only the maximum eigenvalue's eigenvector. The scalar value is definced as the eigenvectors times the hyperdimensional unitary vector, the maximum absolute value is selected.
5. Multiplication of the maximum scalar and the corresponding eigenvalue's eigenvector. This is donde using the indices in each iteration. Pay attention to the variables `max_idx_scalar` in `pca_analysis.py`

### OUTPUTS: The principal outputs of the program will be the second and first periods from the iterations. It is possible to obtain other values from the residual data.

1. `period_final1`, `period_final2`, correspond to the best period found using the method.
2. The first three eigenvalues, in a `.npy` output format.
3. The first three scalars (hyperdiagonal unitary vector times the eigenvector), in a `.npy` output format.
3. The defined merit functions, exaplained in METHODS 5. in a `.csv` format. Please read this value using the package `import astropy.io as ascii`, then `data = ascii.read('file.csv')`